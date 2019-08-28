/*!
* @file graph.cpp
*
* @brief Graph class source file
*/

#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "thread_pool/thread_pool.hpp"
#include "ram/ram.hpp"

#include "pile.hpp"
#include "graph.hpp"

namespace raven {

std::unique_ptr<Graph> createGraph(std::shared_ptr<thread_pool::ThreadPool> thread_pool) {

    if (thread_pool == nullptr) {
        throw std::invalid_argument("[raven::createGraph] error: "
            "thread_pool is nullptr!");
    }

    return std::unique_ptr<Graph>(new Graph(thread_pool));
}

Graph::Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool)
        : thread_pool_(thread_pool), piles_() {
}

Graph::~Graph() {
}

void Graph::construct(std::vector<std::unique_ptr<ram::Sequence>>& sequences) {

    if (sequences.empty()) {
        return;
    }

    // ram::Overlap helper functions
    auto overlap_reverse = [] (const ram::Overlap& o) -> ram::Overlap {
        return ram::Overlap(o.t_id, o.t_begin, o.t_end, o.q_id, o.q_begin, o.q_end,
            o.strand, o.matches);
    };
    auto overlap_length = [] (const ram::Overlap& o) -> std::uint32_t {
        return std::max(o.t_end - o.t_begin, o.q_end - o.q_begin);
    };
    auto overlap_update = [&] (ram::Overlap& o) -> bool {
        if (piles_[o.q_id]->is_invalid() || piles_[o.t_id]->is_invalid()) {
            return false;
        }
        if (o.q_begin >= piles_[o.q_id]->end() || o.q_end <= piles_[o.q_id]->begin() ||
            o.t_begin >= piles_[o.t_id]->end() || o.t_end <= piles_[o.t_id]->begin()) {
            return false;
        }

        std::uint32_t q_begin = o.q_begin + (o.strand ?
            (o.t_begin < piles_[o.t_id]->begin() ? piles_[o.t_id]->begin() - o.t_begin : 0) :
            (o.t_end > piles_[o.t_id]->end() ? o.t_end - piles_[o.t_id]->end() : 0));
        std::uint32_t q_end = o.q_end - (o.strand ?
            (o.t_end > piles_[o.t_id]->end() ? o.t_end - piles_[o.t_id]->end() : 0) :
            (o.t_begin < piles_[o.t_id]->begin() ? piles_[o.t_id]->begin() - o.t_begin : 0));

        std::uint32_t t_begin = o.t_begin + (o.strand ?
            (o.q_begin < piles_[o.q_id]->begin() ? piles_[o.q_id]->begin() - o.q_begin : 0) :
            (o.q_end > piles_[o.q_id]->end() ? o.q_end - piles_[o.q_id]->end() : 0));
        std::uint32_t t_end = o.t_end - (o.strand ?
            (o.q_end > piles_[o.q_id]->end() ? o.q_end - piles_[o.q_id]->end() : 0) :
            (o.q_begin < piles_[o.q_id]->begin() ? piles_[o.q_id]->begin() - o.q_begin : 0));

        if (q_begin >= piles_[o.q_id]->end() || q_end <= piles_[o.q_id]->begin() ||
            t_begin >= piles_[o.t_id]->end() || t_end <= piles_[o.t_id]->begin()) {
            return false;
        }

        q_begin = std::max(q_begin, piles_[o.q_id]->begin());
        q_end = std::min(q_end, piles_[o.q_id]->end());
        t_begin = std::max(t_begin, piles_[o.t_id]->begin());
        t_end = std::min(t_end, piles_[o.t_id]->end());

        if (q_begin >= q_end || q_end - q_begin < 84 ||
            t_begin >= t_end || t_end - t_begin < 84) {
            return false;
        }

        o.q_begin = q_begin;
        o.q_end = q_end;
        o.t_begin = t_begin;
        o.t_end = t_end;

        return true;
    };
    auto overlap_type = [&] (const ram::Overlap& o) -> std::uint32_t {
        std::uint32_t q_length = piles_[o.q_id]->end() - piles_[o.q_id]->begin();
        std::uint32_t q_begin = o.q_begin - piles_[o.q_id]->begin();
        std::uint32_t q_end = o.q_end - piles_[o.q_id]->begin();

        std::uint32_t t_length = piles_[o.t_id]->end() - piles_[o.t_id]->begin();
        std::uint32_t t_begin = o.strand ?
            o.t_begin - piles_[o.t_id]->begin() :
            t_length - (o.t_end - piles_[o.t_id]->begin());
        std::uint32_t t_end = o.strand ?
            o.t_end - piles_[o.t_id]->begin():
            t_length - (o.t_begin - piles_[o.t_id]->begin());

        std::uint32_t overhang = std::min(q_begin, t_begin) +
            std::min(q_length - q_end, t_length - t_end);

        if (q_end - q_begin < (q_end - q_begin + overhang) * 0.875 ||
            t_end - t_begin < (t_end - t_begin + overhang) * 0.875) {
            return 0; // internal
        }
        if (q_begin <= t_begin && (q_length - q_end) <= (t_length - t_end)) {
            return 1; // q contained
        }
        if (t_begin <= q_begin && (t_length - t_end) <= (q_length - q_end)) {
            return 2; // t contained
        }
        if (q_begin > t_begin) {
            return 3; // q -> t
        }
        return 4; // t -> q
    };

    std::vector<std::vector<ram::Overlap>> overlaps(sequences.size());

    // connected components on overlaps
    auto connected_components = [&] () -> std::vector<std::vector<std::uint32_t>> {
        std::vector<std::vector<std::uint32_t>> connections(sequences.size());
        for (const auto& it: overlaps) {
            for (const auto& jt: it) {
                if (overlap_type(jt) > 2) {
                    connections[jt.q_id].emplace_back(jt.t_id);
                    connections[jt.t_id].emplace_back(jt.q_id);
                }
            }
        }

        std::vector<std::vector<std::uint32_t>> dst;
        std::vector<char> is_visited(sequences.size(), false);
        for (std::uint32_t i = 0; i < connections.size(); ++i) {
            if (piles_[i]->is_invalid() || is_visited[i]) {
                continue;
            }

            dst.resize(dst.size() + 1);

            std::deque<std::uint32_t> que = { i };
            while (!que.empty()) {
                std::uint32_t j = que.front();
                que.pop_front();

                if (is_visited[j]) {
                    continue;
                }
                is_visited[j] = true;
                dst.back().emplace_back(j);

                for (const auto& it: connections[j]) {
                    que.emplace_back(it);
                }
            }
        }

        return dst;
    };

    // find overlaps and create piles
    for (const auto& it: sequences) {
        piles_.emplace_back(createPile(it->id, it->data.size()));
    }

    auto minimizer_engine = ram::createMinimizerEngine(15, 5, thread_pool_);

    for (std::uint32_t i = 0, j = 0, bytes = 0; i < sequences.size(); ++i) {
        bytes += sequences[i]->data.size();
        if (i != sequences.size() - 1 && bytes < (1U << 30)) {
            continue;
        }

        minimizer_engine->minimize(sequences.begin() + j, sequences.begin() + i + 1);
        minimizer_engine->filter(0.001);

        std::vector<std::uint32_t> num_overlaps(overlaps.size());
        for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
            num_overlaps[k] = overlaps[k].size();
        }

        {
            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t k = 0; k < i + 1; ++k) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                        return minimizer_engine->map(sequences[i], true, true);
                    }
                , k));
            }
            for (auto& it: thread_futures) {
                it.wait();
                auto overlaps_part = it.get();
                for (const auto& jt: overlaps_part) {
                    overlaps[jt.q_id].emplace_back(jt);
                    overlaps[jt.t_id].emplace_back(overlap_reverse(jt));
                }
            }
        }

        {
            std::vector<std::future<void>> thread_futures;
            for (std::uint32_t k = 0; k < piles_.size(); ++k) {
                if (overlaps[k].empty()) {
                    continue;
                }
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> void {
                        piles_[i]->add_layers(overlaps[i].begin() + num_overlaps[i],
                            overlaps[i].end());

                        if (overlaps.size() < 16) {
                            return;
                        }

                        std::sort(overlaps[i].begin(), overlaps[i].end(),
                            [&] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
                                return overlap_length(lhs) > overlap_length(rhs);
                            }
                        );
                        overlaps[i].resize(16, ram::Overlap(0, 0, 0, 0, 0, 0, 0, 0));
                    }
                , k));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
        }

        bytes = 0;
        j = i + 1;
    }

    // trim and annotate piles
    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint32_t i) -> void {
                piles_[i]->find_valid_region(4);
                if (piles_[i]->is_invalid()) {
                    std::vector<ram::Overlap>().swap(overlaps[i]);
                } else {
                    piles_[i]->find_median();
                    piles_[i]->find_chimeric_regions();
                }
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    // find contained reads and resolve chimeric reads
    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        std::uint32_t k = 0;
        for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
            if (!overlap_update(overlaps[i][j])) {
                continue;
            }
            switch (overlap_type(overlaps[i][j])) {
                case 1:
                    if (!piles_[overlaps[i][j].t_id]->has_chimeric_region()) {
                        piles_[i]->set_contained();
                    }
                    break;
                case 2:
                    if (!piles_[i]->has_chimeric_region()) {
                        piles_[overlaps[i][j].t_id]->set_contained();
                    }
                    break;
                default: overlaps[i][k++] = overlaps[i][j]; break;
            }
        }
        overlaps[i].resize(k, ram::Overlap(0, 0, 0, 0, 0, 0, 0, 0));
    }
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
        if (piles_[i]->is_contained()) {
            piles_[i]->set_invalid();
            std::vector<ram::Overlap>().swap(overlaps[i]);
        }
    }

    while (true) {
        auto components = connected_components();
        for (const auto& it: components) {

            std::vector<std::uint32_t> medians;
            for (const auto& jt: it) {
                medians.emplace_back(piles_[jt]->median());
            }
            std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                medians.end());
            std::uint32_t median = medians[medians.size() / 2];

            for (const auto& jt: it) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> void {
                        piles_[i]->resolve_chimeric_regions(median);
                        if (piles_[i]->is_invalid()) {
                            std::vector<ram::Overlap>().swap(overlaps[i]);
                        }
                    }
                , jt));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
            thread_futures.clear();
        }

        bool is_changed = false;
        for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
            std::uint32_t k = 0;
            for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
                if (overlap_update(overlaps[i][j])) {
                    overlaps[i][k++] = overlaps[i][j];
                } else {
                    is_changed = true;
                }
            }
            overlaps[i].resize(k, ram::Overlap(0, 0, 0, 0, 0, 0, 0, 0));
        }

        if (!is_changed) {
            for (const auto& it: overlaps) {
                for (const auto& jt: it) {
                    switch (overlap_type(jt)) {
                        case 1:
                            piles_[jt.q_id]->set_contained();
                            piles_[jt.q_id]->set_invalid();
                            break;
                        case 2:
                            piles_[jt.t_id]->set_contained();
                            piles_[jt.t_id]->set_invalid();
                            break;
                        default: break;
                    }
                }
            }
            overlaps.clear();
            break;
        }
    }

    // update piles with more sensitive overlaps
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
        if (piles_[i]->is_invalid()) {
            continue;
        }
        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint32_t i) -> void {
                piles_[i]->clear_valid_region();
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<ram::Sequence>& lhs,
            const std::unique_ptr<ram::Sequence>& rhs) -> bool {
            return piles_[lhs->id]->is_invalid() < piles_[rhs->id]->is_invalid();
        }
    );
    std::uint32_t s = 0;
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        if (piles_[sequences[i]->id]->is_invalid()) {
            s = i;
            break;
        }
    }

    overlaps.emplace_back(std::vector<ram::Overlap>()); // valid x valid
    overlaps.emplace_back(std::vector<ram::Overlap>()); // invalid x valid
    for (std::uint32_t i = 0, j = 0, bytes = 0; i < s; ++i) {
        bytes += sequences[i]->data.size();
        if (i != s - 1 && bytes < (1U << 30)) {
            continue;
        }

        minimizer_engine->minimize(sequences.begin() + j, sequences.begin() + i + 1);

        { // map valid reads to each other
            minimizer_engine->filter(0.001);
            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t k = 0; k < s; ++k) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                        return minimizer_engine->map(sequences[i], true, true);
                    }
                , k));
            }
            for (auto& it: thread_futures) {
                it.wait();
                auto overlaps_part = it.get();
                for (auto& jt: overlaps_part) {
                    if (!overlap_update(jt) || overlap_type(jt) < 3) {
                        continue;
                    }
                    if (overlaps.front().size() &&
                        overlaps.front().back().q_id == jt.q_id &&
                        overlaps.front().back().t_id == jt.t_id) {
                        if (overlap_length(overlaps.front().back()) < overlap_length(jt)) {
                            overlaps.front().back() = jt;
                        }
                    } else {
                        overlaps.front().emplace_back(jt);
                    }
                }
            }
        }

        { // map invalid reads to valid reads
            minimizer_engine->filter(0.00001);
            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t k = s; k < sequences.size(); ++k) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                        return minimizer_engine->map(sequences[i], true, false);
                    }
                , k));
            }
            for (auto& it: thread_futures) {
                it.wait();
                auto overlaps_part = it.get();
                overlaps.back().insert(overlaps.back().end(),
                    overlaps_part.begin(), overlaps_part.end());
            }
        }

        {
            std::vector<std::future<void>> thread_futures;
            for (std::uint32_t k = 0; k < s; ++k) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> void {
                        piles_[i]->add_layers(overlaps.back().begin(),
                            overlaps.back().end());
                    }
                , sequences[k]->id));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
        }
        overlaps.back().clear();

        bytes = 0;
        j = i + 1;
    }

    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
        if (piles_[i]->is_invalid()) {
            continue;
        }
        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint32_t i) -> void {
                piles_[i]->clear_invalid_region();
                piles_[i]->find_median();
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    std::cerr << overlaps.front().size() << std::endl;

    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<ram::Sequence>& lhs,
            const std::unique_ptr<ram::Sequence>& rhs) -> bool {
            return lhs->id < rhs->id;
        }
    );

    // annotate and resolve repetitive regions
    while (true) {
        auto components = connected_components();
        for (const auto& it: components) {

            std::vector<std::uint32_t> medians;
            for (const auto& jt: it) {
                medians.emplace_back(piles_[jt]->median());
            }
            std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                medians.end());
            std::uint32_t median = medians[medians.size() / 2];

            for (const auto& jt: it) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> void {
                        piles_[i]->find_repetitive_regions(median);
                    }
                , jt));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
            thread_futures.clear();
        }

        for (const auto& it: overlaps.front()) {
            piles_[it.q_id]->resolve_repetitive_regions(it);
            piles_[it.t_id]->resolve_repetitive_regions(it);
        }

        bool is_changed = false;
        std::uint32_t j = 0;
        for (std::uint32_t i = 0; i < overlaps.front().size(); ++i) {
            const auto& it = overlaps.front()[i];
            if (piles_[it.q_id]->is_false_overlap(it) ||
                piles_[it.t_id]->is_false_overlap(it)) {
                is_changed = true;
            } else {
                overlaps.front()[j++] = it;
            }
        }
        overlaps.front().resize(j, ram::Overlap(0, 0, 0, 0, 0, 0, 0, 0));

        if (!is_changed) {
            break;
        }

        for (const auto& it: components) {
            for (const auto& jt: it) {
                piles_[jt]->clear_repetitive_regions();
            }
        }
    }

    // construct assembly graph
}

void Graph::print_json(const std::string& path) const {

    std::ofstream os(path);
    os << "{\"piles\":{";
    bool is_first = true;

    for (const auto& it: piles_) {
        if (it->is_invalid() || !it->has_repetitive_region()) {
            continue;
        }
        if (!is_first) {
            os << ",";
        }
        is_first = false;
        os << it->to_json();
    }

    os << "}}";
    os.close();
}

}
