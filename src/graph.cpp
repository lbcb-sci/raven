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

void Graph::construct(const std::vector<std::unique_ptr<ram::Sequence>>& sequences) {

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
        if (!piles_[o.q_id]->is_valid() || !piles_[o.t_id]->is_valid()) {
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

    // find overlaps and create piles
    for (const auto& it: sequences) {
        piles_.emplace_back(createPile(it->id, it->data.size()));
    }

    std::vector<std::vector<ram::Overlap>> overlaps(sequences.size());
    auto minimizer_engine = ram::createMinimizerEngine(15, 5, thread_pool_);

    for (std::uint32_t i = 0, j = 0, bytes = 0; i < sequences.size(); ++i) {
        bytes += sequences[i]->data.size();
        if (i != sequences.size() - 1 && bytes < (2U << 30)) {
            continue;
        }

        minimizer_engine->minimize(sequences.begin() + j, sequences.begin() + i + 1,
            0.001);

        std::vector<std::uint32_t> num_overlaps(overlaps.size());
        for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
            num_overlaps[k] = overlaps[k].size();
        }

        {
            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t k = 0; k < sequences.size(); ++k) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t id) -> std::vector<ram::Overlap> {
                        return minimizer_engine->map(sequences[id], true, true);
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
            for (std::uint32_t k = 0; k < sequences.size(); ++k) {
                if (overlaps[k].empty()) {
                    continue;
                }
                thread_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t id) -> void {
                        piles_[id]->add_layers(overlaps[id].begin() + num_overlaps[id],
                            overlaps[id].end());

                        if (overlaps.size() < 16) {
                            return;
                        }

                        std::sort(overlaps[id].begin(), overlaps[id].end(),
                            [&] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
                                return overlap_length(lhs) > overlap_length(rhs);
                            }
                        );
                        overlaps[id].resize(16, ram::Overlap(0, 0, 0, 0, 0, 0, 0, 0));
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
            [&] (std::uint32_t id) -> void {
                piles_[id]->find_valid_region(4);
                if (piles_[id]->is_valid()) {
                    piles_[id]->find_median();
                    piles_[id]->find_chimeric_regions();
                } else {
                    std::vector<ram::Overlap>().swap(overlaps[id]);
                }
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    // update overlaps
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit(
            [&] (std::uint32_t id) -> void {
                std::uint32_t k = 0;
                for (std::uint32_t j = 0; j < overlaps[id].size(); ++j) {
                    if (overlap_update(overlaps[id][j])) {
                        overlaps[id][k++] = overlaps[id][j];
                    }
                }
                overlaps[id].resize(k, ram::Overlap(0, 0, 0, 0, 0, 0, 0, 0));
            }
        , i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();
}

}
