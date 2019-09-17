/*!
* @file graph.cpp
*
* @brief Graph class source file
*/

#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"
#include "ram/ram.hpp"

#include "pile.hpp"
#include "graph.hpp"

namespace raven {

struct Graph::Node {
    Node(std::uint32_t id, std::uint32_t sequence, const std::string& name,
        const std::string& data);
    Node(std::uint32_t id, Node* begin, Node* end);
    Node(const Node&) = delete;
    const Node& operator=(const Node&) = delete;

    ~Node() = default;

    bool is_rc() const {
        return id & 1;
    }

    std::uint32_t length() const {
        return data.size();
    }

    std::uint32_t indegree() const {
        return inedges.size();
    }
    std::uint32_t outdegree() const {
        return outedges.size();
    }

    bool is_junction() const {
        return outdegree() > 1 || indegree() > 1;
    }
    bool is_tip() const {
        return outdegree() > 0 && indegree() == 0 && sequences.size() < 6;
    }

    std::uint32_t id;
    std::string name;
    std::string data;
    std::uint32_t state;
    std::vector<std::uint32_t> sequences;
    std::vector<Edge*> inedges;
    std::vector<Edge*> outedges;
    std::unordered_set<std::uint32_t> transitive;
    Node* pair;
};

struct Graph::Edge {
    Edge(std::uint32_t id, Node* begin, Node* end, std::uint32_t length);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    ~Edge() = default;

    std::string label() const {
        return begin->data.substr(0, length);
    }

    std::uint32_t id;
    std::uint32_t length;
    std::uint32_t state;
    double weight;
    Node* begin;
    Node* end;
    Edge* pair;
};

std::unique_ptr<Graph> createGraph(std::shared_ptr<thread_pool::ThreadPool> thread_pool) {

    if (thread_pool == nullptr) {
        throw std::invalid_argument("[raven::createGraph] error: "
            "thread_pool is nullptr!");
    }

    return std::unique_ptr<Graph>(new Graph(thread_pool));
}

Graph::Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool)
        : minimizer_engine_(ram::createMinimizerEngine(15, 5, thread_pool)),
        thread_pool_(thread_pool), piles_(), nodes_(), edges_(),
        marked_edges_() {
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

    logger::Logger logger;

    // find overlaps and create piles
    for (const auto& it: sequences) {
        piles_.emplace_back(createPile(it->id, it->data.size()));
    }

    for (std::uint32_t i = 0, j = 0, bytes = 0; i < sequences.size(); ++i) {
        bytes += sequences[i]->data.size();
        if (i != sequences.size() - 1 && bytes < (1U << 30)) {
            continue;
        }
        bytes = 0;

        logger.log();

        minimizer_engine_->minimize(sequences.begin() + j, sequences.begin() + i + 1);
        minimizer_engine_->filter(0.001);

        logger.log("[raven::Graph::construct] minimized " +
            std::to_string(j) + " - " + std::to_string(i + 1) +
            " / " + std::to_string(sequences.size()));
        logger.log();

        std::vector<std::uint32_t> num_overlaps(overlaps.size());
        for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
            num_overlaps[k] = overlaps[k].size();
        }

        std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;

        for (std::uint32_t k = 0; k < i + 1; ++k) {
            thread_futures.emplace_back(thread_pool_->submit(
                [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                    return minimizer_engine_->map(sequences[i], true, true);
                }
            , k));

            bytes += sequences[k]->data.size();
            if (k != i && bytes < (1U << 30)) {
                continue;
            }
            bytes = 0;

            for (auto& it: thread_futures) {
                for (const auto& jt: it.get()) {
                    overlaps[jt.q_id].emplace_back(jt);
                    overlaps[jt.t_id].emplace_back(overlap_reverse(jt));
                }
            }
            thread_futures.clear();

            std::vector<std::future<void>> void_futures;
            for (const auto& it: piles_) {
                if (overlaps[it->id()].empty()) {
                    continue;
                }
                void_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> void {
                        piles_[i]->add_layers(overlaps[i].begin() + num_overlaps[i],
                            overlaps[i].end());

                        if (overlaps[i].size() < 16) {
                            return;
                        }

                        std::sort(overlaps[i].begin(), overlaps[i].end(),
                            [&] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
                                return overlap_length(lhs) > overlap_length(rhs);
                            }
                        );

                        std::vector<ram::Overlap> o;
                        o.insert(o.end(), overlaps[i].begin(), overlaps[i].begin() + 16);
                        o.swap(overlaps[i]);
                    }
                , it->id()));
            }
            for (const auto& it: void_futures) {
                it.wait();
            }
        }

        logger.log("[raven::Graph::construct] mapped sequences");

        j = i + 1;
    }

    logger.log();

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

    logger.log("[raven::Graph::construct] annotated piles");
    logger.log();

    // find contained reads and resolve chimeric reads
    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        std::uint32_t k = 0;
        for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
            if (!overlap_update(overlaps[i][j])) {
                continue;
            }
            std::uint32_t type = overlap_type(overlaps[i][j]);
            if (type == 1 && !piles_[overlaps[i][j].t_id]->has_chimeric_region()) {
                piles_[i]->set_contained();
            } else if (type == 2 && !piles_[i]->has_chimeric_region()) {
                piles_[overlaps[i][j].t_id]->set_contained();
            } else {
                overlaps[i][k++] = overlaps[i][j];
            }
        }
        overlaps[i].resize(k);
    }
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
        if (piles_[i]->is_contained()) {
            piles_[i]->set_invalid();
            std::vector<ram::Overlap>().swap(overlaps[i]);
        }
    }

    logger.log("[raven::Graph::construct] removed contained sequences");
    logger.log();

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
            overlaps[i].resize(k);
        }

        if (!is_changed) {
            for (const auto& it: overlaps) {
                for (const auto& jt: it) {
                    std::uint32_t type = overlap_type(jt);
                    if (type == 1) {
                        piles_[jt.q_id]->set_contained();
                        piles_[jt.q_id]->set_invalid();
                    } else if (type == 2) {
                        piles_[jt.t_id]->set_contained();
                        piles_[jt.t_id]->set_invalid();
                    }
                }
            }
            overlaps.clear();
            break;
        }
    }

    logger.log("[raven::Graph::construct] removed chimeric sequences");
    logger.log();

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

    logger.log("[raven::Graph::construct] cleared piles");
    logger.log();

    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<ram::Sequence>& lhs,
            const std::unique_ptr<ram::Sequence>& rhs) -> bool {
            return piles_[lhs->id]->is_invalid() < piles_[rhs->id]->is_invalid() ||
                (piles_[lhs->id]->is_invalid() == piles_[rhs->id]->is_invalid() && lhs->id < rhs->id);
        }
    );
    std::uint32_t s = 0;
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        if (piles_[sequences[i]->id]->is_invalid()) {
            s = i;
            break;
        }
    }

    logger.log("[raven::Graph::construct] rearranged sequences");

    overlaps.resize(sequences.size() + 1);
    for (std::uint32_t i = 0, j = 0, bytes = 0; i < s; ++i) {
        bytes += sequences[i]->data.size();
        if (i != s - 1 && bytes < (1U << 30)) {
            continue;
        }
        bytes = 0;

        logger.log();

        minimizer_engine_->minimize(sequences.begin() + j, sequences.begin() + i + 1);

        logger.log("[raven::Graph::construct] minimized " +
            std::to_string(j) + " - " + std::to_string(i + 1) +
            " / " + std::to_string(s));
        logger.log();

        // map valid reads to each other
        std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
        minimizer_engine_->filter(0.001);
        for (std::uint32_t k = 0; k < i + 1; ++k) {
            thread_futures.emplace_back(thread_pool_->submit(
                [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                    return minimizer_engine_->map(sequences[i], true, true);
                }
            , k));
        }
        for (auto& it: thread_futures) {
            for (auto& jt: it.get()) {
                if (!overlap_update(jt)) {
                    continue;
                }
                std::uint32_t type = overlap_type(jt);
                if (type == 1) {
                    continue;
                } else if (type == 1) {
                    piles_[jt.q_id]->set_contained();
                } else if (type == 2) {
                    piles_[jt.t_id]->set_contained();
                } else {
                    if (overlaps.back().size() &&
                        overlaps.back().back().q_id == jt.q_id &&
                        overlaps.back().back().t_id == jt.t_id) {
                        if (overlap_length(overlaps.back().back()) < overlap_length(jt)) {
                            overlaps.back().back() = jt;
                        }
                    } else {
                        overlaps.back().emplace_back(jt);
                    }
                }
            }
        }
        thread_futures.clear();

        logger.log("[raven::Graph::construct] mapped valid sequences");
        logger.log();

        // map invalid reads to valid reads
        minimizer_engine_->filter(0.00001);
        for (std::uint32_t k = s; k < sequences.size(); ++k) {
            thread_futures.emplace_back(thread_pool_->submit(
                [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                    return minimizer_engine_->map(sequences[i], true, false);
                }
            , k));

            bytes += sequences[k]->data.size();
            if (k != sequences.size() - 1 && bytes < (1U << 30)) {
                continue;
            }
            bytes = 0;

            for (auto& it: thread_futures) {
                for (const auto& jt: it.get()) {
                    overlaps[jt.t_id].emplace_back(jt);
                }
            }
            thread_futures.clear();

            std::vector<std::future<void>> void_futures;
            for (std::uint32_t k = j; k < i + 1; ++k) {
                if (overlaps[sequences[k]->id].empty()) {
                    continue;
                }
                void_futures.emplace_back(thread_pool_->submit(
                    [&] (std::uint32_t i) -> void {
                        piles_[i]->add_layers(overlaps[i].begin(), overlaps[i].end());
                        std::vector<ram::Overlap>().swap(overlaps[i]);
                    }
                , sequences[k]->id));
            }
            for (const auto& it: void_futures) {
                it.wait();
            }
        }

        logger.log("[raven::Graph::construct] mapped invalid sequences");

        j = i + 1;
    }

    logger.log();

    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
        if (piles_[i]->is_contained()) {
            piles_[i]->set_invalid();
            continue;
        }
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

    logger.log("[raven::Graph::construct] updated piles");
    logger.log();

    {
        std::uint32_t k = 0;
        for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
            if (overlap_update(overlaps.back()[i])) {
                overlaps.back()[k++] = overlaps.back()[i];
            }
        }
        overlaps.back().resize(k);
    }

    logger.log("[raven::Graph::construct] updated overlaps");
    logger.log();

    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<ram::Sequence>& lhs,
            const std::unique_ptr<ram::Sequence>& rhs) -> bool {
            return lhs->id < rhs->id;
        }
    );

    logger.log("[raven::Graph::construct] rearranged sequences");
    logger.log();

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

        for (const auto& it: overlaps.back()) {
            piles_[it.q_id]->resolve_repetitive_regions(it);
            piles_[it.t_id]->resolve_repetitive_regions(it);
        }

        bool is_changed = false;
        std::uint32_t j = 0;
        for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
            const auto& it = overlaps.back()[i];
            if (piles_[it.q_id]->is_false_overlap(it) ||
                piles_[it.t_id]->is_false_overlap(it)) {
                is_changed = true;
            } else {
                overlaps.back()[j++] = it;
            }
        }
        overlaps.back().resize(j);

        if (!is_changed) {
            break;
        }

        for (const auto& it: components) {
            for (const auto& jt: it) {
                piles_[jt]->clear_repetitive_regions();
            }
        }
    }

    logger.log("[raven::Graph::construct] removed false overlaps");
    logger.log();

    // construct assembly graph
    auto reverse_complement = [] (const std::string& src) -> std::string {
        std::string dst;
        for (const auto& it: src) {
            switch (it) {
                case 'A': case 'a': dst += 'T'; break;
                case 'T': case 't': dst += 'A'; break;
                case 'G': case 'g': dst += 'C'; break;
                case 'C': case 'c': dst += 'G'; break;
                default: dst += it; break;
            }
        }
        std::reverse(dst.begin(), dst.end());
        return dst;
    };

    std::vector<std::int32_t> sequence_to_node(piles_.size(), -1);
    std::uint32_t node_id = 0;
    for (const auto& it: piles_) {
        if (it->is_invalid()) {
            continue;
        }

        sequence_to_node[it->id()] = node_id;

        const auto& name = sequences[it->id()]->name;
        const auto& data = sequences[it->id()]->data.substr(it->begin(),
            it->end() - it->begin());

        std::unique_ptr<Node> n(new Node(node_id++, it->id(), name, data));
        std::unique_ptr<Node> nc(new Node(node_id++, it->id(), name,
            reverse_complement(data)));

        n->pair = nc.get();
        nc->pair = n.get();

        nodes_.emplace_back(std::move(n));
        nodes_.emplace_back(std::move(nc));
    }

    logger.log("[raven::Graph::construct] stored nodes");
    logger.log();

    std::uint32_t edge_id = 0;
    for (const auto& it: overlaps.back()) {
        Node* q = nodes_[sequence_to_node[it.q_id]].get();
        Node* t = nodes_[sequence_to_node[it.t_id] + 1 - it.strand].get();

        std::uint32_t q_length = piles_[it.q_id]->end() - piles_[it.q_id]->begin();
        std::uint32_t q_begin = it.q_begin - piles_[it.q_id]->begin();
        std::uint32_t q_end = it.q_end - piles_[it.q_id]->begin();

        std::uint32_t t_length = piles_[it.t_id]->end() - piles_[it.t_id]->begin();
        std::uint32_t t_begin = it.strand ?
            it.t_begin - piles_[it.t_id]->begin() :
            t_length - (it.t_end - piles_[it.t_id]->begin());
        std::uint32_t t_end = it.strand ?
            it.t_end - piles_[it.t_id]->begin():
            t_length - (it.t_begin - piles_[it.t_id]->begin());

        std::uint32_t type = overlap_type(it);

        if (type == 3) {
            std::unique_ptr<Edge> e(new Edge(edge_id++, q, t, q_begin - t_begin));
            std::unique_ptr<Edge> ec(new Edge(edge_id++, t->pair, q->pair,
                (t_length - t_end) - (q_length - q_end)));

            e->pair = ec.get();
            ec->pair = e.get();

            q->outedges.emplace_back(e.get());
            q->pair->inedges.emplace_back(ec.get());
            t->inedges.emplace_back(e.get());
            t->pair->outedges.emplace_back(ec.get());

            edges_.emplace_back(std::move(e));
            edges_.emplace_back(std::move(ec));
        } else if (type == 4) {
            std::unique_ptr<Edge> e(new Edge(edge_id++, t, q, t_begin - q_begin));
            std::unique_ptr<Edge> ec(new Edge(edge_id++, q->pair, t->pair,
                (q_length - q_end) - (t_length - t_end)));

            e->pair = ec.get();
            ec->pair = e.get();

            t->outedges.emplace_back(e.get());
            t->pair->inedges.emplace_back(ec.get());
            q->inedges.emplace_back(e.get());
            q->pair->outedges.emplace_back(ec.get());

            edges_.emplace_back(std::move(e));
            edges_.emplace_back(std::move(ec));
        }
    }

    logger.log("[raven::Graph::construct] stored edges");
    logger.total("[raven::Graph::construct]");
}

void Graph::assemble(std::vector<std::unique_ptr<ram::Sequence>>& dst) {

    logger::Logger logger;
    logger.log();

    remove_transitive_edges();

    logger.log("[raven::Graph::assemble] removed transitive edges");
    logger.log();

    while (true) {
        std::uint32_t num_changes = remove_tips();
        num_changes += remove_bubbles();
        if (num_changes == 0) {
            break;
        }
    }

    logger.log("[raven::Graph::assemble] removed tips and bubbles");
    logger.log();

    create_unitigs(42);
    for (std::uint32_t i = 0; i < 5; ++i) {
        create_force_directed_layout();
        remove_long_edges();
        remove_tips();
    }

    logger.log("[raven::Graph::assemble] removed long edges");

    while (true) {
        std::uint32_t num_changes = remove_tips();
        num_changes += remove_bubbles();
        if (num_changes == 0) {
            break;
        }
    }

    extract_unitigs(dst);

    logger.total("[raven::Graph::assemble]");
}

std::uint32_t Graph::remove_transitive_edges() {

    std::uint32_t num_transitive = 0;
    std::vector<Edge*> candidate(nodes_.size(), nullptr);

    auto comparable = [] (double a, double b, double eps) -> bool {
        return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
               (b >= a * (1 - eps) && b <= a * (1 + eps));
    };

    for (const auto& node_a: nodes_) {
        if (node_a == nullptr) {
            continue;
        }

        for (const auto& edge_ab: node_a->outedges) {
            candidate[edge_ab->end->id] = edge_ab;
        }

        for (const auto& edge_ab: node_a->outedges) {
            const auto& node_b = nodes_[edge_ab->end->id];

            for (const auto& edge_bc: node_b->outedges) {
                uint64_t c = edge_bc->end->id;

                if (candidate[c] != nullptr && !(candidate[c]->state & 1)) {
                    if (comparable(edge_ab->length + edge_bc->length,
                        candidate[c]->length, 0.12)) {

                        candidate[c]->state |= 1;
                        candidate[c]->pair->state |= 1;
                        marked_edges_.emplace(candidate[c]->id);
                        marked_edges_.emplace(candidate[c]->pair->id);
                        ++num_transitive;
                    }
                }
            }
        }

        for (const auto& edge_ab: node_a->outedges) {
            candidate[edge_ab->end->id] = nullptr;
        }
    }

    for (const auto& it: marked_edges_) {
        if (it & 1) {
            nodes_[(edges_[it]->begin->id >> 1) << 1]->transitive.emplace(
                (edges_[it]->end->id >> 1) << 1);
            nodes_[(edges_[it]->end->id >> 1) << 1]->transitive.emplace(
                (edges_[it]->begin->id >> 1) << 1);
        }
    }

    remove_marked_objects();

    return num_transitive;
}

std::uint32_t Graph::remove_tips() {

    std::uint32_t num_tip_edges = 0;
    std::vector<char> is_visited(nodes_.size(), 0);

    for (const auto& it: nodes_) {
        if (it == nullptr || is_visited[it->id] || !it->is_tip()) {
            continue;
        }

        bool is_circular = false;
        std::uint32_t num_sequences = 0;

        auto end = it.get();
        while (!end->is_junction()) {
            num_sequences += end->sequences.size();
            is_visited[end->id] = 1;
            is_visited[end->pair->id] = 1;
            if (end->outdegree() == 0 ||
                end->outedges[0]->end->is_junction()) {
                break;
            }
            end = end->outedges[0]->end;
            if (end->id == it->id) {
                is_circular = true;
                break;
            }
        }

        if (is_circular || end->outdegree() == 0 || num_sequences > 5) {
            continue;
        }

        std::uint32_t num_removed_edges = 0;

        for (const auto& edge: end->outedges) {
            if (edge->end->indegree() > 1) {
                edge->state |= 1;
                edge->pair->state |= 1;
                marked_edges_.emplace(edge->id);
                marked_edges_.emplace(edge->pair->id);
                ++num_removed_edges;
            }
        }

        if (num_removed_edges == end->outedges.size()) {
            auto curr = it.get();
            while (curr->id != end->id) {
                curr->outedges[0]->state |= 1;
                curr->outedges[0]->pair->state |= 1;
                marked_edges_.emplace(curr->outedges[0]->id);
                marked_edges_.emplace(curr->outedges[0]->pair->id);
                curr = curr->outedges[0]->end;
            }
        }

        num_tip_edges += num_removed_edges;

        remove_marked_objects(true);
    }

    return num_tip_edges;
}

std::uint32_t Graph::remove_bubbles() {

    std::vector<std::uint32_t> distance(nodes_.size(), 0);
    std::vector<std::uint32_t> visited(nodes_.size(), 0);
    std::uint32_t visited_length = 0;
    std::vector<std::int32_t> predecessor(nodes_.size(), -1);
    std::deque<std::uint32_t> node_queue;

    // helper functions
    auto path_extract = [&] (std::vector<std::uint32_t>& dst,
        std::uint32_t source, std::uint32_t sink) -> void {

        std::uint32_t curr_id = sink;
        while (curr_id != source) {
            dst.emplace_back(curr_id);
            curr_id = predecessor[curr_id];
        }
        dst.emplace_back(source);
        std::reverse(dst.begin(), dst.end());
    };
    auto path_type = [&] (const std::vector<std::uint32_t>& path) -> bool {
        if (path.empty()) {
            return false;
        }
        for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
            if (nodes_[path[i]]->is_junction()) {
                return false;
            }
        }
        return true;
    };
    auto bubble_type = [&] (const std::vector<std::uint32_t>& path,
        const std::vector<std::uint32_t>& other_path) -> bool {

        if (path.empty() || other_path.empty()) {
            return false;
        }

        std::unordered_set<std::uint32_t> node_set;
        for (const auto& it: path) {
            node_set.emplace(it);
        }
        for (const auto& it: other_path) {
            node_set.emplace(it);
        }
        if (path.size() + other_path.size() - 2 != node_set.size()) {
            return false;
        }
        for (const auto& it: path) {
            std::uint32_t pair_id = nodes_[it]->pair->id;
            if (node_set.count(pair_id) != 0) {
                return false;
            }
        }

        if (path_type(path) && path_type(other_path)) {
            return true;
        }

        auto sequence = std::unique_ptr<ram::Sequence>(new ram::Sequence());
        auto other_sequence = std::unique_ptr<ram::Sequence>(new ram::Sequence());

        for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
            for (const auto& it: nodes_[path[i]]->outedges) {
                if (it->end->id == path[i + 1]) {
                    sequence->data += it->label();
                    break;
                }
            }
        }
        sequence->data += nodes_[path.back()]->data;

        for (std::uint32_t i = 0; i < other_path.size() - 1; ++i) {
            for (const auto& it: nodes_[other_path[i]]->outedges) {
                if (it->end->id == other_path[i + 1]) {
                    other_sequence->data += it->label();
                    break;
                }
            }
        }
        other_sequence->data += nodes_[other_path.back()]->data;

        if (std::min(sequence->data.size(), other_sequence->data.size()) <
            std::max(sequence->data.size(), other_sequence->data.size()) * 0.8) {
            return false;
        }

        auto overlaps = minimizer_engine_->map(sequence, other_sequence);

        std::uint32_t matches = 0, length = 0;
        for (const auto& it: overlaps) {
            std::uint32_t l = std::max(it.q_end - it.q_begin, it.t_end - it.t_begin);
            if (length < l) {
                length = l;
                matches = it.matches;
            }
        }

        return static_cast<double>(matches) / length > 0.5;
    };

    std::uint32_t num_bubbles_popped = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->outdegree() < 2) {
            continue;
        }

        bool found_sink = false;
        std::uint32_t sink = 0, sink_other_predecesor = 0;
        std::uint32_t source = node->id;

        // BFS
        node_queue.emplace_back(source);
        visited[visited_length++] = source;
        while (!node_queue.empty() && !found_sink) {
            std::uint32_t v = node_queue.front();
            const auto& curr_node = nodes_[v];

            node_queue.pop_front();

            for (const auto& edge: curr_node->outedges) {
                std::uint32_t w = edge->end->id;

                if (w == source) {
                    continue; // Cycle
                }

                if (distance[v] + edge->length > 500000) {
                    continue; // Out of reach
                }

                distance[w] = distance[v] + edge->length;
                visited[visited_length++] = w;
                    node_queue.emplace_back(w);

                if (predecessor[w] != -1) {
                    sink = w;
                    sink_other_predecesor = v;
                    found_sink = true;
                    break;
                }

                predecessor[w] = v;
            }
        }

        if (found_sink) {
            std::vector<std::uint32_t> path;
            path_extract(path, source, sink);

            std::vector<std::uint32_t> other_path(1, sink);
            path_extract(other_path, source, sink_other_predecesor);

            if (bubble_type(path, other_path)) {
                std::uint32_t path_num_reads = 0;
                for (const auto& it: path) {
                    path_num_reads += nodes_[it]->sequences.size();
                }

                std::uint32_t other_path_num_reads = 0;
                for (const auto& it: other_path) {
                    other_path_num_reads += nodes_[it]->sequences.size();
                }

                std::vector<std::uint32_t> edges_for_removal;
                if (path_num_reads > other_path_num_reads) {
                    find_removable_edges(edges_for_removal, other_path);
                    if (edges_for_removal.empty()) {
                        find_removable_edges(edges_for_removal, path);
                    }
                } else {
                    find_removable_edges(edges_for_removal, path);
                    if (edges_for_removal.empty()) {
                        find_removable_edges(edges_for_removal, other_path);
                    }
                }

                for (const auto& edge_id: edges_for_removal) {
                    edges_[edge_id]->state |= 1;
                    edges_[edge_id]->pair->state |= 1;
                    marked_edges_.emplace(edge_id);
                    marked_edges_.emplace(edges_[edge_id]->pair->id);
                }
                if (!edges_for_removal.empty()) {
                    remove_marked_objects(true);
                    ++num_bubbles_popped;
                }
            }
        }

        node_queue.clear();
        for (std::uint32_t i = 0; i < visited_length; ++i) {
            distance[visited[i]] = 0;
            predecessor[visited[i]] = -1;
        }
        visited_length = 0;
    }

    return num_bubbles_popped;
}

std::uint32_t Graph::remove_long_edges() {

    std::uint32_t num_long_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->outedges.size() < 2){
            continue;
        }

        for (const auto& edge: node->outedges) {
            for (const auto& other_edge: node->outedges) {
                if (edge->id == other_edge->id || edge->state & 1 || other_edge->state & 1) {
                    continue;
                }
                if (edge->weight * 2.0 < other_edge->weight) {
                    other_edge->state |= 1;
                    other_edge->pair->state |= 1;
                    marked_edges_.emplace(other_edge->id);
                    marked_edges_.emplace(other_edge->pair->id);
                    ++num_long_edges;
                }
            }
        }
    }

    remove_marked_objects();

    return num_long_edges;
}

void Graph::create_force_directed_layout(const std::string& path) {

    std::ofstream os;
    bool is_first = true;
    if (!path.empty()) {
        os.open(path);
        os << "{" << std::endl;
        os << "  \"assembly\": {" << std::endl;
    }

    std::vector<std::unordered_set<std::uint32_t>> components;
    std::vector<char> is_visited(piles_.size(), 0);
    for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i] == nullptr || is_visited[i]) {
            continue;
        }

        components.resize(components.size() + 1);

        std::deque<std::uint32_t> que = { i };
        while (!que.empty()) {
            std::uint32_t j = que.front();
            que.pop_front();

            if (is_visited[j]) {
                continue;
            }
            const auto& node = nodes_[j];
            is_visited[node->id] = 1;
            is_visited[node->pair->id] = 1;
            components.back().emplace((node->id >> 1) << 1);

            for (const auto& it: node->inedges) {
                que.emplace_back(it->begin->id);
            }
            for (const auto& it: node->outedges) {
                que.emplace_back(it->end->id);
            }
        }
    }
    std::vector<char>().swap(is_visited);

    std::sort(components.begin(), components.end(),
        [] (const std::unordered_set<std::uint32_t>& lhs,
            const std::unordered_set<std::uint32_t>& rhs) {
            return lhs.size() > rhs.size();
        }
    );

    std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<> distribution(0., 1.);

    using Point = std::pair<double, double>;

    auto point_add = [] (const Point& x, const Point& y) -> Point {
        return std::make_pair(x.first + y.first, x.second + y.second);
    };
    auto point_substract = [] (const Point& x, const Point& y) -> Point {
       return std::make_pair(x.first - y.first, x.second - y.second);
    };
    auto point_multiply = [] (const Point& x, double s) -> Point {
        return std::make_pair(x.first * s, x.second * s);
    };
    auto point_norm = [] (const Point& x) -> double {
        return std::sqrt(x.first * x.first + x.second * x.second);
    };

    std::uint32_t c = 0;
    for (const auto& component: components) {

        if (component.size() < 6) {
            continue;
        }

        bool has_junctions = false;
        for (const auto& it: component) {
            if (nodes_[it]->is_junction()) {
                has_junctions = true;
                break;
            }
        }
        if (has_junctions == false) {
            continue;
        }

        // update transitive edges
        for (const auto& n: component) {
            std::unordered_set<std::uint32_t> valid;
            for (const auto& m: nodes_[n]->transitive) {
                if (component.find(m) != component.end()) {
                    valid.emplace(m);
                }
            }
            nodes_[n]->transitive.swap(valid);
        }

        std::uint32_t num_iterations = 100;
        double k = sqrt(1. / static_cast<double>(component.size()));
        double t = 0.1;
        double dt = t / static_cast<double>(num_iterations + 1);

        std::vector<Point> points(nodes_.size());
        for (const auto& it: component) {
            points[it].first = distribution(generator);
            points[it].second = distribution(generator);
        }

        for (std::uint32_t i = 0; i < num_iterations; ++i) {
            std::vector<std::future<void>> thread_futures;
            std::vector<Point> displacements(nodes_.size());

            auto thread_task = [&](std::uint32_t n) -> void {
                Point displacement = {0., 0.};
                for (const auto& m: component) {
                    if (n == m) {
                        continue;
                    }
                    auto delta = point_substract(points[n], points[m]);
                    auto distance = point_norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = point_add(displacement,
                        point_multiply(delta, (k * k) / (distance * distance)));
                }
                for (const auto& e: nodes_[n]->inedges) {
                    auto m = (e->begin->id >> 1) << 1;
                    auto delta = point_substract(points[n], points[m]);
                    auto distance = point_norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = point_add(displacement,
                        point_multiply(delta, -1. * distance / k));
                }
                for (const auto& e: nodes_[n]->outedges) {
                    auto m = (e->end->id >> 1) << 1;
                    auto delta = point_substract(points[n], points[m]);
                    auto distance = point_norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = point_add(displacement,
                        point_multiply(delta, -1. * distance / k));
                }
                for (const auto& m: nodes_[n]->transitive) {
                    auto delta = point_substract(points[n], points[m]);
                    auto distance = point_norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = point_add(displacement,
                        point_multiply(delta, -1. * distance / k));
                }
                auto length = point_norm(displacement);
                if (length < 0.01) {
                    length = 0.1;
                }
                displacements[n] = point_add(displacements[n],
                    point_multiply(displacement, t / length));
                return;
            };

            for (const auto& n: component) {
                thread_futures.emplace_back(thread_pool_->submit(thread_task, n));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
            for (const auto& n: component) {
                points[n] = point_add(points[n], displacements[n]);
            }

            t -= dt;
            ++i;
        }

        for (const auto& it: edges_) {
            if (it == nullptr || it->id & 1) {
                continue;
            }
            auto n = (it->begin->id >> 1) << 1;
            auto m = (it->end->id >> 1) << 1;

            if (component.find(n) != component.end() &&
                component.find(m) != component.end()) {
                it->weight = point_norm(point_substract(points[n], points[m]));
                it->pair->weight = it->weight;
            }
        }

        if (!path.empty()) {
            if (!is_first) {
                os << "," << std::endl;
            }
            is_first = false;

            os << "    \"component_" << c++ << "\": {" << std::endl;

            bool is_first_node = true;
            os << "      \"nodes\": {" << std::endl;
            for (const auto& it: component) {
                if (!is_first_node) {
                    os << "," << std::endl;
                }
                is_first_node = false;
                os << "        \"" << it << "\": [";
                os << points[it].first << ", ";
                os << points[it].second << ", ";
                os << (nodes_[it]->is_junction() ? 1 : 0) << ", ";
                os << nodes_[it]->sequences.size() << "]";
            }
            os << std::endl << "      }," << std::endl;

            bool is_first_edge = true;
            os << "      \"edges\": [" << std::endl;
            for (const auto& it: component) {
                for (const auto& e: nodes_[it]->inedges) {
                    auto o = (e->begin->id >> 1) << 1;
                    if (it < o) {
                        continue;
                    }
                    if (!is_first_edge) {
                        os << "," << std::endl;
                    }
                    is_first_edge = false;
                    os << "        [\"" << it << "\", \"" << o << "\", 0]";
                }
                for (const auto& e: nodes_[it]->outedges) {
                    auto o = (e->end->id >> 1) << 1;
                    if (it < o) {
                        continue;
                    }
                    if (!is_first_edge) {
                        os << "," << std::endl;
                    }
                    is_first_edge = false;
                    os << "        [\"" << it << "\", \"" << o << "\", 0]";
                }
                for (const auto& o: nodes_[it]->transitive) {
                    if (it < o) {
                        continue;
                    }
                    if (!is_first_edge) {
                        os << "," << std::endl;
                    }
                    is_first_edge = false;
                    os << "        [\"" << it << "\", \"" << o << "\", 1]";
                }
            }
            os << std::endl << "      ]" << std::endl;

            os << "    }";
        }
    }

    if (!path.empty()) {
        os << std::endl << "  }";
        os << std::endl << "}";
        os << std::endl;
        os.close();
    }
}

std::uint32_t Graph::create_unitigs(std::uint32_t epsilon) {

    std::vector<char> is_visited(nodes_.size(), 0);
    std::vector<std::uint32_t> node_updates(nodes_.size(), 0);

    std::uint32_t node_id = nodes_.size();
    std::vector<std::unique_ptr<Node>> unitigs;

    std::uint32_t edge_id = edges_.size();
    std::vector<std::unique_ptr<Edge>> unitig_edges;

    std::uint32_t num_unitigs_created = 0;

    for (const auto& it: nodes_) {
        if (it == nullptr || is_visited[it->id] || it->is_junction()) {
            continue;
        }

        std::uint32_t extension = 1;

        bool is_circular = false;
        auto begin = it.get();
        while (!begin->is_junction()) {
            is_visited[begin->id] = 1;
            is_visited[begin->pair->id] = 1;
            if (begin->indegree() == 0 || begin->inedges[0]->begin->is_junction()) {
                break;
            }
            begin = begin->inedges[0]->begin;
            ++extension;
            if (begin->id == it->id) {
                is_circular = true;
                break;
            }
        }

        auto end = it.get();
        while (!end->is_junction()) {
            is_visited[end->id] = 1;
            is_visited[end->pair->id] = 1;
            if (end->outdegree() == 0 || end->outedges[0]->end->is_junction()) {
                break;
            }
            end = end->outedges[0]->end;
            ++extension;
            if (end->id == it->id) {
                is_circular = true;
                break;
            }
        }

        if (!is_circular && begin == end) {
            continue;
        }
        if (!is_circular && extension < 2 * epsilon + 2) {
            continue;
        }

        if (begin != end) {
            // update begin_node
            for (std::uint32_t i = 0; i < epsilon; ++i) {
                begin = begin->outedges[0]->end;
            }
            // update end_node
            for (std::uint32_t i = 0; i < epsilon; ++i) {
                end = end->inedges[0]->begin;
            }
        }

        std::unique_ptr<Node> u(new Node(node_id++, begin, end));
        std::unique_ptr<Node> uc(new Node(node_id++, end->pair, begin->pair));

        u->pair = uc.get();
        uc->pair = u.get();

        // update node identifiers for transitive edges
        auto node = begin;
        while (node != end) {
            node_updates[(node->id >> 1) << 1] = node_id - 2;
            u->transitive.insert(
                nodes_[(node->id >> 1) << 1]->transitive.begin(),
                nodes_[(node->id >> 1) << 1]->transitive.end());
            node = node->outedges[0]->end;
        }

        if (begin != end) {
            if (begin->indegree() != 0) {
                const auto& edge = begin->inedges[0];

                edge->state |= 1;
                edge->pair->state |= 1;
                marked_edges_.emplace(edge->id);
                marked_edges_.emplace(edge->pair->id);

                std::unique_ptr<Edge> e(new Edge(edge_id++, edge->begin,
                    u.get(), edge->length));
                std::unique_ptr<Edge> ec(new Edge(edge_id++, uc.get(),
                    edge->pair->end, edge->pair->length + uc->length() -
                    begin->pair->length()));

                e->pair = ec.get();
                ec->pair = e.get();

                edge->begin->outedges.emplace_back(e.get());
                edge->pair->end->inedges.emplace_back(ec.get());
                u->inedges.emplace_back(e.get());
                uc->outedges.emplace_back(ec.get());

                unitig_edges.emplace_back(std::move(e));
                unitig_edges.emplace_back(std::move(ec));
            }
            if (end->outdegree() != 0) {
                const auto& edge = end->outedges[0];

                edge->state |= 1;
                edge->pair->state |= 1;
                marked_edges_.emplace(edge->id);
                marked_edges_.emplace(edge->pair->id);

                std::unique_ptr<Edge> e(new Edge(edge_id++, u.get(), edge->end,
                    edge->length + u->length() - end->length()));
                std::unique_ptr<Edge> ec(new Edge(edge_id++, edge->pair->begin,
                    uc.get(), edge->pair->length));

                e->pair = ec.get();
                ec->pair = e.get();

                u->outedges.emplace_back(e.get());
                uc->inedges.emplace_back(ec.get());
                edge->end->inedges.emplace_back(e.get());
                edge->pair->begin->outedges.emplace_back(ec.get());

                unitig_edges.emplace_back(std::move(e));
                unitig_edges.emplace_back(std::move(ec));
            }
        }

        unitigs.emplace_back(std::move(u));
        unitigs.emplace_back(std::move(uc));

        ++num_unitigs_created;

        // mark edges for deletion
        node = begin;
        while (true) {
            const auto& edge = node->outedges[0];

            edge->state |= 1;
            edge->pair->state |= 1;
            marked_edges_.emplace(edge->id);
            marked_edges_.emplace(edge->pair->id);

            node = edge->end;
            if (node == end) {
                break;
            }
        }
    }

    for (std::uint32_t i = 0; i < unitigs.size(); ++i) {
        nodes_.emplace_back(std::move(unitigs[i]));
    }
    for (std::uint32_t i = 0; i < unitig_edges.size(); ++i) {
        edges_.emplace_back(std::move(unitig_edges[i]));
    }

    remove_marked_objects(true);

    // update transitive edges
    for (const auto& it: nodes_) {
        if (it == nullptr) {
            continue;
        }
        std::unordered_set<std::uint32_t> valid;
        for (const auto& jt: it->transitive) {
            valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
        }
        it->transitive.swap(valid);
    }

    return num_unitigs_created;
}

void Graph::extract_unitigs(std::vector<std::unique_ptr<ram::Sequence>>& dst) {

    create_unitigs();

    ram::Sequence::num_objects = 0;

    std::uint32_t contig_id = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->is_rc()) {
            continue;
        }
        if (node->sequences.size() < 6 || node->length() < 10000) {
            continue;
        }

        std::string name = "Ctg" + std::to_string(contig_id);
        name += " RC:i:" + std::to_string(node->sequences.size());
        name += " LN:i:" + std::to_string(node->data.size());

        dst.emplace_back(std::unique_ptr<ram::Sequence>(new ram::Sequence(name, node->data)));
        ++contig_id;
    }
}

void Graph::remove_marked_objects(bool remove_nodes) {

    auto delete_edges = [&](std::vector<Edge*>& edges) -> void {
        std::uint32_t j = 0;
        for (std::uint32_t i = 0; i < edges.size(); ++i) {
            if (edges[i]->state & 1) {
                continue;
                edges[i] = nullptr;
            }
            edges[j++] = edges[i];
        }
        edges.resize(j);
    };

    std::unordered_set<std::uint32_t> marked_nodes;
    for (const auto& it: marked_edges_) {
        if (remove_nodes) {
            marked_nodes.emplace(edges_[it]->begin->id);
            marked_nodes.emplace(edges_[it]->end->id);
        }
        delete_edges(edges_[it]->begin->outedges);
        delete_edges(edges_[it]->end->inedges);
    }

    if (remove_nodes) {
        for (const auto& it: marked_nodes) {
            if (nodes_[it]->outdegree() == 0 && nodes_[it]->indegree() == 0) {
                nodes_[it].reset();
            }
        }
    }

    for (const auto& it: marked_edges_) {
        edges_[it].reset();
    }
    marked_edges_.clear();
}

void Graph::find_removable_edges(std::vector<std::uint32_t>& dst,
    const std::vector<std::uint32_t>& path) {

    if (path.empty()) {
        return;
    }

    auto find_edge = [&] (std::uint32_t src, std::uint32_t dst) -> std::uint32_t {
        std::uint32_t edge_id = 0;
        bool found_edge = false;
        for (const auto& edge: nodes_[src]->outedges) {
            if (edge->end->id == dst) {
                edge_id = edge->id;
                found_edge = true;
                break;
            }
        }

        if (!found_edge) {
            throw std::logic_error("[raven::Graph::find_removable_edges] error: "
                "missing edge between nodes");
        }

        return edge_id;
    };

    // find first node with multiple in edges
    std::int32_t pref = -1;
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
        if (nodes_[path[i]]->indegree() > 1) {
            pref = i;
            break;
        }
    }
    // find last node with multiple out edges
    std::int32_t suff = -1;
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
        if (nodes_[path[i]]->outdegree() > 1) {
            suff = i;
        }
    }

    if (pref == -1 && suff == -1) {
        // remove whole path
        for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
        return;
    }

    if (pref != -1 && nodes_[path[pref]]->outdegree() > 1) {
        return;
    }
    if (suff != -1 && nodes_[path[suff]]->indegree() > 1) {
        return;
    }

    if (pref == -1) {
        // remove everything after last suff node
        for (std::uint32_t i = suff; i < path.size() - 1; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
    } else if (suff == -1) {
        // remove everything before first pref node
        for (std::int32_t i = 0; i < pref; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
    } else if (suff < pref) {
        // remove everything between last suff and first pref node
        for (std::int32_t i = suff; i < pref; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
    }
}

void Graph::print_json(const std::string& path) const {

    std::ofstream os(path);
    os << "{\"piles\":{";
    bool is_first = true;

    for (const auto& it: piles_) {
        if (it->is_invalid()) {
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

void Graph::print_csv(const std::string& path) const {

    std::ofstream os(path);

    for (const auto& it: nodes_) {
        if (it == nullptr || it->is_rc() ||
            (it->outdegree() == 0 && it->indegree() == 0)) {
            continue;
        }
        os << it->id << "[" << it->sequences.front() << "] LN:i:";
        os << it->length() << " RC:i:" << it->sequences.size() << ",";
        os << it->pair->id << "[" << it->pair->sequences.front() << "] LN:i:";
        os << it->pair->length() << " RC:i:" << it->pair->sequences.size() << ",";
        os << "0,-";
        os << std::endl;
    }

    for (const auto& it: edges_) {
        if (it == nullptr) {
            continue;
        }
        os << it->begin->id << "[" << it->begin->sequences.front() << "] LN:i:";
        os << it->begin->length() << " RC:i:" << it->begin->sequences.size() << ",";
        os << it->end->id << "[" << it->end->sequences.front() << "] LN:i:";
        os << it->end->length() << " RC:i:" << it->end->sequences.size() << ",";
        os << "1,";
        os << it->id << " ";
        os << it->length << " ";
        os << it->weight;
        os << std::endl;
    }

    os.close();
}

Graph::Node::Node(std::uint32_t id, std::uint32_t sequence, const std::string& name,
    const std::string& data)
        : id(id), name(name), data(data), state(0), sequences(1, sequence),
        inedges(), outedges(), pair(nullptr) {
    state |= (id & 1) << 1;
    state |= (id & 1) << 2;
}

Graph::Node::Node(std::uint32_t id, Node* begin, Node* end)
        : id(id), name(), data(), state(0), sequences(), inedges(), outedges(),
        transitive(), pair(nullptr) {

    if (begin == nullptr) {
        throw std::invalid_argument("[raven::Graph::Node::Node] error: "
            "begin node is nullptr!");
    }
    if (end == nullptr) {
        throw std::invalid_argument("[raven::Graph::Node::Node] error: "
            "end node is nullptr!");
    }

    state |= (begin->state & (1U << 1));

    auto node = begin;
    while (true) {
        auto edge = node->outedges[0];

        data += edge->label();
        sequences.insert(sequences.end(), node->sequences.begin(),
            node->sequences.end());
        state |= (node->state & (1U << 2));

        node = edge->end;
        if (node == end) {
            break;
        }
    }

    if (begin != end) {
        data += end->data;
        sequences.insert(sequences.end(), end->sequences.begin(),
            end->sequences.end());
        state |= (end->state & (1U << 2));
    }
}

Graph::Edge::Edge(std::uint32_t id, Node* begin, Node* end, std::uint32_t length)
        : id(id), length(length), state(0), weight(0), begin(begin), end(end),
        pair(nullptr) {
}

}
