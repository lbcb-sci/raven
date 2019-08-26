/*!
* @file graph.cpp
*
* @brief Graph class source file
*/

#include <exception>
#include <iostream>

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

    /*!
     * @brief ram::Overlap helper functions
     */
    auto overlap_reverse = [] (const ram::Overlap& o) -> ram::Overlap {
        return ram::Overlap(o.t_id, o.t_begin, o.t_end, o.q_id, o.q_begin, o.q_end,
            o.strand, o.matches);
    };
    auto overlap_length = [] (const ram::Overlap& o) -> std::uint32_t {
        return std::max(o.t_end - o.t_begin, o.q_end - o.q_begin);
    };

    for (const auto& it: sequences) {
        piles_.emplace_back(createPile(it->id, it->data.size()));
    }

    auto minimizer_engine = ram::createMinimizerEngine(15, 5, thread_pool_);

    for (std::uint32_t i = 0, j = 0, bytes = 0; i < sequences.size(); ++i) {
        bytes += sequences[i]->data.size();
        if (i != sequences.size() - 1 && bytes < (2U << 30)) {
            continue;
        }

        minimizer_engine->minimize(sequences.begin() + j, sequences.begin() + i + 1,
            0.001);

        std::vector<std::vector<ram::Overlap>> overlaps(sequences.size());
        {
            std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
            for (std::uint32_t k = 0; k < sequences.size(); ++k) {
                thread_futures.emplace_back(thread_pool_->submit(
                    [&minimizer_engine, &sequences] (std::uint32_t id) -> std::vector<ram::Overlap> {
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
                    [this, &overlaps] (std::uint32_t id) -> void {
                        piles_[id]->add_layers(overlaps[id]);
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
}

}
