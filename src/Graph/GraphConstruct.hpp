#pragma once

#include <thread_pool/thread_pool.hpp>

#include "./Graph/Graph.hpp"

namespace raven {

// Overlap phase
// - split chimeric sequences
// - remove contained sequences
// - remove overlaps not spanning bridged repeats at sequence ends

void constructGraph(Graph& graph,
        std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
        std::shared_ptr<thread_pool::ThreadPool>& threadPool,
        bool checkpoints,
        std::uint8_t kmerLen = 15,
        std::uint8_t windowLen = 5,
        double freq = 0.001);
}

