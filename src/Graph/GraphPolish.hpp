#pragma once

#include "./Graph/Graph.hpp"
#include "thread_pool/thread_pool.hpp"

namespace raven {

// Consensus phase
// - utilize Racon

void polish(std::shared_ptr<thread_pool::ThreadPool> threadPool,
        Graph& graph,
        bool checkpoints,
        const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
        std::uint8_t match,
        std::uint8_t mismatch,
        std::uint8_t gap,
        std::uint32_t cuda_poa_batches,
        bool cuda_banded_alignment,
        std::uint32_t cuda_alignment_batches,
        std::uint32_t num_rounds);
}


