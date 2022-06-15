#ifndef RAVEN_GRAPH_POLISH_H_
#define RAVEN_GRAPH_POLISH_H_

#include "raven/graph/graph.h"
#include "raven/export.h"
#include "thread_pool/thread_pool.hpp"

namespace raven {

// Consensus phase
// - utilize Racon

struct AlignCfg {
  std::int8_t match = 3;
  std::int8_t mismatch = -5;
  std::int8_t gap = -4;
};

struct CudaCfg {
  std::uint32_t poa_batches = 0U;
  std::uint32_t alignment_batches = 0U;
  bool banded_alignment = false;
};

struct PolishCfg {
  AlignCfg align_cfg;
  CudaCfg cuda_cfg;
  std::uint32_t num_rounds = 2U;
};

RAVEN_EXPORT void Polish(std::shared_ptr<thread_pool::ThreadPool> thread_pool, Graph& graph,
            bool checkpoints,
            const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
            PolishCfg cfg);
}  // namespace raven

#endif  // RAVEN_GRAPH_POLISH_H_
