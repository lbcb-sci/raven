#include "raven/graph/graph.h"
#include "thread_pool/thread_pool.hpp"

namespace raven {

// Overlap phase
// - split chimeric sequences
// - remove contained sequences
// - remove overlaps not spanning bridged repeats at sequence ends

struct OverlapPhaseCfg {
  std::uint8_t kmer_len = 15;
  std::uint8_t window_len = 5;
  double freq = 0.001;
};

void ConstructGraph(
    Graph& graph, std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::shared_ptr<thread_pool::ThreadPool>& thread_pool, bool checkpoints,
    OverlapPhaseCfg const cfg);

}  // namespace raven
