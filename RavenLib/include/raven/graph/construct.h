#include "ram/minimizer_engine.hpp"
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
  double identity = 0.999;
};

void FindOverlapsAndCreatePiles(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    ram::MinimizerEngine& minimizer_engine,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    double freq, std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps);

void TrimAndAnnotatePiles(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps);

void ResolveContainedReads(
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps);

void ResolveContainedReads(
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps);

void ResolveChimericSequences(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);

void FindOverlapsAndRepetetiveRegions(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    ram::MinimizerEngine& minimizerEngine, double freq, std::uint8_t kmer_len,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);

void ResolveRepeatInducedOverlaps(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);

void ConstructAssemblyGraph(
    Graph& graph, const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);

void ConstructGraph(
    Graph& graph, std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::shared_ptr<thread_pool::ThreadPool>& thread_pool, bool checkpoints,
    OverlapPhaseCfg const cfg);

}  // namespace raven
