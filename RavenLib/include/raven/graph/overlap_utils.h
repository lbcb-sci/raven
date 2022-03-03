#ifndef RAVEN_GRAPH_OVERLAP_UTILS_H_
#define RAVEN_GRAPH_OVERLAP_UTILS_H_

#include <deque>

#include "biosoup/overlap.hpp"
#include "biosoup/timer.hpp"
#include "raven/pile.h"

namespace raven {

biosoup::Overlap OverlapReverse(const biosoup::Overlap& o);

std::uint32_t GetOverlapLength(const biosoup::Overlap& o);

bool OverlapUpdate(biosoup::Overlap& o,
                   const std::vector<std::unique_ptr<Pile>>& piles);

std::uint32_t GetOverlapType(const biosoup::Overlap& o,
                             const std::vector<std::unique_ptr<Pile>>& piles);

bool OverlapFinalize(biosoup::Overlap& o,
                     const std::vector<std::unique_ptr<Pile>>& piles);

std::vector<std::vector<std::uint32_t>> ConnectedComponents(
    const std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::vector<std::unique_ptr<Pile>>& piles);

}  // namespace raven

#endif  // RAVEN_GRAPH_OVERLAP_UTILS_H_
