#ifndef RAVEN_GRAPH_OVERLAP_UTILS_H_
#define RAVEN_GRAPH_OVERLAP_UTILS_H_

#include <deque>

#include "biosoup/overlap.hpp"
#include "biosoup/timer.hpp"
#include "raven/export.h"
#include "raven/pile.h"

namespace raven {

RAVEN_EXPORT biosoup::Overlap OverlapReverse(const biosoup::Overlap& o);

RAVEN_EXPORT std::uint32_t GetOverlapLength(const biosoup::Overlap& o);

RAVEN_EXPORT bool OverlapUpdate(biosoup::Overlap& o,
                   const std::vector<std::unique_ptr<Pile>>& piles);

RAVEN_EXPORT std::uint32_t GetOverlapType(const biosoup::Overlap& o,
                             const std::vector<std::unique_ptr<Pile>>& piles);

RAVEN_EXPORT bool OverlapFinalize(biosoup::Overlap& o,
                     const std::vector<std::unique_ptr<Pile>>& piles);

RAVEN_EXPORT std::vector<std::vector<std::uint32_t>> ConnectedComponents(
    const std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::vector<std::unique_ptr<Pile>>& piles);

}  // namespace raven

#endif  // RAVEN_GRAPH_OVERLAP_UTILS_H_
