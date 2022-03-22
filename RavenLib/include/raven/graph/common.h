#ifndef RAVEN_GRAPH_COMMON_H_
#define RAVEN_GRAPH_COMMON_H_

#include <cstdint>

#include "raven/export.h"
#include "raven/graph/graph.h"

namespace raven {

// shrink graph by joining paths without branches into a single node
// - ignore nodes less than epsilon away from a node with outdegree > 1
RAVEN_EXPORT std::uint32_t CreateUnitigs(Graph& graph, std::uint32_t epsilon = 0);

RAVEN_EXPORT void RemoveEdges(Graph& graph, const std::unordered_set<std::uint32_t>& indices,
                 bool remove_nodes = false);

RAVEN_EXPORT std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetUnitigs(
    Graph& graph, bool drop_unpolished = false);

}  // namespace raven

#endif  // RAVEN_GRAPH_COMMON_H_
