#pragma once

#include <cstdint>

#include "Graph.hpp"

namespace raven {

// shrink graph by joining paths without branches into a single node
// - ignore nodes less than epsilon away from a node with outdegree > 1
std::uint32_t createUnitigs(Graph& graph, std::uint32_t epsilon = 0);

void removeEdges(Graph& graph, const std::unordered_set<std::uint32_t>& indices,
                 bool remove_nodes = false);

std::vector<std::unique_ptr<biosoup::NucleicAcid>> getUnitigs(
    Graph& graph, bool drop_unpolished = false);

}  // namespace raven
