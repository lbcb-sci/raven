#pragma once

#include "../Graph.hpp"

namespace raven {

void StoreGraphToFile(const Graph& graph);
Graph LoadGraphFromFile();

}  // namespace raven
