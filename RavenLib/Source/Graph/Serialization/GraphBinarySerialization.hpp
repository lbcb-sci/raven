#pragma once

#include "../Graph.hpp"

namespace raven {

void storeGraphToFile(const Graph& graph);
Graph loadGraphFromFile();

}  // namespace raven
