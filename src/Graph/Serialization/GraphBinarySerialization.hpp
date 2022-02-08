#pragma once


#include "Graph/Graph.hpp"

namespace raven {

void storeGraphToFile(const Graph& graph);
Graph loadGraphFromFile();

}
