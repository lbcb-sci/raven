#ifndef RAVEN_GRAPH_SERIALIZATION_GRAPH_REPR_H_
#define RAVEN_GRAPH_SERIALIZATION_GRAPH_REPR_H_

#include <fstream>

#include "cereal/archives/json.hpp"
#include "raven/export.h"
#include "raven/graph/graph.h"

namespace raven {

RAVEN_EXPORT void PrintGfa(const Graph& graph, const std::string& path);
RAVEN_EXPORT void PrintCsv(const Graph& graph, const std::string& path);
RAVEN_EXPORT void PrintJson(const Graph& graph, const std::string& path);

}  // namespace raven

#endif  // RAVEN_GRAPH_SERIALIZATION_GRAPH_REPR_H_
