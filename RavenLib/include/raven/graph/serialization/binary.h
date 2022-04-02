#ifndef RAVEN_GRAPH_SERIALIZATION_H_
#define RAVEN_GRAPH_SERIALIZATION_H_

#include "raven/export.h"
#include "raven/graph/graph.h"

namespace raven {

RAVEN_EXPORT void StoreGraphToFile(const Graph& graph);
Graph LoadGraphFromFile();

}  // namespace raven

#endif  // RAVEN_GRAPH_SERIALIZATION_H_
