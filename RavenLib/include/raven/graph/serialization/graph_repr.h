#ifndef RAVEN_GRAPH_SERIALIZATION_GRAPH_REPR_H_
#define RAVEN_GRAPH_SERIALIZATION_GRAPH_REPR_H_

#include <fstream>
#include <map>

#include "cereal/archives/json.hpp"
#include "raven/export.h"
#include "raven/graph/graph.h"

namespace raven {

RAVEN_EXPORT void PrintGfa(const Graph& graph, const std::string& path);
RAVEN_EXPORT void PrintCsv(const Graph& graph, const std::string& path, bool printSequenceName = false, bool printPileBeginEnd = false, bool printEdgeSimilarity = false);
RAVEN_EXPORT void PrintJson(const Graph& graph, const std::string& path);
RAVEN_EXPORT Graph LoadGfa(const std::string& path);
RAVEN_EXPORT Graph LoadHifiasmGfa(const std::string& path, std::map<std::string, long>& nodeLengths, std::map<std::uint32_t, std::string>& additionalHifiasmEdgeInfo);
RAVEN_EXPORT void PrintHifiasmGfa(const Graph& graph, const std::string& path, std::map<std::string, long>& nodeLengths, std::map<std::uint32_t, std::string>& additionalHifiasmEdgeInfo);
RAVEN_EXPORT std::vector<std::string> getGfa(const Graph& graph);
RAVEN_EXPORT std::vector<std::string> getCsv(const Graph& graph, bool printSequenceName = false, bool printPileBeginEnd = false, bool printEdgeSimilarity = false);

}  // namespace raven

#endif  // RAVEN_GRAPH_SERIALIZATION_GRAPH_REPR_H_
