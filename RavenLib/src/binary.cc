#include "raven/graph/serialization/binary.h"

#include <fstream>

#include "cereal/access.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/unordered_set.hpp"
#include "cereal/types/vector.hpp"

template <class Archive>
void serialize(Archive& archive, raven::Edge& edge) {
  archive(edge.id, edge.length, edge.weight);
}

template <class Archive>
void serialize(Archive& archive, raven::Node& node) {
  archive(node.id, node.sequence, node.count, node.is_unitig, node.is_circular,
          node.is_polished, node.transitive);
}

template <class Archive>
void save(Archive& archive, const raven::Graph& graph) {  // NOLINT

  std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;
  for (const auto& it : graph.edges) {
    if (it && !it->is_rc()) {
      connections.emplace_back(it->tail->id, it->head->id);
    } else {
      connections.emplace_back();  // dummy
    }
  }

  archive(graph.stage, graph.piles, graph.nodes, graph.edges, connections);
}

template <class Archive>
void load(Archive& archive, raven::Graph& graph) {
  std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;

  archive(graph.stage, graph.piles, graph.nodes, graph.edges, connections);

  for (std::uint32_t i = 0; i < graph.nodes.size(); i += 2) {
    if (graph.nodes[i]) {
      graph.nodes[i]->pair = graph.nodes[i + 1].get();
      graph.nodes[i + 1]->pair = graph.nodes[i].get();
    }
  }
  for (std::uint32_t i = 0; i < graph.edges.size(); i += 2) {
    if (graph.edges[i]) {
      graph.edges[i]->pair = graph.edges[i + 1].get();
      graph.edges[i + 1]->pair = graph.edges[i].get();

      graph.edges[i]->tail = graph.nodes[connections[i].first].get();
      graph.edges[i]->head = graph.nodes[connections[i].second].get();
      graph.edges[i + 1]->head = graph.edges[i]->tail->pair;
      graph.edges[i + 1]->tail = graph.edges[i]->head->pair;

      graph.edges[i]->tail->outedges.emplace_back(graph.edges[i].get());
      graph.edges[i]->head->inedges.emplace_back(graph.edges[i].get());
      graph.edges[i + 1]->tail->outedges.emplace_back(graph.edges[i + 1].get());
      graph.edges[i + 1]->head->inedges.emplace_back(graph.edges[i + 1].get());
    }
  }
}

namespace raven {

void StoreGraphToFile(const Graph& graph) {
  std::ofstream os("raven.cereal");
  try {
    cereal::BinaryOutputArchive archive(os);
    save(archive, graph);
  } catch (std::exception&) {
    throw std::logic_error(
        "[raven::Graph::Store] error: unable to store archive");
  }
}

Graph LoadGraphFromFile() {
  Graph graph;
  std::ifstream is("raven.cereal");
  try {
    cereal::BinaryInputArchive archive(is);
    load(archive, graph);
  } catch (std::exception&) {
    throw std::logic_error(
        "[raven::Graph::Load] error: unable to load archive");
  }

  return graph;
}
}  // namespace raven
