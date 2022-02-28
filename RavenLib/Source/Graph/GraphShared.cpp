#include "GraphShared.hpp"

namespace raven {

void removeEdges(Graph& graph, const std::unordered_set<std::uint32_t>& indices,
                 bool remove_nodes) {
  auto erase_remove = [](std::vector<Graph::Edge*>& edges,
                         Graph::Edge* marked) -> void {
    edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
  };

  std::unordered_set<std::uint32_t> node_indices;
  for (auto i : indices) {
    if (remove_nodes) {
      node_indices.emplace(graph.edges[i]->tail->id);
      node_indices.emplace(graph.edges[i]->head->id);
    }
    erase_remove(graph.edges[i]->tail->outedges, graph.edges[i].get());
    erase_remove(graph.edges[i]->head->inedges, graph.edges[i].get());
  }
  if (remove_nodes) {
    for (auto i : node_indices) {
      if (graph.nodes[i]->outdegree() == 0 && graph.nodes[i]->indegree() == 0) {
        graph.nodes[i].reset();
      }
    }
  }
  for (auto i : indices) {
    graph.edges[i].reset();
  }
}

std::uint32_t createUnitigs(Graph& graph, std::uint32_t epsilon) {
  std::unordered_set<std::uint32_t> marked_edges;
  std::vector<std::shared_ptr<Graph::Node>> unitigs;
  std::vector<std::shared_ptr<Graph::Edge>> unitig_edges;
  std::vector<std::uint32_t> node_updates(graph.nodes.size(), 0);
  std::vector<char> is_visited(graph.nodes.size(), 0);

  for (const auto& it : graph.nodes) {
    if (it == nullptr || is_visited[it->id] || it->is_junction()) {
      continue;
    }

    std::uint32_t extension = 1;

    bool is_circular = false;
    auto begin = it.get();
    while (!begin->is_junction()) {  // extend left
      is_visited[begin->id] = 1;
      is_visited[begin->pair->id] = 1;
      if (begin->indegree() == 0 ||
          begin->inedges.front()->tail->is_junction()) {
        break;
      }
      begin = begin->inedges.front()->tail;
      ++extension;
      if (begin == it.get()) {
        is_circular = true;
        break;
      }
    }

    auto end = it.get();
    while (!end->is_junction()) {  // extend right
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 || end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      ++extension;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    if (!is_circular && begin == end) {
      continue;
    }
    if (!is_circular && extension < 2 * epsilon + 2) {
      continue;
    }

    if (begin != end) {  // remove nodes near junctions
      for (std::uint32_t i = 0; i < epsilon; ++i) {
        begin = begin->outedges.front()->head;
      }
      for (std::uint32_t i = 0; i < epsilon; ++i) {
        end = end->inedges.front()->tail;
      }
    }

    auto unitig = std::make_shared<Graph::Node>(begin, end);
    unitigs.emplace_back(unitig);
    unitigs.emplace_back(std::make_shared<Graph::Node>(end->pair, begin->pair));
    unitig->pair = unitigs.back().get();
    unitig->pair->pair = unitig.get();

    if (begin != end) {  // connect unitig to graph
      if (begin->indegree()) {
        marked_edges.emplace(begin->inedges.front()->id);
        marked_edges.emplace(begin->inedges.front()->pair->id);

        auto edge = std::make_shared<Graph::Edge>(
            begin->inedges.front()->tail, unitig.get(),
            begin->inedges.front()->length);
        unitig_edges.emplace_back(edge);
        unitig_edges.emplace_back(std::make_shared<Graph::Edge>(
            unitig->pair, begin->inedges.front()->pair->head,
            begin->inedges.front()->pair->length +
                unitig->pair->sequence.inflated_len -
                begin->pair->sequence.inflated_len));  // NOLINT
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
      if (end->outdegree()) {
        marked_edges.emplace(end->outedges.front()->id);
        marked_edges.emplace(end->outedges.front()->pair->id);

        auto edge = std::make_shared<Graph::Edge>(
            unitig.get(), end->outedges.front()->head,
            end->outedges.front()->length + unitig->sequence.inflated_len -
                end->sequence.inflated_len);  // NOLINT
        unitig_edges.emplace_back(edge);
        unitig_edges.emplace_back(std::make_shared<Graph::Edge>(
            end->outedges.front()->pair->tail, unitig->pair,
            end->outedges.front()->pair->length));
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
    }

    auto jt = begin;
    while (true) {
      marked_edges.emplace(jt->outedges.front()->id);
      marked_edges.emplace(jt->outedges.front()->pair->id);

      // update transitive edges
      node_updates[jt->id & ~1UL] = unitig->id;
      unitig->transitive.insert(graph.nodes[jt->id & ~1UL]->transitive.begin(),
                                graph.nodes[jt->id & ~1UL]->transitive.end());

      if ((jt = jt->outedges.front()->head) == end) {
        break;
      }
    }
  }

  graph.nodes.insert(graph.nodes.end(), unitigs.begin(), unitigs.end());
  graph.edges.insert(graph.edges.end(), unitig_edges.begin(),
                     unitig_edges.end());
  removeEdges(graph, marked_edges, true);

  for (const auto& it : graph.nodes) {  // update transitive edges
    if (it) {
      std::unordered_set<std::uint32_t> valid;
      for (auto jt : it->transitive) {
        valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
      }
      it->transitive.swap(valid);
    }
  }

  return unitigs.size() / 2;
}

std::vector<std::unique_ptr<biosoup::NucleicAcid>> getUnitigs(
    Graph& graph, bool drop_unpolished) {
  createUnitigs(graph);

  biosoup::NucleicAcid::num_objects = 0;

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;

  for (const auto& it : graph.nodes) {
    if (it == nullptr || it->is_rc() || !it->is_unitig) {
      continue;
    }
    if (drop_unpolished && !it->is_polished) {
      continue;
    }

    std::string name = it->sequence.name +
                       " LN:i:" + std::to_string(it->sequence.inflated_len) +
                       " RC:i:" + std::to_string(it->count) +
                       " XO:i:" + std::to_string(it->is_circular);

    dst.emplace_back(
        new biosoup::NucleicAcid(name, it->sequence.InflateData()));
  }

  return dst;
}

}  // namespace raven
