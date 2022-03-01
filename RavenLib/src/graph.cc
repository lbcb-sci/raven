#include "raven/graph/graph.h"

#include <string>

namespace raven {

Node::Node(const biosoup::NucleicAcid& sequence)
    : id(),
      sequence(sequence),
      count(1),
      is_unitig(),
      is_circular(),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {}

Node::Node(Node* begin, Node* end)
    : id(),
      sequence(),
      count(),
      is_unitig(),
      is_circular(begin == end),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {
  std::string data{};

  auto it = begin;
  while (true) {
    data += it->outedges.front()->Label();
    count += it->count;
    if ((it = it->outedges.front()->head) == end) {
      break;
    }
  }
  if (begin != end) {
    data += end->sequence.InflateData();
    count += end->count;
  }

  is_unitig = count > 5 && data.size() > 9999;

  sequence = biosoup::NucleicAcid(
      (is_unitig ? "Utg" : "Ctg") + std::to_string(id & (~1UL)), data);
}

Edge::Edge(Node* tail, Node* head, std::uint32_t length)
    : id(), length(length), weight(0), tail(tail), head(head), pair() {
  tail->outedges.emplace_back(this);
  head->inedges.emplace_back(this);
}

Graph::Graph() : stage(-5), piles(), nodes(), edges() {}

}  // namespace raven
