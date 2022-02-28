#pragma once

#include <biosoup/nucleic_acid.hpp>
#include <unordered_set>

#include "../pile.hpp"

namespace biosoup {

template <class Archive>
void serialize(Archive& archive, NucleicAcid& sequence) {  // NOLINT
  archive(sequence.id, sequence.name, sequence.deflated_data,
          sequence.block_quality, sequence.inflated_len,
          sequence.is_reverse_complement);
}

}  // namespace biosoup

namespace raven {

struct Graph {
  Graph();
  struct Edge;

  struct Node {
    Node() = default;  // needed for cereal

    explicit Node(const biosoup::NucleicAcid& sequence);
    Node(Node* begin, Node* end);

    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;

    Node(Node&&) = default;
    Node& operator=(Node&&) = default;

    ~Node() = default;

    std::uint32_t indegree() const { return inedges.size(); }
    std::uint32_t outdegree() const { return outedges.size(); }

    bool is_rc() const { return id & 1; }
    bool is_junction() const { return outdegree() > 1 || indegree() > 1; }
    bool is_tip() const {
      return outdegree() > 0 && indegree() == 0 && count < 6;
    }

    template <class Archive>
    void serialize(Archive& archive) {  // NOLINT
      archive(id, sequence, count, is_unitig, is_circular, is_polished,
              transitive);
    }

    static std::atomic<std::uint32_t> num_objects;

    std::uint32_t id;
    biosoup::NucleicAcid sequence;
    std::uint32_t count;
    bool is_unitig;
    bool is_circular;
    bool is_polished;
    std::unordered_set<std::uint32_t> transitive;
    std::vector<Edge*> inedges;
    std::vector<Edge*> outedges;
    Node* pair;
  };

  struct Edge {
    Edge() = default;  // needed for cereal
    Edge(Node* tail, Node* head, std::uint32_t length);

    Edge(const Edge&) = delete;
    Edge& operator=(const Edge&) = delete;

    ~Edge() = default;

    std::string Label() const { return tail->sequence.InflateData(0, length); }

    bool is_rc() const { return id & 1; }

    template <class Archive>
    void serialize(Archive& archive) {  // NOLINT
      archive(id, length, weight);
    }

    static std::atomic<std::uint32_t> num_objects;

    std::uint32_t id;
    std::uint32_t length;
    double weight;
    Node* tail;
    Node* head;
    Edge* pair;
  };

  int stage = -5;

  std::vector<std::unique_ptr<Pile>> piles;

  std::vector<std::shared_ptr<Node>> nodes;
  std::vector<std::shared_ptr<Edge>> edges;
};
}  // namespace raven
