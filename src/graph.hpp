// Copyright (c) 2020 Robert Vaser

#ifndef RAVEN_GRAPH_HPP_
#define RAVEN_GRAPH_HPP_

#include <cstdint>
#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <utility>

#include <iostream>

#include "biosoup/sequence.hpp"
#include "cereal/access.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/unordered_set.hpp"
#include "cereal/types/vector.hpp"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

#include "pile.hpp"

namespace raven {

class Graph {
 public:
  Graph(
      bool weaken,
      bool checkpoints,
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph() = default;

  int stage() const {
    return stage_;
  }

  // break chimeric sequences, remove contained sequences and overlaps not
  // spanning bridged repeats at sequence ends
  void Construct(
      std::vector<std::unique_ptr<biosoup::Sequence>>& sequences,  // NOLINT
      bool split = false,
      std::string notations_path = "");

  // simplify with transitive reduction, tip prunning and bubble popping
  void Assemble();

  // Racon wrapper
  void Polish(
      const std::vector<std::unique_ptr<biosoup::Sequence>>& sequences,
      std::uint8_t match,
      std::uint8_t mismatch,
      std::uint8_t gap,
      std::uint32_t cuda_poa_batches,
      bool cuda_banded_alignment,
      std::uint32_t cuda_alignment_batches,
      std::uint32_t num_rounds);

  // ignore nodes that are less than epsilon away from any junction node
  std::uint32_t CreateUnitigs(std::uint32_t epsilon = 0);

  std::vector<std::unique_ptr<biosoup::Sequence>> GetUnitigs(
      bool drop_unpolished = false);

  // draw with misc/plotter.py
  void PrintJson(const std::string& path) const;

  // draw with Cytoscape
  void PrintCsv(const std::string& path) const;

  // draw with Bandage
  void PrintGfa(const std::string& path) const;

  // cereal load wrapper
  void Load();

  // cereal store wrapper
  void Store() const;

 private:
  // inspired by (Myers 1995) & (Myers 2005)
  std::uint32_t RemoveTransitiveEdges();

  std::uint32_t RemoveTips();

  std::uint32_t RemoveBubbles();

  // remove long edges in force directed layout
  std::uint32_t RemoveLongEdges(std::uint32_t num_round);

  friend cereal::access;

  Graph() = default;  // needed for cereal

  template<class Archive>
  void save(Archive& archive) const {  // NOLINT
    std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;
    for (const auto& it : edges_) {
      if (it && !it->is_rc()) {
        connections.emplace_back(it->tail->id, it->head->id);
      } else {
        connections.emplace_back();  // dummy
      }
    }

    archive(stage_, piles_, nodes_, edges_, connections);
  }

  template<class Archive>
  void load(Archive& archive) {  // NOLINT
    std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;

    archive(stage_, piles_, nodes_, edges_, connections);

    for (std::uint32_t i = 0; i < nodes_.size(); i += 2) {
      if (nodes_[i]) {
        nodes_[i]->pair = nodes_[i + 1].get();
        nodes_[i + 1]->pair = nodes_[i].get();
      }
    }
    for (std::uint32_t i = 0; i < edges_.size(); i += 2) {
      if (edges_[i]) {
        edges_[i]->pair = edges_[i + 1].get();
        edges_[i + 1]->pair = edges_[i].get();

        edges_[i]->tail = nodes_[connections[i].first].get();
        edges_[i]->head = nodes_[connections[i].second].get();
        edges_[i + 1]->head = edges_[i]->tail->pair;
        edges_[i + 1]->tail = edges_[i]->head->pair;

        edges_[i]->tail->outedges.emplace_back(edges_[i].get());
        edges_[i]->head->inedges.emplace_back(edges_[i].get());
        edges_[i + 1]->tail->outedges.emplace_back(edges_[i + 1].get());
        edges_[i + 1]->head->inedges.emplace_back(edges_[i + 1].get());
      }
    }

    Node::num_objects = nodes_.size();
    Edge::num_objects = edges_.size();
  }

  struct Node;
  struct Edge;

  struct Node {
   public:
    Node() = default;  // needed for cereal

    explicit Node(const biosoup::Sequence& sequence);
    Node(Node* begin, Node* end);

    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;

    Node(Node&&) = default;
    Node& operator=(Node&&) = default;

    ~Node() = default;

    std::uint32_t indegree() const {
      return inedges.size();
    }
    std::uint32_t outdegree() const {
      return outedges.size();
    }

    bool is_rc() const {
      return id & 1;
    }
    bool is_junction() const {
      return outdegree() > 1 || indegree() > 1;
    }
    bool is_tip() const {
      return outdegree() > 0 && indegree() == 0 && count < 6;
    }
    bool is_unitig() const {
      return count > 5 && data.size() > 9999;
    }

    template<class Archive>
    void serialize(Archive& archive) {  // NOLINT
      archive(id, name, data, count, is_circular, is_polished, transitive);
    }

    static std::atomic<std::uint32_t> num_objects;

    std::uint32_t id;
    std::string name;
    std::string data;
    std::uint32_t count;
    bool is_circular;
    bool is_polished;
    std::unordered_set<std::uint32_t> transitive;
    std::vector<Edge*> inedges;
    std::vector<Edge*> outedges;
    Node* pair;
  };
  struct Edge {
   public:
    Edge() = default;  // needed for cereal

    Edge(Node* tail, Node* head, std::uint32_t length);

    Edge(const Edge&) = delete;
    Edge& operator=(const Edge&) = delete;

    ~Edge() = default;

    bool is_rc() const {
      return id & 1;
    }

    std::string Label() const {
      return tail->data.substr(0, length);
    }

    template<class Archive>
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

  std::unordered_set<std::uint32_t> FindRemovableEdges(
      const std::vector<Node*>& path);

  void RemoveEdges(
      const std::unordered_set<std::uint32_t>& indices,
      bool remove_nodes = false);

  // use (Fruchterman & Reingold 1991) with (Barnes & Hut 1986) approximation
  // (draw with misc/plotter.py)
  void CreateForceDirectedLayout(const std::string& path = "");

  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  ram::MinimizerEngine minimizer_engine_;

  int stage_;
  bool checkpoints_;
  std::vector<std::unique_ptr<Pile>> piles_;
  std::vector<std::shared_ptr<Node>> nodes_;
  std::vector<std::shared_ptr<Edge>> edges_;
};

}  // namespace raven

#endif  // RAVEN_GRAPH_HPP_
