#ifndef RAVEN_GRAPH_GRAPH_H_
#define RAVEN_GRAPH_GRAPH_H_

#include <type_traits>
#include <unordered_set>

#include "biosoup/nucleic_acid.hpp"
#include "tsl/robin_set.h"
#include "raven/export.h"
#include "raven/pile.h"

namespace raven {

namespace detail {

template <class T>
class RAVEN_EXPORT IndexManager {
 public:
  using ValueType = T;
  static_assert(
      std::conjunction_v<std::is_integral<ValueType>,
                         std::negation<std::is_same<ValueType, bool>>>);

  explicit IndexManager(const ValueType init_val) : next_idx_val_(init_val) {}

  IndexManager(const IndexManager&) = delete;
  IndexManager& operator=(const IndexManager&) = delete;

  IndexManager(IndexManager&&) = default;
  IndexManager& operator=(IndexManager&&) = default;

  void Reset(const ValueType val) noexcept { next_idx_val_.store(val); }

  ValueType FetchAndIncrement() noexcept { return next_idx_val_++; }

  ValueType Fetch() const noexcept { return next_idx_val_; }

 private:
  ValueType next_idx_val_;
};

template <class I, class T>
class RAVEN_EXPORT IndexedFactory {
 private:
  using IndexManagerT = IndexManager<I>;

 public:
  using IndexType = typename IndexManagerT::ValueType;
  using ValueType = T;

  explicit IndexedFactory() : idx_manager_(0) {}

  explicit IndexedFactory(const IndexType init_index)
      : idx_manager_(init_index) {}

  IndexType NextIndex() const noexcept { return idx_manager_.Fetch(); }

  template <class... Args>
  ValueType Make(Args... args) {
    ValueType dst(std::forward<Args>(args)...);
    dst.id = idx_manager_.FetchAndIncrement();

    return dst;
  }

  template <class... Args>
  std::unique_ptr<ValueType> MakeUnique(Args... args) {
    auto dst = std::make_unique<ValueType>(idx_manager_.FetchAndIncrement(),
                                           std::forward<Args>(args)...);

    return dst;
  }

  void Reset() { idx_manager_.Reset(0); }

 private:
  IndexManagerT idx_manager_;
};

}  // namespace detail

struct RAVEN_EXPORT Edge;
struct RAVEN_EXPORT Node;

struct Node {
  Node() = default;  // needed for cereal

  explicit Node(const biosoup::NucleicAcid& sequence);
  explicit Node(Node* begin, Node* end);

  explicit Node(const std::uint32_t id, const biosoup::NucleicAcid& sequence);
  explicit Node(const std::uint32_t id, Node* begin, Node* end);

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

  std::uint32_t id;
  biosoup::NucleicAcid sequence;
  std::uint32_t count;
  bool is_unitig;
  bool is_circular;
  bool is_polished;
  tsl::robin_set<std::uint32_t> transitive;
  std::vector<Edge*> inedges;
  std::vector<Edge*> outedges;
  Node* pair;
};

struct Edge {
  Edge() = default;  // needed for cereal
  Edge(Node* tail, Node* head, std::uint32_t length);
  Edge(const std::uint32_t id, Node* tail, Node* head, std::uint32_t length);

  Edge(const Edge&) = delete;
  Edge& operator=(const Edge&) = delete;

  ~Edge() = default;

  std::string Label() const { return tail->sequence.InflateData(0, length); }

  bool is_rc() const { return id & 1; }

  std::uint32_t id;
  std::uint32_t length;
  double weight;
  Node* tail;
  Node* head;
  Edge* pair;
};

using NodeFactory = detail::IndexedFactory<std::uint32_t, Node>;
using EdgeFactory = detail::IndexedFactory<std::uint32_t, Edge>;

struct RAVEN_EXPORT Graph {
  Graph();

  int stage = -5;

  NodeFactory node_factory;
  EdgeFactory edge_factory;

  std::vector<std::unique_ptr<Pile>> piles;
  std::vector<std::unique_ptr<Node>> nodes;
  std::vector<std::unique_ptr<Edge>> edges;
};

}  // namespace raven

#endif  // RAVEN_GRAPH_GRAPH_H_
