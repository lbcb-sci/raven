#include "raven/graph//assemble.h"

#include <algorithm>
#include <deque>
#include <exception>
#include <fstream>
#include <future>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread_pool/thread_pool.hpp>

#include "raven/graph/common.h"
#include "raven/graph/serialization/binary.h"
#include "biosoup/timer.hpp"
#include "edlib.h"  // NOLINT
#include "ram/minimizer_engine.hpp"
#include "raven/graph/graph.h"

namespace raven {

static std::uint32_t RemoveTransitiveEdges(Graph& graph) {
  biosoup::Timer timer;
  timer.Start();

  auto is_comparable = [](double a, double b) -> bool {
    double eps = 0.12;
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
           (b >= a * (1 - eps) && b <= a * (1 + eps));
  };

  std::vector<Edge*> candidate(graph.nodes.size(), nullptr);
  std::unordered_set<std::uint32_t> marked_edges;

  for (const auto& it : graph.nodes) {
    if (it == nullptr) {
      continue;
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = jt;
    }
    for (auto jt : it->outedges) {
      for (auto kt : jt->head->outedges) {
        if (candidate[kt->head->id] &&
            is_comparable(jt->length + kt->length,
                          candidate[kt->head->id]->length)) {  // NOLINT
          marked_edges.emplace(candidate[kt->head->id]->id);
          marked_edges.emplace(candidate[kt->head->id]->pair->id);
        }
      }
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = nullptr;
    }
  }

  for (auto i : marked_edges) {  // store for force directed layout
    if (i & 1) {
      auto lhs = graph.edges[i]->tail->id & ~1UL;
      auto rhs = graph.edges[i]->head->id & ~1UL;
      graph.nodes[lhs]->transitive.emplace(rhs);
      graph.nodes[rhs]->transitive.emplace(lhs);
    }
  }

  RemoveEdges(graph, marked_edges);

  std::cerr << "[raven::Graph::Assemble] removed transitive edges "
            << std::fixed << timer.Stop() << "s" << std::endl;

  return marked_edges.size() / 2;
}

static std::uint32_t RemoveTips(Graph& graph) {
  std::uint32_t num_tips = 0;
  std::vector<char> is_visited(graph.nodes.size(), 0);

  for (const auto& it : graph.nodes) {
    if (it == nullptr || is_visited[it->id] || !it->is_tip()) {
      continue;
    }
    bool is_circular = false;
    std::uint32_t num_sequences = 0;

    auto end = it.get();
    while (!end->is_junction()) {
      num_sequences += end->count;
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 || end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    if (is_circular || end->outdegree() == 0 || num_sequences > 5) {
      continue;
    }

    std::unordered_set<std::uint32_t> marked_edges;
    for (auto jt : end->outedges) {
      if (jt->head->indegree() > 1) {
        marked_edges.emplace(jt->id);
        marked_edges.emplace(jt->pair->id);
      }
    }
    if (marked_edges.size() / 2 == end->outedges.size()) {  // delete whole
      auto begin = it.get();
      while (begin != end) {
        marked_edges.emplace(begin->outedges.front()->id);
        marked_edges.emplace(begin->outedges.front()->pair->id);
        begin = begin->outedges.front()->head;
      }
      ++num_tips;
    }

    RemoveEdges(graph, marked_edges, true);
  }

  return num_tips;
}

static std::unordered_set<std::uint32_t> FindRemovableEdges(
    const std::vector<Node*>& path) {
  if (path.empty()) {
    return std::unordered_set<std::uint32_t>{};
  }

  auto find_edge = [](Node* tail, Node* head) -> Edge* {
    for (auto it : tail->outedges) {
      if (it->head == head) {
        return it;
      }
    }
    return nullptr;
  };

  // find first node with multiple in edges
  std::int32_t pref = -1;
  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->indegree() > 1) {
      pref = i;
      break;
    }
  }

  // find last node with multiple out edges
  std::int32_t suff = -1;
  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->outdegree() > 1) {
      suff = i;
    }
  }

  std::unordered_set<std::uint32_t> dst;
  if (pref == -1 && suff == -1) {  // remove whole path
    for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
    return dst;
  }

  if (pref != -1 && path[pref]->outdegree() > 1) {  // complex path
    return dst;                                     // empty
  }
  if (suff != -1 && path[suff]->indegree() > 1) {  // complex path
    return dst;                                    // empty
  }

  if (pref == -1) {  // remove everything after last suffix node
    for (std::uint32_t i = suff; i < path.size() - 1; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  } else if (suff == -1) {  // remove everything before first prefix node
    for (std::int32_t i = 0; i < pref; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  } else if (suff < pref) {  // remove everything in between
    for (std::int32_t i = suff; i < pref; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  }
  return dst;  // empty
}

static std::uint32_t RemoveBubbles(Graph& graph) {
  std::vector<std::uint32_t> distance(graph.nodes.size(), 0);
  std::vector<Node*> predecessor(graph.nodes.size(), nullptr);

  // path helper functions
  auto path_extract = [&](Node* begin, Node* end) -> std::vector<Node*> {
    std::vector<Node*> dst;
    while (end != begin) {
      dst.emplace_back(end);
      end = predecessor[end->id];
    }
    dst.emplace_back(begin);
    std::reverse(dst.begin(), dst.end());
    return dst;
  };
  auto path_is_simple = [](const std::vector<Node*>& path) -> bool {
    if (path.empty()) {
      return false;
    }
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
      if (path[i]->is_junction()) {
        return false;  // complex
      }
    }
    return true;  // without branches
  };
  auto path_sequence = [](const std::vector<Node*>& path) -> std::string {
    std::string data{};
    for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
      for (auto it : path[i]->outedges) {
        if (it->head == path[i + 1]) {
          data += it->Label();
          break;
        }
      }
    }
    data += path.back()->sequence.InflateData();
    return data;
  };
  auto bubble_pop =
      [&](const std::vector<Node*>& lhs,
          const std::vector<Node*>& rhs) -> std::unordered_set<std::uint32_t> {
    if (lhs.empty() || rhs.empty()) {
      return std::unordered_set<std::uint32_t>{};
    }

    // check BFS result
    std::unordered_set<Node*> bubble;
    bubble.insert(lhs.begin(), lhs.end());
    bubble.insert(rhs.begin(), rhs.end());
    if (lhs.size() + rhs.size() - 2 != bubble.size()) {
      return std::unordered_set<std::uint32_t>{};
    }
    for (const auto& it : lhs) {
      if (bubble.count(it->pair) != 0) {
        return std::unordered_set<std::uint32_t>{};
      }
    }

    if (!path_is_simple(lhs) || !path_is_simple(rhs)) {  // complex path(s)
      // check poppability
      if (FindRemovableEdges(lhs).empty() && FindRemovableEdges(rhs).empty()) {
        return std::unordered_set<std::uint32_t>{};
      }

      // check similarity
      auto l = path_sequence(lhs);
      auto r = path_sequence(rhs);
      if (std::min(l.size(), r.size()) < std::max(l.size(), r.size()) * 0.8) {
        return std::unordered_set<std::uint32_t>{};
      }

      EdlibAlignResult result = edlibAlign(l.c_str(), l.size(), r.c_str(),
                                           r.size(), edlibDefaultAlignConfig());
      double score = 0;
      if (result.status == EDLIB_STATUS_OK) {
        score = 1 - result.editDistance /
                        static_cast<double>(std::max(l.size(), r.size()));
        edlibFreeAlignResult(result);
      }
      if (score < 0.8) {
        return std::unordered_set<std::uint32_t>{};
      }
    }

    std::uint32_t lhs_count = 0;
    for (auto jt : lhs) {
      lhs_count += jt->count;
    }
    std::uint32_t rhs_count = 0;
    for (auto jt : rhs) {
      rhs_count += jt->count;
    }
    auto marked_edges = FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
    if (marked_edges.empty()) {
      marked_edges = FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
    }
    return marked_edges;
  };
  // path helper functions

  std::uint32_t num_bubbles = 0;
  for (const auto& it : graph.nodes) {
    if (it == nullptr || it->outdegree() < 2) {
      continue;
    }

    // BFS
    Node* begin = it.get();
    Node* end = nullptr;
    Node* other_end = nullptr;
    std::deque<Node*> que{begin};
    std::vector<Node*> visited{1, begin};
    while (!que.empty() && !end) {
      auto jt = que.front();
      que.pop_front();

      for (auto kt : jt->outedges) {
        if (kt->head == begin) {  // cycle
          continue;
        }
        if (distance[jt->id] + kt->length > 500000) {  // out of reach
          continue;
        }
        distance[kt->head->id] = distance[jt->id] + kt->length;
        visited.emplace_back(kt->head);
        que.emplace_back(kt->head);

        if (predecessor[kt->head->id]) {  // found bubble
          end = kt->head;
          other_end = jt;
          break;
        }

        predecessor[kt->head->id] = jt;
      }
    }

    std::unordered_set<std::uint32_t> marked_edges;
    if (end) {
      auto lhs = path_extract(begin, end);
      auto rhs = path_extract(begin, other_end);
      rhs.emplace_back(end);
      marked_edges = bubble_pop(lhs, rhs);
    }

    for (auto jt : visited) {
      distance[jt->id] = 0;
      predecessor[jt->id] = nullptr;
    }

    RemoveEdges(graph, marked_edges, true);
    num_bubbles += 1 - marked_edges.empty();
  }

  return num_bubbles;
}

static void CreateForceDirectedLayout(
    std::shared_ptr<thread_pool::ThreadPool>& threadPool, const Graph& graph,
    const std::string& path = "") {
  std::ofstream os;
  bool is_first = true;
  if (!path.empty()) {
    os.open(path);
    os << "{" << std::endl;
  }

  std::vector<std::unordered_set<std::uint32_t>> components;
  std::vector<char> is_visited(graph.nodes.size(), 0);
  for (std::uint32_t i = 0; i < graph.nodes.size(); ++i) {
    if (graph.nodes[i] == nullptr || is_visited[i]) {
      continue;
    }

    components.resize(components.size() + 1);

    std::deque<std::uint32_t> que = {i};
    while (!que.empty()) {
      std::uint32_t j = que.front();
      que.pop_front();

      if (is_visited[j]) {
        continue;
      }
      const auto& node = graph.nodes[j];
      is_visited[node->id] = 1;
      is_visited[node->pair->id] = 1;
      components.back().emplace((node->id >> 1) << 1);

      for (auto it : node->inedges) {
        que.emplace_back(it->tail->id);
      }
      for (auto it : node->outedges) {
        que.emplace_back(it->head->id);
      }
    }
  }
  std::vector<char>().swap(is_visited);

  std::sort(components.begin(), components.end(),
            [](const std::unordered_set<std::uint32_t>& lhs,
               const std::unordered_set<std::uint32_t>& rhs) {
              return lhs.size() > rhs.size();
            });

  static std::uint64_t seed = 21;
  seed <<= 1;

  std::mt19937 generator(seed);
  std::uniform_real_distribution<> distribution(0., 1.);

  struct Point {
    Point() = default;
    Point(double x, double y) : x(x), y(y) {}

    bool operator==(const Point& other) const {
      return x == other.x && y == other.y;
    }
    Point operator+(const Point& other) const {
      return Point(x + other.x, y + other.y);
    }
    Point& operator+=(const Point& other) {
      x += other.x;
      y += other.y;
      return *this;
    }
    Point operator-(const Point& other) const {
      return Point(x - other.x, y - other.y);
    }
    Point operator*(double c) const { return Point(x * c, y * c); }
    Point& operator/=(double c) {
      x /= c;
      y /= c;
      return *this;
    }
    double Norm() const { return sqrt(x * x + y * y); }

    double x;
    double y;
  };

  struct Quadtree {
    Quadtree(Point nucleus, double width)
        : nucleus(nucleus), width(width), center(0, 0), mass(0), subtrees() {}

    bool Add(const Point& p) {
      if (nucleus.x - width > p.x || p.x > nucleus.x + width ||
          nucleus.y - width > p.y || p.y > nucleus.y + width) {
        return false;
      }
      ++mass;
      if (mass == 1) {
        center = p;
      } else if (subtrees.empty()) {
        if (center == p) {
          return true;
        }
        double w = width / 2;
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y - w), w);
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y - w), w);
        for (auto& it : subtrees) {
          if (it.Add(center)) {
            break;
          }
        }
      }
      for (auto& it : subtrees) {
        if (it.Add(p)) {
          break;
        }
      }
      return true;
    }

    void Centre() {
      if (subtrees.empty()) {
        return;
      }
      center = Point(0, 0);
      for (auto& it : subtrees) {
        it.Centre();
        center += it.center * it.mass;
      }
      center /= mass;
    }

    Point Force(const Point& p, double k) const {
      auto delta = p - center;
      auto distance = delta.Norm();
      if (width * 2 / distance < 1) {
        return delta * (mass * (k * k) / (distance * distance));
      }
      delta = Point(0, 0);
      for (const auto& it : subtrees) {
        delta += it.Force(p, k);
      }
      return delta;
    }

    Point nucleus;
    double width;
    Point center;
    std::uint32_t mass;
    std::vector<Quadtree> subtrees;
  };

  std::uint32_t c = 0;
  for (const auto& component : components) {
    if (component.size() < 6) {
      continue;
    }

    bool has_junctions = false;
    for (const auto& it : component) {
      if (graph.nodes[it]->is_junction()) {
        has_junctions = true;
        break;
      }
    }
    if (has_junctions == false) {
      continue;
    }

    // update transitive edges
    for (const auto& n : component) {
      std::unordered_set<std::uint32_t> valid;
      for (const auto& m : graph.nodes[n]->transitive) {
        if (component.find(m) != component.end()) {
          valid.emplace(m);
        }
      }
      graph.nodes[n]->transitive.swap(valid);
    }

    std::uint32_t num_iterations = 100;
    double k = sqrt(1. / static_cast<double>(component.size()));
    double t = 0.1;
    double dt = t / static_cast<double>(num_iterations + 1);

    std::vector<Point> points(graph.nodes.size());
    for (const auto& it : component) {
      points[it].x = distribution(generator);
      points[it].y = distribution(generator);
    }

    for (std::uint32_t i = 0; i < num_iterations; ++i) {
      Point x = {0, 0}, y = {0, 0};
      for (const auto& n : component) {
        x.x = std::min(x.x, points[n].x);
        x.y = std::max(x.y, points[n].x);
        y.x = std::min(y.x, points[n].y);
        y.y = std::max(y.y, points[n].y);
      }
      double w = (x.y - x.x) / 2, h = (y.y - y.x) / 2;

      Quadtree tree(Point(x.x + w, y.x + h), std::max(w, h) + 0.01);
      for (const auto& n : component) {
        tree.Add(points[n]);
      }
      tree.Centre();

      std::vector<std::future<void>> thread_futures;
      std::vector<Point> displacements(graph.nodes.size(), Point(0, 0));

      auto thread_task = [&](std::uint32_t n) -> void {
        auto displacement = tree.Force(points[n], k);
        for (auto e : graph.nodes[n]->inedges) {
          auto m = (e->tail->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (auto e : graph.nodes[n]->outedges) {
          auto m = (e->head->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (const auto& m : graph.nodes[n]->transitive) {
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        auto length = displacement.Norm();
        if (length < 0.01) {
          length = 0.1;
        }
        displacements[n] = displacement * (t / length);
        return;
      };

      for (const auto& n : component) {
        thread_futures.emplace_back(threadPool->Submit(thread_task, n));
      }
      for (const auto& it : thread_futures) {
        it.wait();
      }
      for (const auto& n : component) {
        points[n] += displacements[n];
      }

      t -= dt;
    }

    for (const auto& it : graph.edges) {
      if (it == nullptr || it->id & 1) {
        continue;
      }
      auto n = (it->tail->id >> 1) << 1;
      auto m = (it->head->id >> 1) << 1;

      if (component.find(n) != component.end() &&
          component.find(m) != component.end()) {
        it->weight = (points[n] - points[m]).Norm();
        it->pair->weight = it->weight;
      }
    }

    if (!path.empty()) {
      if (!is_first) {
        os << "," << std::endl;
      }
      is_first = false;

      os << "    \"component_" << c++ << "\": {" << std::endl;

      bool is_first_node = true;
      os << "      \"nodes\": {" << std::endl;
      for (const auto& it : component) {
        if (!is_first_node) {
          os << "," << std::endl;
        }
        is_first_node = false;
        os << "        \"" << it << "\": [";
        os << points[it].x << ", ";
        os << points[it].y << ", ";
        os << (graph.nodes[it]->is_junction() ? 1 : 0) << ", ";
        os << graph.nodes[it]->count << "]";
      }
      os << std::endl << "      }," << std::endl;

      bool is_first_edge = true;
      os << "      \"edges\": [" << std::endl;
      for (const auto& it : component) {
        for (auto e : graph.nodes[it]->inedges) {
          auto o = (e->tail->id >> 1) << 1;
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 0]";
        }
        for (auto e : graph.nodes[it]->outedges) {
          auto o = (e->head->id >> 1) << 1;
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 0]";
        }
        for (const auto& o : graph.nodes[it]->transitive) {
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 1]";
        }
      }
      os << std::endl << "      ]" << std::endl;
      os << "    }";
    }
  }

  if (!path.empty()) {
    os << std::endl << "}";
    os << std::endl;
    os.close();
  }
}

static std::uint32_t RemoveLongEdges(
    std::shared_ptr<thread_pool::ThreadPool>& threadPool, Graph& graph,
    std::uint32_t num_rounds) {
  std::uint32_t num_long_edges = 0;

  for (std::uint32_t i = 0; i < num_rounds; ++i) {
    CreateForceDirectedLayout(threadPool, graph);

    std::unordered_set<std::uint32_t> marked_edges;
    for (const auto& it : graph.nodes) {
      if (it == nullptr || it->outdegree() < 2) {
        continue;
      }
      for (auto jt : it->outedges) {
        for (auto kt : it->outedges) {
          if (jt != kt && jt->weight * 2.0 < kt->weight) {  // TODO(rvaser)
            marked_edges.emplace(kt->id);
            marked_edges.emplace(kt->pair->id);
          }
        }
      }
    }

    RemoveEdges(graph, marked_edges);
    num_long_edges += marked_edges.size() / 2;

    RemoveTips(graph);
  }

  return num_long_edges;
}

static std::uint32_t SalvagePlasmids(
    std::shared_ptr<thread_pool::ThreadPool>& threadPool, Graph& graph) {
  CreateUnitigs(graph);

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> plasmids;
  for (const auto& it : graph.nodes) {
    if (!it || it->is_rc() || it->is_unitig || !it->is_circular) {
      continue;
    }
    plasmids.emplace_back(new biosoup::NucleicAcid(it->sequence));
  }

  // remove duplicates within plasmids
  std::sort(plasmids.begin(), plasmids.end(),
            [](const std::unique_ptr<biosoup::NucleicAcid>& lhs,
               const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
              return lhs->inflated_len < rhs->inflated_len;
            });
  for (std::uint32_t i = 0; i < plasmids.size(); ++i) {
    plasmids[i]->id = i;
  }
  ram::MinimizerEngine minimizer_engine{threadPool};
  minimizer_engine.Minimize(plasmids.begin(), plasmids.end());
  minimizer_engine.Filter(0.001);
  for (auto& it : plasmids) {
    if (!minimizer_engine.Map(it, true, true).empty()) {
      it.reset();
    }
  }
  plasmids.erase(std::remove(plasmids.begin(), plasmids.end(), nullptr),
                 plasmids.end());

  if (plasmids.empty()) {
    return 0;
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
  for (const auto& it : graph.nodes) {
    if (!it || it->is_rc() || !it->is_unitig) {
      continue;
    }
    unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
  }

  // remove duplicates between plasmids and unitigs
  minimizer_engine.Minimize(unitigs.begin(), unitigs.end(), true);
  minimizer_engine.Filter(0.001);
  for (auto& it : plasmids) {
    if (!minimizer_engine.Map(it, false, false).empty()) {
      it.reset();
    }
  }
  plasmids.erase(std::remove(plasmids.begin(), plasmids.end(), nullptr),
                 plasmids.end());

  // update nodes
  for (const auto& it : plasmids) {
    const auto& node = graph.nodes[std::atoi(&it->name[3])];
    node->is_unitig = node->pair->is_unitig = true;
    node->sequence.name[0] = node->pair->sequence.name[0] = 'U';
  }

  return plasmids.size();
}

static void RemoveTipsAndBubbles(Graph& graph) {
  biosoup::Timer timer;
  timer.Start();

  while (true) {
    std::uint32_t num_changes = RemoveTips(graph);
    num_changes += RemoveBubbles(graph);
    if (num_changes == 0) {
      break;
    }
  }

  std::cerr << "[raven::Graph::Assemble] removed tips and bubbles "
            << std::fixed << timer.Stop() << "s" << std::endl;
}

static void RemoveLongEdgesStage(
    Graph& graph, std::shared_ptr<thread_pool::ThreadPool>& threadPool) {
  biosoup::Timer timer;
  timer.Start();

  CreateUnitigs(graph, 42);  // speed up force directed layout
  RemoveLongEdges(threadPool, graph, 16);

  std::cerr << "[raven::Graph::Assemble] removed long edges " << std::fixed
            << timer.Stop() << "s" << std::endl;

  timer.Start();

  while (true) {
    std::uint32_t num_changes = RemoveTips(graph);
    num_changes += RemoveBubbles(graph);
    if (num_changes == 0) {
      break;
    }
  }

  SalvagePlasmids(threadPool, graph);

  timer.Stop();
}

template <typename Fun, typename... Args>
static void StageExecution(Graph& graph, bool checkpoints, Fun fun,
                           Args&... args) {
  biosoup::Timer timer;
  timer.Start();

  fun(graph, args...);

  ++graph.stage;

  if (checkpoints) {
    timer.Start();
    raven::StoreGraphToFile(graph);
    std::cerr << "[raven::Graph::Assemble] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }
}

void Assemble(std::shared_ptr<thread_pool::ThreadPool> threadPool, Graph& graph,
              bool checkpoints) {
  if (graph.stage < -3 || graph.stage > -1) {
    return;
  }

  biosoup::Timer timer;
  timer.Start();

  if (graph.stage == -3) {  // remove transitive edges
    StageExecution(graph, checkpoints, RemoveTransitiveEdges);
  }

  if (graph.stage == -2) {  // remove tips and bubbles
    StageExecution(graph, checkpoints, RemoveTipsAndBubbles);
  }

  if (graph.stage == -1) {
    StageExecution(graph, checkpoints, RemoveLongEdgesStage, threadPool);
  }

  std::cerr << "[raven::Graph::Assemble] " << std::fixed << timer.Stop() << "s"
            << std::endl;
}

std::uint32_t RemoveTransitiveEdgesFromGraph(Graph& graph) {
  return RemoveTransitiveEdges(graph);
}

void RemoveTipsAndBubblesFromGraph(Graph& graph) {
  RemoveTipsAndBubbles(graph);
}

void RemoveLongEdgesFromGraph(Graph& graph, std::shared_ptr<thread_pool::ThreadPool>& threadPool) {
  RemoveLongEdgesStage(graph, threadPool);
}

}  // namespace raven
