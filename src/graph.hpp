// Copyright (c) 2020 Robert Vaser

#ifndef RAVEN_GRAPH_HPP_
#define RAVEN_GRAPH_HPP_

#include <cstdint>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <unordered_set>

#include "biosoup/sequence.hpp"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

#include "pile.hpp"

namespace raven {

class Graph {
 public:
  explicit Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);  // NOLINT

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph();

  /*!
   * @brief Constructs the assembly graph by breaking chimeric sequences,
   * removing contained sequences and removing overlaps between repetitive
   * sequences.
   */
  void Construct(std::vector<std::unique_ptr<biosoup::Sequence>>& sequences);

  /*!
   * @brief Simplify the assembly graph via transitive reduction, tip
   * pruning and popping bubble-like structures.
   */
  void Assemble();

  /*!
   * @brief Removes transitive edge (inspired by Myers 1995 & 2005).
   */
  std::uint32_t RemoveTransitiveEdges();

  /*!
   * @brief Removes nodes which are dead ends in graph.
   */
  std::uint32_t RemoveTips();

  /*!
   * @brief Removes bubble-like strucutres.
   */
  std::uint32_t RemoveBubbles();

  /*!
   * @brief Removes long edges in the force directed layout which function
   * needs to be called before.
   */
  std::uint32_t RemoveLongEdges();

  /*!
   * @brief Calculate edge lengths in a force directed layout (Fruchterman & Reingold 1991)
   * by applying (Barnes & Hut 1986) approximation (can be drawn with misc/plotter.py).
   */
  void CreateForceDirectedLayout(const std::string& path = "");

  /*!
   * @brief Employs Racon module on all eligible unitigs (Vaser et al 2017).
   */
  void Polish(const std::vector<std::unique_ptr<biosoup::Sequence>>& sequences,
      std::uint8_t match, std::uint8_t mismatch, std::uint8_t gap,
      std::uint32_t cuda_poa_batches, bool cuda_banded_alignment,
      std::uint32_t cuda_alignment_batches, std::uint32_t num_rounds);

  /*!
   * @brief Creates unitigs which are at least epsilon away from junction
   * nodes. This can be used to speed up creation of the force directed layout.
   */
  std::uint32_t CreateUnitigs(std::uint32_t epsilon = 0);

  void GetUnitigs(std::vector<std::unique_ptr<biosoup::Sequence>>& dst,
      bool drop_unpolished = false);

  /*!
   * @brief Prints all valid read piles in JSON format (can be drawn
   * with misc/plotter.py).
   */
  void PrintJson(const std::string& path) const;

  /*!
   * @brief Prints the assembly graph in csv format.
   */
  void PrintCsv(const std::string& path) const;

  /*!
   * @brief Prints the assembly graph in gfa format.
   */
  void PrintGfa(const std::string& path) const;

 private:
  void RemoveMarkedObjects(bool remove_nodes = false);

  void FindRemovableEdges(std::vector<std::uint32_t>& dst,
      const std::vector<std::uint32_t>& path);

  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  ram::MinimizerEngine minimizer_engine_;

  std::vector<std::unique_ptr<Pile>> piles_;

  struct Node;
  std::vector<std::unique_ptr<Node>> nodes_;

  struct Edge;
  std::vector<std::unique_ptr<Edge>> edges_;
  std::unordered_set<std::uint32_t> marked_edges_;
};

}  // namespace raven

#endif  // RAVEN_GRAPH_HPP_
