/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <unordered_set>

namespace thread_pool {
    class ThreadPool;
}

namespace ram {
    struct Sequence;
    class MinimizerEngine;
}

namespace raven {

class Pile;

class Graph;
std::unique_ptr<Graph> createGraph(std::shared_ptr<thread_pool::ThreadPool> thread_pool);

class Graph {
public:
    ~Graph();

    /*!
     * @brief Constructs the assembly graph by breaking chimeric sequences,
     * removing contained sequences and removing overlaps between repetitive
     * sequences.
     */
    void construct(std::vector<std::unique_ptr<ram::Sequence>>& sequences);

    /*!
     * @brief Simplify the assembly graph via transitive reduction, tip
     * pruning and popping bubble-like structures, and return contigs afterwards.
     */
    void assemble(std::vector<std::unique_ptr<ram::Sequence>>& dst);

    /*!
     * @brief Removes transitive edge (inspired by Myers 1995 & 2005).
     */
    std::uint32_t remove_transitive_edges();

    /*!
     * @brief Removes nodes which are dead ends in graph.
     */
    std::uint32_t remove_tips();

    /*!
     * @brief Removes bubble-like strucutres.
     */
    std::uint32_t remove_bubbles();

    /*!
     * @brief Removes long edges in the force directed layout which function
     * needs to be called before.
     */
    std::uint32_t remove_long_edges();

    /*!
     * @brief Calculate edge lengths in a force directed layout (Fruchterman & Reingold 1991)
     * (can be drawn with misc/plotter.py).
     */
    void create_force_directed_layout(const std::string& path = "");

    /*!
     * @brief Creates unitigs which are at least epsilon away from junction
     * nodes. This can be used to speed up creation of the force directed layout.
     */
    std::uint32_t create_unitigs(std::uint32_t epsilon);

    /*!
     * @brief Creates unitigs by merging chains of overlapping sequences.
     */
    std::uint32_t create_unitigs();

    void extract_unitigs(std::vector<std::unique_ptr<ram::Sequence>>& dst);

    /*!
     * @brief Prints all valid read piles in JSON format (can be drawn
     * with misc/plotter.py).
     */
    void print_json(const std::string& path) const;

    /*!
     * @brief Prints the assembly graph in csv format.
     */
    void print_csv(const std::string& path) const;

    friend std::unique_ptr<Graph> createGraph(
        std::shared_ptr<thread_pool::ThreadPool> thread_pool);
private:
    Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    void remove_marked_objects(bool remove_nodes = false);

    void find_removable_edges(std::vector<std::uint32_t>& dst,
        const std::vector<std::uint32_t>& path);

    std::unique_ptr<ram::MinimizerEngine> minimizer_engine_;

    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

    std::vector<std::unique_ptr<Pile>> piles_;

    struct Node;
    std::vector<std::unique_ptr<Node>> nodes_;

    struct Edge;
    std::vector<std::unique_ptr<Edge>> edges_;
    std::unordered_set<std::uint32_t> marked_edges_;

    std::vector<std::pair<std::uint32_t, std::uint32_t>> transitive_edges_;
};

}
