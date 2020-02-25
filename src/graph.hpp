/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <cstdint>
#include <string>
#include <memory>
#include <vector>
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
     * pruning and popping bubble-like structures.
     */
    void assemble();

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
     * by applying (Barnes & Hut 1986) approximation (can be drawn with misc/plotter.py).
     */
    void create_force_directed_layout(const std::string& path = "");

    /*!
     * @brief Employs Racon module on all eligible unitigs (Vaser et al 2017).
     */
    void polish(const std::vector<std::unique_ptr<ram::Sequence>>& sequences,
        std::uint8_t match, std::uint8_t mismatch, std::uint8_t gap,
        std::uint32_t cuda_poa_batches, bool cuda_banded_alignment,
        std::uint32_t cuda_alignment_batches, std::uint32_t num_rounds);

    /*!
     * @brief Creates unitigs which are at least epsilon away from junction
     * nodes. This can be used to speed up creation of the force directed layout.
     */
    std::uint32_t create_unitigs(std::uint32_t epsilon = 0);

    void get_unitigs(std::vector<std::unique_ptr<ram::Sequence>>& dst,
        bool drop_unpolished = false);

    /*!
     * @brief Prints all valid read piles in JSON format (can be drawn
     * with misc/plotter.py).
     */
    void print_json(const std::string& path) const;

    /*!
     * @brief Prints the assembly graph in csv format.
     */
    void print_csv(const std::string& path) const;

    /*!
     * @brief Prints the assembly graph in gfa format.
     */
    void print_gfa(const std::string& path) const;

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
};

}
