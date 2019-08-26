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

namespace thread_pool {
    class ThreadPool;
}

namespace ram {
    class Sequence;
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
     * sequences
     */
    void construct(const std::vector<std::unique_ptr<ram::Sequence>>& sequences);

    /*!
    * @brief Prints all unresolved graph junctions in JSON format (can be drawn
    * with misc/plotter.py)
    */
    void print_json(const std::string& path) const;

    friend std::unique_ptr<Graph> createGraph(
        std::shared_ptr<thread_pool::ThreadPool> thread_pool);
private:
    Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
    std::vector<std::unique_ptr<Pile>> piles_;
};

}
