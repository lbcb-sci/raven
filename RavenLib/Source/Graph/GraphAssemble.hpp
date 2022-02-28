#pragma once

#include "Graph.hpp"

namespace raven {

// Layout phase
// - remove transitive edges
// - remove tips
// - remove bubbles
// - remove elongated edges in 2D layout
void assemble(std::shared_ptr<thread_pool::ThreadPool> threadPool, Graph& graph,
              bool checkpoints);

}  // namespace raven
