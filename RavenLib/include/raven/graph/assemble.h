#ifndef RAVEN_GRAPH_ASSEMBLE_H_
#define RAVEN_GRAPH_ASSEMBLE_H_

#include "raven/graph/graph.h"
#include "thread_pool/thread_pool.hpp"

namespace raven {

// Layout phase
// - remove transitive edges
// - remove tips
// - remove bubbles
// - remove elongated edges in 2D layout
void Assemble(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
              Graph& graph, bool checkpoints);

}  // namespace raven

#endif  // RAVEN_GRAPH_ASSEMBLE_H_