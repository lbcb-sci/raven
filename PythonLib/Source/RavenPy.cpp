#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <utility>
#include <memory>

#include <thread_pool/thread_pool.hpp>
#include <bioparser/fasta_parser.hpp>
#include <bioparser/fastq_parser.hpp>

#include <Graph/Graph.hpp>
#include <Graph/GraphConstruct.hpp>
#include <Graph/GraphAssemble.hpp>
#include <Graph/GraphPolish.hpp>

namespace py = pybind11;

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

struct Raven {
    raven::Graph graph;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    std::shared_ptr<thread_pool::ThreadPool> threadPool;
};

static std::shared_ptr<Raven> initializeRaven(std::string sequenceFile, int numberOfThreads) {

    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(sequenceFile);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences = p->Parse(-1);

    return std::shared_ptr<Raven>( new Raven { {},  std::move(sequences), std::make_shared<thread_pool::ThreadPool>(numberOfThreads)});
}

static void constructGraph(Raven& raven, bool checkpoints) {
    raven::constructGraph(raven.graph, raven.sequences, raven.threadPool, checkpoints);
}

static void assemble(Raven& raven, bool checkpoints) {
    raven::assemble(raven.threadPool, raven.graph, checkpoints);
}

static void polish(Raven& raven, bool checkpoints) {
    raven::polish(raven.threadPool, raven.graph, checkpoints, raven.sequences, 3, -5, -4, 0, false, 0, 2);
}

PYBIND11_MODULE(RavenPy, m) {

    m.doc() = "Raven";

    m.def("initializeRaven", &initializeRaven, py::return_value_policy::take_ownership);

    py::class_<Raven, std::shared_ptr<Raven>>(m, "Raven")
        .def("constructGraph", &constructGraph)
        .def("assemble", &assemble)
        .def("polish", &polish);
}
