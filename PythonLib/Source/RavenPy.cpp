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

std::uint8_t kmerLength                 = 10;
std::uint8_t windowLength               = 5;
double       frequency                  = 0.001;

std::int32_t numberOfPolishingRounds    =  2;
std::int8_t  polishingMatch             =  3;
std::int8_t  polishingMismatch          = -5;
std::int8_t  polishingGap               = -4;

std::uint32_t cudaPoaBatches            = 0;
std::uint32_t cudaAlignmentBatches      = 0;
bool          cudaBandedAlignment       = false;


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


static std::uint8_t getKmerLength() {
    return kmerLength;
}

static std::uint8_t getWindowLength() {
    return windowLength;
}

static double getFrequency() {
    return frequency;
}

static std::int32_t getNumberOfPolishingRounds() {
    return numberOfPolishingRounds;
}

static std::int8_t getPolishingMatch() {
    return polishingMatch;
}

static std::int8_t getPolishingMismatch() {
    return polishingMismatch;
}

static std::int8_t getPolishingGap() {
    return polishingGap;
}

static std::uint32_t getCudaPoaBatches() {
    return cudaPoaBatches;
}

static std::uint32_t getCudaAlignmentBatches() {
    return cudaAlignmentBatches;
}

static bool getCudaBandedAlignment() {
    return cudaBandedAlignment;
}

static void setKmerLength(std::uint8_t newKmerLength) {
    kmerLength = newKmerLength;
}

static void setWindowLength(std::uint8_t newWindowLength) {
    windowLength = newWindowLength;
}

static void setFrequency(double newFrequency) {
    frequency = newFrequency;
}

static void setNumberOfPolishingRounds(std::int32_t newNumberOfPolishingRounds) {
    numberOfPolishingRounds = newNumberOfPolishingRounds;
}

static void setPolishingMatch(std::int8_t newPolishingMatch) {
    polishingMatch = newPolishingMatch;
}

static void setPolishingMismatch(std::int8_t newPolishingMismatch) {
    polishingMismatch = newPolishingMismatch;
}

static void setPolishingGap(std::int8_t newPolishingGap) {
    polishingGap = newPolishingGap;
}

static void setCudaPoaBatches(std::uint32_t newCudaPoaBatches) {
    cudaPoaBatches = newCudaPoaBatches;
}

static void setCudaAlignmentBatches(std::uint32_t newCudaAlignmentBatches) {
    cudaAlignmentBatches = newCudaAlignmentBatches;
}

static void setCudaBandedAlignment(bool newCudaBandedAlignment) {
    cudaBandedAlignment = newCudaBandedAlignment;
}

static void constructGraph(Raven& raven, bool checkpoints) {
    raven::ConstructGraph(raven.graph, raven.sequences, raven.threadPool, checkpoints, kmerLength, windowLength, frequency);
}

static void assemble(Raven& raven, bool checkpoints) {
    raven::Assemble(raven.threadPool, raven.graph, checkpoints);
}

static void polish(Raven& raven, bool checkpoints) {

    raven::Polish(raven.threadPool, raven.graph, checkpoints, raven.sequences, polishingMatch, polishingMismatch, polishingGap, cudaPoaBatches, cudaBandedAlignment, cudaAlignmentBatches, numberOfPolishingRounds);

}

PYBIND11_MODULE(RavenPy, m) {

    m.doc() = "Raven";

    m.def("initializeRaven", &initializeRaven, py::return_value_policy::take_ownership);

    m.def("getKmerLength", &getKmerLength, "Get k-mer length");
    m.def("getWindowLength", &getWindowLength, "Get window length");
    m.def("getFrequency", &getFrequency, "Get frequency");
    m.def("getNumberOfPolishingRounds", &getNumberOfPolishingRounds, "Get number of polishing rounds");
    m.def("getPolishingMatch", &getPolishingMatch, "Get match for polishing");
    m.def("getPolishingMismatch", &getPolishingMismatch, "Get mismatch for polishing");
    m.def("getPolishingGap", &getPolishingGap, "Get gap for polishing");
    m.def("getCudaPoaBatches", &getCudaPoaBatches, "Get cuda poa batches");
    m.def("getCudaAlignmentBatches", &getCudaAlignmentBatches, "Get cuda alignment batches");
    m.def("getCudaBandedAlignment", &getCudaBandedAlignment, "Get set cuda banded alignment");

    m.def("setKmerLength", &setKmerLength, "Set k-mer length");
    m.def("setWindowLength", &setWindowLength, "Set window length");
    m.def("setFrequency", &setFrequency, "Set frequency");
    m.def("setNumberOfPolishingRounds", &setNumberOfPolishingRounds, "Set number of polishing rounds");
    m.def("setPolishingMatch", &setPolishingMatch, "Set match for polishing");
    m.def("setPolishingMismatch", &setPolishingMismatch, "Set mismatch for polishing");
    m.def("setPolishingGap", &setPolishingGap, "Set gap for polishing");
    m.def("setCudaPoaBatches", &setCudaPoaBatches, "Set cuda poa batches");
    m.def("setCudaAlignmentBatches", &setCudaAlignmentBatches, "Set cuda alignment batches");
    m.def("setCudaBandedAlignment", &setCudaBandedAlignment, "Set set cuda banded alignment");

    py::class_<Raven, std::shared_ptr<Raven>>(m, "Raven")
        .def("constructGraph", &constructGraph)
        .def("assemble", &assemble)
        .def("polish", &polish);
}
