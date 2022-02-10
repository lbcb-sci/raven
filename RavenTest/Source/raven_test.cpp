// Copyright (c) 2020 Robert Vaser

#include <gtest/gtest.h>

#include <bioparser/fasta_parser.hpp>
#include <bioparser/fastq_parser.hpp>
#include <edlib.h>  // NOLINT

#include <Graph/Graph.hpp>
#include <Graph/GraphConstruct.hpp>
#include <Graph/GraphAssemble.hpp>
#include <Graph/GraphPolish.hpp>
#include <Graph/GraphShared.hpp>
#include <Graph/Serialization/GraphBinarySerialization.hpp>
#include <Graph/Serialization/GraphOutputs.hpp>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace raven {

namespace test {

class RavenTest : public ::testing::Test {
public:
    void SetUp() {
        biosoup::NucleicAcid::num_objects = 0;
        auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(  // NOLINT
                TEST_DATA + std::string("ERA476754.fastq.gz"));
        s = p->Parse(-1);

        biosoup::NucleicAcid::num_objects = 0;
        p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(  // NOLINT
                TEST_DATA + std::string("NC_001416.fasta.gz"));
        r = std::move(p->Parse(-1).front());
    }

    static std::uint32_t EditDistance(
            const std::string& lhs,
            const std::string& rhs) {
        EdlibAlignResult result = edlibAlign(
                lhs.c_str(), lhs.size(),
                rhs.c_str(), rhs.size(),
                edlibDefaultAlignConfig());
        std::uint32_t ed = result.editDistance;
        edlibFreeAlignResult(result);
        return ed;
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> s;
    std::unique_ptr<biosoup::NucleicAcid> r;
};

TEST_F(RavenTest, Assemble) {
    Graph graph;

    auto threadPool = std::make_shared<thread_pool::ThreadPool>(1);

    raven::constructGraph(graph, s, threadPool, false);
    raven::assemble(threadPool, graph, false);
    raven::polish(threadPool, graph, false, s, 3, -5, -4, 0, false, 0, 2);

    auto u = std::move(getUnitigs(graph).front());

    u->ReverseAndComplement();

    EXPECT_EQ(1137, EditDistance(u->InflateData(), r->InflateData()));
}

TEST_F(RavenTest, Checkpoints) {
    Graph graph;

    auto threadPool = std::make_shared<thread_pool::ThreadPool>(1);

    raven::constructGraph(graph, s, threadPool, false);
    raven::assemble(threadPool, graph, false);
    raven::polish(threadPool, graph, false, s, 2, -5, -2, 0, false, 0, 2);

    auto u = std::move(getUnitigs(graph).front());

    SetUp();

    graph = {};
    raven::constructGraph(graph, s, threadPool, true);

    graph = raven::loadGraphFromFile();
    raven::assemble(threadPool, graph, true);


    graph = raven::loadGraphFromFile();
    raven::polish(threadPool, graph, true, s, 2, -5, -2, 0, false, 0, 1);

    graph = raven::loadGraphFromFile();
    raven::polish(threadPool, graph, true, s, 2, -5, -2, 0, false, 0, 2);


    EXPECT_EQ(u->deflated_data, getUnitigs(graph, true).front()->deflated_data);
}

}  // namespace test
}  // namespace raven
