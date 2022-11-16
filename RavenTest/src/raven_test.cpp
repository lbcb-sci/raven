// Copyright (c) 2020 Robert Vaser

#include <edlib.h>  // NOLINT
#include <gtest/gtest.h>

#include <bioparser/fasta_parser.hpp>
#include <bioparser/fastq_parser.hpp>

#include "raven/graph/assemble.h"
#include "raven/graph/common.h"
#include "raven/graph/construct.h"
#include "raven/graph/graph.h"
#include "raven/graph/polish.hpp"
#include "raven/graph/serialization/binary.h"
#include "raven/graph/serialization/graph_repr.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace raven::test {

class RavenTest : public ::testing::Test {
 public:
  void SetUp() {
    biosoup::NucleicAcid::num_objects = 0;
    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastqParser>(  // NOLINT
        TEST_DATA + std::string("ERA476754.fastq.gz"));
    s = p->Parse(-1);

    biosoup::NucleicAcid::num_objects = 0;
    p = bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastaParser>(  // NOLINT
        TEST_DATA + std::string("NC_001416.fasta.gz"));
    r = std::move(p->Parse(-1).front());
  }

  static std::uint32_t EditDistance(const std::string& lhs,
                                    const std::string& rhs) {
    EdlibAlignResult result = edlibAlign(lhs.c_str(), lhs.size(), rhs.c_str(),
                                         rhs.size(), edlibDefaultAlignConfig());
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

  raven::OverlapPhaseCfg overlapPhaseCfg = OverlapPhaseCfg{};
  overlapPhaseCfg.useMinhash = true;

  ConstructGraph(graph, s, threadPool, false, overlapPhaseCfg);
  Assemble(threadPool, graph, false);
  Polish(threadPool, graph, false, s, PolishCfg{});

  auto u = std::move(GetUnitigs(graph).front());

  u->ReverseAndComplement();

  EXPECT_EQ(1137, EditDistance(u->InflateData(), r->InflateData()));
}

TEST_F(RavenTest, Checkpoints) {
  Graph graph;

  auto threadPool = std::make_shared<thread_pool::ThreadPool>(1);

  ConstructGraph(graph, s, threadPool, false, OverlapPhaseCfg{});
  Assemble(threadPool, graph, false);
  Polish(threadPool, graph, false, s, PolishCfg{});

  auto u = std::move(GetUnitigs(graph).front());

  SetUp();

  graph = {};
  raven::ConstructGraph(graph, s, threadPool, true, OverlapPhaseCfg{});

  graph = raven::LoadGraphFromFile();
  raven::Assemble(threadPool, graph, true);

  graph = raven::LoadGraphFromFile();
  raven::Polish(threadPool, graph, true, s, raven::PolishCfg{});

  graph = raven::LoadGraphFromFile();
  raven::Polish(threadPool, graph, true, s, raven::PolishCfg{});

  EXPECT_EQ(u->deflated_data, GetUnitigs(graph, true).front()->deflated_data);
}

}  // namespace raven::test
