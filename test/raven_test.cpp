// Copyright (c) 2020 Robert Vaser

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "edlib.h"  // NOLINT
#include "gtest/gtest.h"

#include "graph.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace raven {
namespace test {

class RavenTest: public ::testing::Test {
 public:
  void SetUp() {
    biosoup::Sequence::num_objects = 0;
    auto p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(  // NOLINT
        RAVEN_DATA_PATH + std::string("ERA476754.fastq.gz"));
    s = p->Parse(-1);

    biosoup::Sequence::num_objects = 0;
    p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(  // NOLINT
        RAVEN_DATA_PATH + std::string("NC_001416.fasta.gz"));
    r = std::move(p->Parse(-1).front());
  }

  static std::uint32_t EditDistance(
      const std::string& lhs,
      const std::string &rhs) {
    EdlibAlignResult result = edlibAlign(
        lhs.c_str(), lhs.size(),
        rhs.c_str(), rhs.size(),
        edlibDefaultAlignConfig());
    std::uint32_t ed = result.editDistance;
    edlibFreeAlignResult(result);
    return ed;
  }

  std::vector<std::unique_ptr<biosoup::Sequence>> s;
  std::unique_ptr<biosoup::Sequence> r;
};

TEST_F(RavenTest, Assemble) {
  Graph g{false, false, nullptr};
  g.Construct(s);
  g.Assemble();
  g.Polish(s, 3, -5, -4, 0, false, 0, 2);
  auto u = std::move(g.GetUnitigs(true).front());
  u->ReverseAndComplement();
  EXPECT_EQ(1128, EditDistance(u->data, r->data));
}

TEST_F(RavenTest, Checkpoints) {
  Graph g{false, false, nullptr};
  g.Construct(s);
  g.Assemble();
  g.Polish(s, 2, -5, -2, 0, false, 0, 2);
  auto u = std::move(g.GetUnitigs(true).front());

  SetUp();

  g = {false, true, nullptr};
  g.Construct(s);

  g = {false, true, nullptr};
  g.Load();
  g.Assemble();

  g = {false, true, nullptr};
  g.Load();
  g.Polish(s, 2, -5, -2, 0, false, 0, 1);

  g = {false, true, nullptr};
  g.Load();
  g.Polish(s, 2, -5, -2, 0, false, 0, 2);

  EXPECT_EQ(u->data, g.GetUnitigs(true).front()->data);
}

}  // namespace test
}  // namespace raven
