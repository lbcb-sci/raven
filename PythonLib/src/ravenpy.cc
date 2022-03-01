#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <utility>

// pybind11 includes
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

// thread pool and bioparser
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "thread_pool/thread_pool.hpp"
// raven includes
#include "raven/graph/assemble.h"
#include "raven/graph/common.h"
#include "raven/graph/construct.h"
#include "raven/graph/graph.h"
#include "raven/graph/polish.hpp"
#include "raven/graph/serialization/graph_repr.h"
#include "raven/io.h"

namespace py = pybind11;

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

static std::vector<std::unique_ptr<biosoup::NucleicAcid>> LoadSequences(
    const std::vector<std::string>& paths) {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto parse_file = [](const std::string& path)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
    auto sparser = raven::CreateParser(path);
    if (!sparser) {
      throw std::runtime_error(
          "[ravenpy::LoadSequences] failed to create parser for : " + path);
    }

    return sparser->Parse(-1);
  };

  for (const auto& path : paths) {
    auto chunk = parse_file(path);
    if (dst.empty()) {
      std::swap(dst, chunk);
    } else {
      std::move(chunk.begin(), chunk.end(), std::back_inserter(dst));
    }
  }

  dst.shrink_to_fit();
  return dst;
}

struct SequencesHandle {
  SequencesHandle(const std::vector<std::string>& paths)
      : seqs(LoadSequences(paths)) {}

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> seqs;
};

PYBIND11_MODULE(ravenpy, m) {
  m.doc() = "Raven";

  py::class_<raven::AlignCfg>(m, "AlignCfg")
      .def(py::init<std::int8_t, std::int8_t, std::int8_t>())
      .def_readwrite("match", &raven::AlignCfg::match)
      .def_readwrite("mismatch", &raven::AlignCfg::mismatch)
      .def_readwrite("gap", &raven::AlignCfg::gap);

  py::class_<raven::CudaCfg>(m, "CudaCfg")
      .def(py::init<std::uint32_t, std::uint32_t, bool>())
      .def_readwrite("poa_batches", &raven::CudaCfg::poa_batches)
      .def_readwrite("alignment_batches", &raven::CudaCfg::alignment_batches)
      .def_readwrite("banded_alignment", &raven::CudaCfg::banded_alignment);

  py::class_<raven::PolishCfg>(m, "PolishCfg")
      .def(py::init<raven::AlignCfg, raven::CudaCfg, std::uint32_t>())
      .def_readwrite("align_cfg", &raven::PolishCfg::align_cfg)
      .def_readwrite("cuda_cfg", &raven::PolishCfg::cuda_cfg)
      .def_readwrite("num_rounds", &raven::PolishCfg::num_rounds);

  py::class_<raven::OverlapPhaseCfg>(m, "OverlapPhaseCfg")
      .def(py::init<std::uint8_t, std::uint8_t, double>())
      .def_readwrite("kmer_len", &raven::OverlapPhaseCfg::kmer_len)
      .def_readwrite("window_len", &raven::OverlapPhaseCfg::window_len)
      .def_readwrite("freq", &raven::OverlapPhaseCfg::freq);

  py::class_<thread_pool::ThreadPool, std::shared_ptr<thread_pool::ThreadPool>>(
      m, "ThreadPool")
      .def(py::init<std::uint32_t>());

  py::class_<SequencesHandle, std::shared_ptr<SequencesHandle>>(
      m, "SequencesHandle")
      .def(py::init<const std::vector<std::string>&>());

  py::class_<raven::Graph>(m, "Graph").def(py::init<>());

  m.def("construct_graph",
        [](raven::Graph& graph, std::shared_ptr<SequencesHandle> seqs_handle,
           std::shared_ptr<thread_pool::ThreadPool> thread_pool,
           bool checkpoints, raven::OverlapPhaseCfg ovlp_cfg) -> void {
          raven::ConstructGraph(graph, seqs_handle->seqs, thread_pool,
                                checkpoints, ovlp_cfg);
        });

  m.def("assemble_graph", raven::Assemble);

  m.def("polish_graph",
        [](std::shared_ptr<thread_pool::ThreadPool> thread_pool,
           raven::Graph& graph, bool checkpoints,
           std::shared_ptr<SequencesHandle> seqs_handle,
           raven::PolishCfg cfg) -> void {
          raven::Polish(thread_pool, graph, checkpoints, seqs_handle->seqs,
                        cfg);
        });

  m.def("graph_print_gfa", raven::PrintGfa);
  m.def("graph_print_unitgs",
        [](raven::Graph& graph,
           const std::uint32_t num_polishing_rounds) -> void {
          for (const auto& it :
               raven::GetUnitigs(graph, num_polishing_rounds > 0)) {
            std::cout << ">" << it->name << std::endl;
            std::cout << it->InflateData() << std::endl;
          }
        });
}
