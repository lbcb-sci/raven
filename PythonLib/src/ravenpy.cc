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
using namespace pybind11::literals;

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

struct OverlapsHandle {
  OverlapsHandle(std::shared_ptr<SequencesHandle> seqs_handle)
      : overlaps() { overlaps.resize(seqs_handle->seqs.size()); }

  std::vector<std::vector<biosoup::Overlap>> overlaps;
};

PYBIND11_MODULE(ravenpy, m) {
  m.doc() = "raven";

  py::class_<biosoup::NucleicAcid>(m, "NucleicAcid")
      .def(py::init<const std::string&, const std::string&>())
      .def(py::init<const std::string&, const std::string&,
                    const std::string&>())
      .def(
          "__deepcopy__",
          [](const biosoup::NucleicAcid& seq,
             py::dict) -> biosoup::NucleicAcid { return seq; },
          "memo"_a)
      .def("code", &biosoup::NucleicAcid::Code)
      .def("score", &biosoup::NucleicAcid::Score)
      .def("inflate_data", &biosoup::NucleicAcid::InflateData,
           py::arg("i") = 0U,
           py::arg("len") = std::numeric_limits<std::uint32_t>::max())
      .def("inflate_quality", &biosoup::NucleicAcid::InflateQuality,
           py::arg("i") = 0U,
           py::arg("len") = std::numeric_limits<std::uint32_t>::max())
      .def("reverse_and_complement",
           &biosoup::NucleicAcid::ReverseAndComplement)
      .def_readwrite("id", &biosoup::NucleicAcid::id)
      .def_readwrite("name", &biosoup::NucleicAcid::name)
      .def_readonly("__len__", &biosoup::NucleicAcid::inflated_len);

  py::class_<biosoup::Overlap>(m, "Overlap")
      .def(py::init<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
                    std::uint32_t, std::uint32_t, std::uint32_t, bool>(),
           py::arg("lhs_id"), py::arg("lhs_begin"), py::arg("lhs_end"),
           py::arg("rhs_id"), py::arg("rhs_begin"), py::arg("rhs_end"),
           py::arg("score"), py::arg("strand") = true)
      .def(py::init<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
                    std::uint32_t, std::uint32_t, std::uint32_t,
                    const std::string&, bool>(),
           py::arg("lhs_id"), py::arg("lhs_begin"), py::arg("lhs_end"),
           py::arg("rhs_id"), py::arg("rhs_begin"), py::arg("rhs_end"),
           py::arg("score"), py::arg("alignment"), py::arg("strand") = true)
      .def(
          "__deepcopy__",
          [](const biosoup::Overlap& ovlp, py::dict) -> biosoup::Overlap {
            return ovlp;
          },
          "memo"_a)
      .def_readwrite("lhs_id", &biosoup::Overlap::lhs_id)
      .def_readwrite("lhs_begin", &biosoup::Overlap::lhs_begin)
      .def_readwrite("lhs_end", &biosoup::Overlap::lhs_end)
      .def_readwrite("rhs_id", &biosoup::Overlap::rhs_id)
      .def_readwrite("rhs_begin", &biosoup::Overlap::rhs_begin)
      .def_readwrite("rhs_end", &biosoup::Overlap::rhs_end)
      .def_readwrite("score", &biosoup::Overlap::score)
      .def_readwrite("strand", &biosoup::Overlap::strand)
      .def_readwrite("alignment", &biosoup::Overlap::alignment);


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
      .def(py::init<std::uint8_t, std::uint8_t, double, double, std::size_t, double>())
      .def_readwrite("kmer_len", &raven::OverlapPhaseCfg::kmer_len)
      .def_readwrite("window_len", &raven::OverlapPhaseCfg::window_len)
      .def_readwrite("freq", &raven::OverlapPhaseCfg::freq)
      .def_readwrite("identity", &raven::OverlapPhaseCfg::identity)
      .def_readwrite("kMaxNumOverlaps", &raven::OverlapPhaseCfg::kMaxNumOverlaps)
      .def_readwrite("weightedMinimizerSampling", &raven::OverlapPhaseCfg::weightedMinimizerSampling);

  py::class_<thread_pool::ThreadPool, std::shared_ptr<thread_pool::ThreadPool>>(
      m, "ThreadPool")
      .def(py::init<std::uint32_t>());

  py::class_<raven::Graph>(m, "Graph").def(py::init<>());

  py::class_<SequencesHandle, std::shared_ptr<SequencesHandle>>(
      m, "SequencesHandle")
      .def(py::init<const std::vector<std::string>&>());

  py::class_<OverlapsHandle, std::shared_ptr<OverlapsHandle>>(
      m, "OverlapsHandle")
      .def(py::init<std::shared_ptr<SequencesHandle>>())
      .def_readonly("overlaps", &OverlapsHandle::overlaps);

  m.def("reset_seq_id_cnt",
        []() -> void { biosoup::NucleicAcid::num_objects.store(0); });

  m.def("load_sequences",
        [](std::vector<std::string>& paths) -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
          return LoadSequences(paths);
        });

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

  m.def("graph_print_csv", raven::PrintCsv);
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

  m.def("graph_get_csv", raven::getCsv);
  m.def("graph_get_gfa", raven::getGfa);

  m.def("graph_load_gfa", raven::LoadGfa);

  py::class_<ram::MinimizerEngine>(m, "MinimizerEngine")
      .def(py::init<std::shared_ptr<thread_pool::ThreadPool>, std::uint32_t, std::uint32_t>());
 
  m.def("find_overlaps_and_create_piles",
      [](std::shared_ptr<thread_pool::ThreadPool> thread_pool, ram::MinimizerEngine& minimizer_engine, std::shared_ptr<SequencesHandle> seqs_handle,
          double freq, raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle, std::size_t kMaxNumOverlaps) -> void {
        raven::FindOverlapsAndCreatePiles(thread_pool, minimizer_engine, seqs_handle->seqs, freq, graph.piles, overlaps_handle->overlaps, kMaxNumOverlaps);
      });

  m.def("trim_and_annotate_piles",
      [](std::shared_ptr<thread_pool::ThreadPool> thread_pool,raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle) -> void {
        raven::TrimAndAnnotatePiles(thread_pool, graph.piles, overlaps_handle->overlaps);
      });

  m.def("resolve_contained_reads",
      [](raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle, std::shared_ptr<SequencesHandle> seqs_handle,
          std::shared_ptr<thread_pool::ThreadPool> thread_pool, double identity) -> void {
        raven::ResolveContainedReads(graph.piles, overlaps_handle->overlaps, seqs_handle->seqs, thread_pool,  identity);
      });

  m.def("resolve_chimeric_sequences",
      [](std::shared_ptr<thread_pool::ThreadPool> thread_pool,raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle,
          std::shared_ptr<SequencesHandle> seqs_handle) -> void {
        raven::ResolveChimericSequences(thread_pool, graph.piles, overlaps_handle->overlaps, seqs_handle->seqs);
      });

  m.def("find_overlaps_and_repetetive_regions",
      [](std::shared_ptr<thread_pool::ThreadPool> thread_pool, ram::MinimizerEngine& minimizer_engine, double freq, std::uint8_t kmer_len,
          double identity, raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle, std::shared_ptr<SequencesHandle> seqs_handle) -> void {
        raven::FindOverlapsAndRepetetiveRegions(thread_pool, minimizer_engine, freq, kmer_len, identity, graph.piles, overlaps_handle->overlaps, seqs_handle->seqs);
      });

  m.def("resolve_repeat_induced_overlaps",
      [](std::shared_ptr<thread_pool::ThreadPool> thread_pool, raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle,
          std::shared_ptr<SequencesHandle> seqs_handle) -> void {
        raven::ResolveRepeatInducedOverlaps(thread_pool,  graph.piles, overlaps_handle->overlaps, seqs_handle->seqs);
      });

  m.def("construct_assembly_graph",
      [](raven::Graph& graph, std::shared_ptr<OverlapsHandle> overlaps_handle,
          std::shared_ptr<SequencesHandle> seqs_handle) -> void {
        raven::ConstructAssemblyGraph(graph, graph.piles, overlaps_handle->overlaps, seqs_handle->seqs);
      });

  m.def("remove_transitive_edges_from_graph",
      [](raven::Graph& graph) -> std::uint32_t {
        return raven::RemoveTransitiveEdgesFromGraph(graph);
      });

  m.def("remove_tips_and_bubbles_from_graph",
      [](raven::Graph& graph) -> void {
        raven::RemoveTipsAndBubblesFromGraph(graph);
      });

  m.def("remove_long_edges_from_graph",
      [](raven::Graph& graph, std::shared_ptr<thread_pool::ThreadPool> thread_pool) -> void {
        raven::RemoveLongEdgesFromGraph(graph, thread_pool);
      });

}
