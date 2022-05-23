// Copyright (c) 2020 Robert Vaser

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"
#include "raven/raven.h"
#include "raven_cfg.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

int main(int argc, char** argv) {
  auto options = cxxopts::Options(
      "raven",
      "Raven is a de novo genome assembler for long uncorrected reads");

  /* clang-format off */
  options.add_options("maping arguments")
    ("k,kmer-len", "length of minimizers used to find overlaps",
      cxxopts::value<std::uint8_t>()->default_value("15"))
    ("w,window-len", "length of sliding window from which minimizers are sampled",
      cxxopts::value<std::uint8_t>()->default_value("5"))
    ("f,frequency", "threshold for ignoring most frequent minimizers",
      cxxopts::value<double>()->default_value("0.001"));
  options.add_options("layout arguments")
    ("i,identity", 
     "threshold for overlap between two reads in order to construct an edge between them", 
      cxxopts::value<double>()->default_value("0"))
    ("o,kMaxNumOverlaps", 
      "maximum number of overlaps that will be taken during FindOverlapsAndCreatePiles stage", 
      cxxopts::value<std::size_t>()->default_value("32"));
  options.add_options("polishing arguments")
    ("p,polishing-rounds", "number of times racon is invoked",
      cxxopts::value<std::uint32_t>()->default_value("2"))
    ("m,match", "score for matching bases",
      cxxopts::value<std::int8_t>()->default_value("3"))
    ("n,mismatch", "score for mismatching bases",
      cxxopts::value<std::int8_t>()->default_value("-5"))
    ("g,gap", "gap penalty (must be negative)",
      cxxopts::value<std::int8_t>()->default_value("-4"));

#ifndef CUDA_ENABLED
  options.add_options("cuda arguments")
    ("c,cuda-poa-batches", "number of batches for CUDA accelerated polishing",
      cxxopts::value<std::uint32_t>()->default_value("0"))
    ("b,cuda-banded-alignment", 
     "use banding approximation for polishing on GPU; (only applicable when -c is used)",
     cxxopts::value<bool>()->default_value("false"))
    ("a,cuda-alignment-batches", "number of batches for CUDA accelerated alignment",
      cxxopts::value<std::uint32_t>()->default_value("0"));
#endif /* CUDA_ENABLED */

  options.add_options("utility arguments")
    ("graphical-fragment-assembly", "prints the assembly graph in GFA format",
      cxxopts::value<std::string>()->default_value("./raven.gfa"))
    ("disable-checkpoints", "disable checkpoint file creation",
      cxxopts::value<bool>()->default_value("false"))
    ("resume", "resume previous run from last checkpoint",
      cxxopts::value<bool>()->default_value("false"))
    ("version", "prints the version number",
      cxxopts::value<bool>()->default_value("false"))
    ("t,threads", "number of threads", 
      cxxopts::value<std::uint32_t>()->default_value("1"))
    ("h,help", "display help message");

  options.add_options("sequence input")
    ("paths", "input fasta/fastq reads", 
      cxxopts::value<std::vector<std::string>>());
  /* clang-format on */

  options.parse_positional("paths");
  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cerr << options.help() << std::endl;
    return EXIT_SUCCESS;
  }

  if (result.count("version")) {
    std::cerr << PROJECT_VER << std::endl;
    return EXIT_SUCCESS;
  }

  auto const kmer_len = result["kmer-len"].as<std::uint8_t>();
  auto const window_len = result["window-len"].as<std::uint8_t>();
  auto const freq = result["frequency"].as<double>();
  auto const identity = result["identity"].as<double>();
  auto const kMaxNumOverlaps = result["kMaxNumOverlaps"].as<std::size_t>();

  auto const num_polishing_rounds =
      result["polishing-rounds"].as<std::uint32_t>();
  auto const m = result["match"].as<std::int8_t>();
  auto const n = result["mismatch"].as<std::int8_t>();
  auto const g = result["gap"].as<std::int8_t>();

  auto const gfa_path = result["graphical-fragment-assembly"].as<std::string>();
  auto const resume = result["resume"].as<bool>();
  auto const checkpoints = !result["disable-checkpoints"].as<bool>();

  auto const num_threads = result["threads"].as<std::uint32_t>();

#ifdef CUDA_ENABLED
  auto const cuda_poa_batches = 0U;
  auto const cuda_alignment_batches = 0U;
  auto const cuda_banded_alignment = false;
#else
  auto const cuda_poa_batches = result["cuda-poa-batches"].as<std::uint32_t>();
  auto const cuda_alignment_batches =
      result["cuda-alignment-batches"].as<std::uint32_t>();
  auto const cuda_banded_alignment = result["cuda-banded-alignment"].as<bool>();
#endif /* CUDA_ENABLED */

  biosoup::Timer timer{};
  timer.Start();

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  raven::Graph graph;

  if (resume) {
    try {
      graph = raven::LoadGraphFromFile();
    } catch (std::exception& exception) {
      std::cerr << exception.what() << std::endl;
      return EXIT_FAILURE;
    }

    std::cerr << "[raven::] loaded previous run " << std::fixed << timer.Stop()
              << "s" << std::endl;

    timer.Start();
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
  if (graph.stage < -3 || num_polishing_rounds > std::max(0, graph.stage)) {
    auto paths_strs = result["paths"].as<std::vector<std::string>>();
    auto paths = std::vector<std::filesystem::path>();

    for (auto const& it : paths_strs) {
      auto it_path = std::filesystem::path(std::move(it));
      if (std::filesystem::is_regular_file(it_path)) {
        paths.emplace_back(std::move(it_path));
      } else if (std::filesystem::is_directory(it_path)) {
        for (auto dir_entry : std::filesystem::recursive_directory_iterator(
                 std::move(it_path))) {
          if (std::filesystem::is_regular_file(dir_entry)) {
            paths.emplace_back(std::move(dir_entry));
          }
        }
      }
    }

    sequences = raven::LoadSequences(thread_pool, paths);

    if (sequences.empty()) {
      std::cerr << "[raven::] error: empty sequences set" << std::endl;
      return EXIT_FAILURE;
    }

    std::cerr << "[raven::] loaded " << sequences.size()
              << " sequences "  // NOLINT
              << std::fixed << timer.Stop() << "s" << std::endl;

    timer.Start();
  }

  raven::ConstructGraph(
      graph, sequences, thread_pool, checkpoints,
      raven::OverlapPhaseCfg{.kmer_len = kmer_len,
                             .window_len = window_len,
                             .freq = freq,
                             .identity = identity,
                             .kMaxNumOverlaps = kMaxNumOverlaps});

  raven::Assemble(thread_pool, graph, checkpoints);

  raven::Polish(
      thread_pool, graph, checkpoints, sequences,
      raven::PolishCfg{
          .align_cfg = raven::AlignCfg{.match = m, .mismatch = n, .gap = g},
          .cuda_cfg =
              raven::CudaCfg{.poa_batches = cuda_poa_batches,
                             .alignment_batches = cuda_alignment_batches,
                             .banded_alignment = cuda_banded_alignment},
          .num_rounds = num_polishing_rounds}

  );

  raven::PrintGfa(graph, gfa_path);

  for (const auto& it : raven::GetUnitigs(graph, num_polishing_rounds > 0)) {
    std::cout << ">" << it->name << std::endl;
    std::cout << it->InflateData() << std::endl;
  }

  timer.Stop();
  std::cerr << "[raven::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return EXIT_SUCCESS;
}
