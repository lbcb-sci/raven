// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/timer.hpp"

#include "graph.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
  {"kmer-len", required_argument, nullptr, 'k'},
  {"window-len", required_argument, nullptr, 'w'},
  {"frequency", required_argument, nullptr, 'f'},
  {"polishing-rounds", required_argument, nullptr, 'p'},
  {"match", required_argument, nullptr, 'm'},
  {"mismatch", required_argument, nullptr, 'n'},
  {"gap", required_argument, nullptr, 'g'},
#ifdef CUDA_ENABLED
  {"cuda-poa-batches", optional_argument, nullptr, 'c'},
  {"cuda-banded-alignment", no_argument, nullptr, 'b'},
  {"cuda-alignment-batches", required_argument, nullptr, 'a'},
#endif
  {"graphical-fragment-assembly", required_argument, nullptr, 'F'},
  {"resume", no_argument, nullptr, 'r'},
  {"disable-checkpoints", no_argument, nullptr, 'd'},
  {"split", no_argument, nullptr, 's'},
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[raven::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: raven [options ...] <sequences> [<sequences> ...]\n"
      "\n"
      "  # default output is to stdout in FASTA format\n"
      "  <sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -k, --kmer-len <int>\n"
      "      default: 15\n"
      "      length of minimizers used to find overlaps\n"
      "    -w, --window-len <int>\n"
      "      default: 5\n"
      "      length of sliding window from which minimizers are sampled\n"
      "    -f, --frequency <double>\n"
      "      default: 0.001\n"
      "      threshold for ignoring most frequent minimizers\n"
      "    -p, --polishing-rounds <int>\n"
      "      default: 2\n"
      "      number of times racon is invoked\n"
      "    -m, --match <int>\n"
      "      default: 3\n"
      "      score for matching bases\n"
      "    -n, --mismatch <int>\n"
      "      default: -5\n"
      "      score for mismatching bases\n"
      "    -g, --gap <int>\n"
      "      default: -4\n"
      "      gap penalty (must be negative)\n"
#ifdef CUDA_ENABLED
      "    -c, --cuda-poa-batches <int>\n"
      "      default: 0\n"
      "      number of batches for CUDA accelerated polishing\n"
      "    -b, --cuda-banded-alignment\n"
      "      use banding approximation for polishing on GPU\n"
      "      (only applicable when -c is used)\n"
      "    -a, --cuda-alignment-batches <int>\n"
      "      default: 0\n"
      "      number of batches for CUDA accelerated alignment\n"
#endif
      "    --graphical-fragment-assembly <string>\n"
      "      prints the assembly graph in GFA format\n"
      "    --resume\n"
      "      resume previous run from last checkpoint\n"
      "    --disable-checkpoints\n"
      "      disable checkpoint file creation\n"
      "    --split\n"
      "      store read information and abort\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::uint8_t kmer_len = 15;
  std::uint8_t window_len = 5;
  double freq = 0.001;

  std::int32_t num_polishing_rounds = 2;
  std::int8_t m = 3;
  std::int8_t n = -5;
  std::int8_t g = -4;

  std::string gfa_path = "";
  bool resume = false;
  bool checkpoints = true;
  bool split = false;

  std::uint32_t num_threads = 1;

  std::uint32_t cuda_poa_batches = 0;
  std::uint32_t cuda_alignment_batches = 0;
  bool cuda_banded_alignment = false;

  std::string optstr = "k:w:f:p:m:n:g:t:h";
#ifdef CUDA_ENABLED
  optstr += "c:ba:";
#endif
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'k': kmer_len = std::atoi(optarg); break;
      case 'w': window_len = std::atoi(optarg); break;
      case 'f': freq = std::atof(optarg); break;
      case 'p': num_polishing_rounds = std::atoi(optarg); break;
      case 'm': m = std::atoi(optarg); break;
      case 'n': n = std::atoi(optarg); break;
      case 'g': g = std::atoi(optarg); break;
#ifdef CUDA_ENABLED
      case 'c':
        cuda_poa_batches = 1;
        // next text entry is not an option, assuming it's the arg for option c
        if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
          cuda_poa_batches = std::atoi(argv[optind++]);
        }
        // optional argument provided in the ususal way
        if (optarg != NULL) {
          cuda_poa_batches = std::atoi(optarg);
        }
        break;
      case 'b':
        cuda_banded_alignment = true;
        break;
      case 'a':
        cuda_alignment_batches = std::atoi(optarg);
        break;
#endif
      case 'F': gfa_path = optarg; break;
      case 'r': resume = true; break;
      case 'd': checkpoints = false; break;
      case 's': split = true; break;
      case 't': num_threads = std::atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc) {
    std::cerr << "[raven::] error: missing input file(s)!" << std::endl;
    return 1;
  }

  biosoup::Timer timer{};
  timer.Start();

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  raven::Graph graph{checkpoints, thread_pool};
  if (resume) {
    try {
      graph.Load();
    } catch (std::exception& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded previous run "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
  if (graph.stage() < -3 || num_polishing_rounds > std::max(0, graph.stage())) {
    for (int i = optind; i < argc; ++i) {
      auto sparser = CreateParser(argv[i]);
      if (sparser == nullptr) {
        return 1;
      }

      std::vector<std::unique_ptr<biosoup::NucleicAcid>> chunk;
      try {
        chunk = sparser->Parse(-1);
      } catch (const std::invalid_argument& exception) {
        std::cerr << exception.what() << " (" << argv[i] << ")"
                  << std::endl;
        return 1;
      }

      if (chunk.empty()) {
        std::cerr << "[raven::] warning: file " << argv[i] << " is empty"
                  << std::endl;
        continue;
      }

      if (sequences.empty()) {
        sequences.swap(chunk);
      } else {
        sequences.insert(
            sequences.end(),
            std::make_move_iterator(chunk.begin()),
            std::make_move_iterator(chunk.end()));
      }
    }

    if (sequences.empty()) {
      std::cerr << "[raven::] error: empty sequences set" << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded " << sequences.size() << " sequences "  // NOLINT
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();
  }

  graph.Construct(sequences, split, kmer_len, window_len, freq);
  graph.Assemble();
  graph.Polish(sequences, m, n, g, cuda_poa_batches, cuda_banded_alignment,
      cuda_alignment_batches, num_polishing_rounds);
  graph.PrintGfa(gfa_path);

  for (const auto& it : graph.GetUnitigs(num_polishing_rounds > 0)) {
    std::cout << ">" << it->name << std::endl;
    std::cout << it->InflateData() << std::endl;
  }

  timer.Stop();
  std::cerr << "[raven::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
