#include <getopt.h>

#include <cstdint>
#include <string>
#include <iostream>
#include <exception>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"
#include "ram/sequence.hpp"

#include "graph.hpp"

static const std::string version = "v0.0.7";

static struct option options[] = {
    {"polishing-rounds", required_argument, nullptr, 'p'},
    {"match", required_argument, nullptr, 'm'},
    {"mismatch", required_argument, nullptr, 'n'},
    {"gap", required_argument, nullptr, 'g'},
#ifdef CUDA_ENABLED
    {"cuda-poa-batches", optional_argument, nullptr, 'c'},
    {"cuda-banded-alignment", no_argument, nullptr, 'b'},
    {"cuda-alignment-batches", required_argument, nullptr, 'a'},
#endif
    {"graphical-fragment-assembly", required_argument, nullptr, 'f'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path);

void help();

int main(int argc, char** argv) {

    std::uint32_t num_polishing_rounds = 2;
    std::int8_t m = 3;
    std::int8_t n = -5;
    std::int8_t g = -4;

    std::string gfa_path = "";

    std::uint32_t num_threads = 1;

    std::uint32_t cuda_poa_batches = 0;
    std::uint32_t cuda_alignment_batches = 0;
    bool cuda_banded_alignment = false;

    std::string optstring = "p:m:n:g:t:h";
#ifdef CUDA_ENABLED
    optstring += "c:b:a:";
#endif

    char argument;
    while ((argument = getopt_long(argc, argv, optstring.c_str(), options, nullptr)) != -1) {
        switch (argument) {
            case 'p': num_polishing_rounds = atoi(optarg); break;
            case 'm': m = atoi(optarg); break;
            case 'n': n = atoi(optarg); break;
            case 'g': g = atoi(optarg); break;
#ifdef CUDA_ENABLED
            case 'c':
                cuda_poa_batches = 1;
                // next text entry is not an option, assuming it's the arg for option 'c'
                if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
                    cuda_poa_batches = atoi(argv[optind++]);
                }
                // optional argument provided in the ususal way
                if (optarg != NULL) {
                    cuda_poa_batches = atoi(optarg);
                }
                break;
            case 'b':
                cuda_banded_alignment = true;
                break;
            case 'a':
                cuda_alignment_batches = atoi(optarg);
                break;
#endif
            case 'f': gfa_path = optarg; break;
            case 't': num_threads = atoi(optarg); break;
            case 'v': std::cout << version << std::endl; return 0;
            case 'h': help(); return 0;
            default: return 1;
        }
    }

    if (optind >= argc) {
        std::cerr << "[raven::] error: missing input file!" << std::endl;
        help();
        return 1;
    }

    auto sparser = createParser(argv[optind]);
    if (sparser == nullptr) {
        return 1;
    }

    logger::Logger logger;
    logger.log();

    std::vector<std::unique_ptr<ram::Sequence>> sequences;
    try {
        sparser->parse(sequences, -1);
    } catch (const std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    if (sequences.empty()) {
        std::cerr << "[raven::] error: empty sequences set" << std::endl;
        return 1;
    }

    logger.log("[raven::] loaded sequences");
    logger.log();

    std::shared_ptr<thread_pool::ThreadPool> thread_pool;
    try {
        thread_pool = thread_pool::createThreadPool(num_threads);
    } catch (const std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    auto graph = raven::createGraph(thread_pool);
    graph->construct(sequences);
    graph->assemble();
    graph->polish(sequences, m, n, g, cuda_poa_batches, cuda_banded_alignment,
        cuda_alignment_batches, num_polishing_rounds);
    graph->print_gfa(gfa_path);

    std::vector<std::unique_ptr<ram::Sequence>> unitigs;
    graph->get_unitigs(unitigs);
    for (const auto& it: unitigs) {
        std::cout << ">" << it->name << std::endl;
        std::cout << it->data << std::endl;
    }

    logger.log("[raven::]");

    return 0;
}

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path) {

    auto is_suffix = [] (const std::string& src, const std::string& suffix) {
        return src.size() < suffix.size() ? false :
            src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
        is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
        try {
            return bioparser::createParser<bioparser::FastaParser, ram::Sequence>(path);
        } catch (const std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            return nullptr;
        }
    }
    if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
        is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
        try {
            return bioparser::createParser<bioparser::FastqParser, ram::Sequence>(path);
        } catch (const std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            return nullptr;
        }
    }

    std::cerr << "[raven::] error: file " << path
              << " has unsupported format extension (valid extensions: .fasta, "
              << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!"
              << std::endl;
    return nullptr;
}

void help() {
    std::cout <<
        "usage: raven [options ...] <sequences>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences\n"
        "\n"
        "    options:\n"
        "        -p, --polishing-rounds <int>\n"
        "            default: 2\n"
        "            number of times racon is invoked\n"
        "        -m, --match <int>\n"
        "            default: 3\n"
        "            score for matching bases\n"
        "        -n, --mismatch <int>\n"
        "            default: -5\n"
        "            score for mismatching bases\n"
        "        -g, --gap <int>\n"
        "            default: -4\n"
        "            gap penalty (must be negative)\n"
        "        --graphical-fragment-assembly <string>\n"
        "            prints the assemblg graph in GFA format\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n"
        "    (only available when built with CUDA)\n"
        "        -c, --cuda-poa-batches\n"
        "            default: 1\n"
        "            number of batches for CUDA accelerated polishing\n"
        "        -b, --cuda-banded-alignment\n"
        "            use banding approximation for polishing on GPU\n"
        "            (only applicable when -c is used)\n"
        "        -a, --cuda-alignment-batches\n"
        "            number of batches for CUDA accelerated alignment\n";
}
