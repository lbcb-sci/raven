#include <getopt.h>

#include <cstdint>
#include <string>
#include <iostream>
#include <exception>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"
#include "ram/sequence.hpp"
#include "racon/racon.hpp"

#include "graph.hpp"

static const std::string version = "v0.0.0";

static struct option options[] = {
    {"polishing-rounds", required_argument, nullptr, 'p'},
#ifdef CUDA_ENABLED
    {"cudapoa-batches", optional_argument, nullptr, 'c'},
    {"cuda-banded-alignment", no_argument, nullptr, 'b'},
    {"cudaaligner-batches", required_argument, nullptr, 'a'},
#endif
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path);

void help();

int main(int argc, char** argv) {

    std::uint32_t num_polishing_rounds = 2;
    std::uint32_t num_threads = 1;

    std::uint32_t cudapoa_batches = 0;
    std::uint32_t cudaaligner_batches = 0;
    bool cuda_banded_alignment = false;

    std::string optstring = "p:t:vh";
#ifdef CUDA_ENABLED
    optstring += "c:b:a:";
#endif

    char argument;
    while ((argument = getopt_long(argc, argv, optstring.c_str(), options, nullptr)) != -1) {
        switch (argument) {
            case 'p': num_polishing_rounds = atoi(optarg); break;
#ifdef CUDA_ENABLED
            case 'c':
                cudapoa_batches = 1;
                // next text entry is not an option, assuming it's the arg for option 'c'
                if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
                    cudapoa_batches = atoi(argv[optind++]);
                }
                // optional argument provided in the ususal way
                if (optarg != NULL) {
                    cudapoa_batches = atoi(optarg);
                }
                break;
            case 'b':
                cuda_banded_alignment = true;
                break;
            case 'a':
                cudaaligner_batches = atoi(optarg);
                break;
#endif
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

    std::vector<std::unique_ptr<ram::Sequence>> contigs;
    graph->assemble(contigs);

    double q = 0;
    for (const auto& it: sequences) {
        if (it->quality.empty()) {
            continue;
        }
        double p = 0;
        for (const auto& jt: it->quality) {
            p += jt - 33;
        }
        q += p / it->quality.size();
    }
    q /= sequences.size();

    auto polisher = racon::createPolisher(q, 0.3, 500, true, 3, -5, -4,
        thread_pool, cudapoa_batches, cuda_banded_alignment,
        cudaaligner_batches);

    for (std::uint32_t i = 0; i < num_polishing_rounds; ++i) {
        std::vector<std::unique_ptr<ram::Sequence>> polished;
        polisher->initialize(contigs, sequences);
        polisher->polish(polished, true);
        contigs.swap(polished);
    }

    for (const auto& it: contigs) {
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
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n"
        "    (only available when built with CUDA)\n"
        "        -c, --cudapoa-batches\n"
        "            default: 1\n"
        "            number of batches for CUDA accelerated polishing\n"
        "        -b, --cuda-banded-alignment\n"
        "            use banding approximation for polishing on GPU\n"
        "            (only applicable when -c is used)\n"
        "        -a, --cudaaligner-batches\n"
        "            number of batches for CUDA accelerated alignment\n";
}
