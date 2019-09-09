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

static const std::string version = "v0.0.0";

static struct option options[] = {
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<ram::Sequence>> createParser(const std::string& path);

void help();

int main(int argc, char** argv) {

    std::uint32_t num_threads = 1;

    char argument;
    while ((argument = getopt_long(argc, argv, "k:w:f:t:h", options, nullptr)) != -1) {
        switch (argument) {
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
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n";
}
