#include <getopt.h>

#include <cstdint>
#include <string>
#include <iostream>
#include <exception>

static const std::string version = "v0.0.0";

static struct option options[] = {
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

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

    return 0;
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
