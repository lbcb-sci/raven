// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <iostream>
#include <cstdlib>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/timer.hpp"
#include "raven/raven.h"
std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
  {"input-graph", required_argument, nullptr, 'i'},
  {"output-graph", required_argument, nullptr, 'o'}
};

void Help() {
  std::cout
      << "usage: graph_cleanup -i <input file in gfa format> -o <output file in gfa format>\n"
         "\n"
         "  options:\n"
         "    -i, --input-graph <string>\n"
         "      path to the file containing input graph in GFA format\n"
         "    -o, --output-graph <string>\n"
         "      path to the file where output graph will be stored in GFA format\n"
         "    -h, --help\n"
         "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::string input_gfa_path = "";
  std::string output_gfa_path = "";

  std::string optstr = "i:o:h";

  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {
    switch (arg) {
      case 'i':
        input_gfa_path = optarg;
        break;
      case 'o':
        output_gfa_path = optarg;
        break;
      case 'h':
        Help();
        return EXIT_SUCCESS;
      default:
        return EXIT_FAILURE;
    }
  }

  if (argc == 1) {
    Help();
    return EXIT_SUCCESS;
  }

  biosoup::Timer timer{};
  timer.Start();

  raven::Graph graph;

  graph = raven::LoadGfa(input_gfa_path);

  raven::RemoveInvalidEdgesFromGraph(graph);

  raven::RemoveInvalidConnectionsFromGraph(graph);

  raven::PrintGfa(graph, output_gfa_path);

  timer.Stop();
  std::cerr << "[graph_cleanup::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return EXIT_SUCCESS;
}
