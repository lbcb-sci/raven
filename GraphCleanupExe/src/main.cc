// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <iostream>
#include <cstdlib>
#include <map>

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

  std::map<std::string, long> nodeLengths;
  std::map<std::uint32_t, std::string> additionalHifiasmEdgeInfo;

  graph = raven::LoadHifiasmGfa(input_gfa_path, nodeLengths, additionalHifiasmEdgeInfo);

  raven::PrintGfa(graph, output_gfa_path + ".raven.gfa");

  raven::RemoveInvalidEdgesFromGraph(graph);

  raven::PrintGfa(graph, output_gfa_path + ".removedInvalidEdges.raven.gfa");
  raven::PrintHifiasmGfa(graph, output_gfa_path + ".removedInvalidEdges.hifiasm.gfa", nodeLengths, additionalHifiasmEdgeInfo);

  raven::RemoveInvalidConnectionsFromGraph(graph);

  raven::PrintGfa(graph, output_gfa_path + ".removedInvalidConnections.raven.gfa");
  raven::PrintHifiasmGfa(graph, output_gfa_path + ".removedInvalidConnections.hifiasm.gfa", nodeLengths, additionalHifiasmEdgeInfo);

  timer.Stop();
  std::cerr << "[graph_cleanup::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return EXIT_SUCCESS;
}
