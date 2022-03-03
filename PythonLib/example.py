import sys
import argparse

import ravenpy

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='ravenpy demo program')
  parser.add_argument(
    '--threads', type=int, nargs='?', default=1,
    help='number of threads available to the program')

  parser.add_argument(
    '--polish', type=int, nargs='?', default=2,
    help='number of polishing rounds done with racon')

  parser.add_argument(
    'paths', type=str, nargs='+', help='input read sequences')

  args = parser.parse_args()

  seq_paths = args.paths
  num_polishing_rounds = args.polish

  thread_pool = ravenpy.ThreadPool(args.threads)
  seqs_handle = ravenpy.SequencesHandle(seq_paths)

  graph = ravenpy.Graph()

  ravenpy.construct_graph(
    graph, seqs_handle, thread_pool, False, 
    ravenpy.OverlapPhaseCfg(15, 5, 0.001))

  ravenpy.assemble_graph(thread_pool, graph, False)

  ravenpy.polish_graph(
    thread_pool, graph, False, seqs_handle,
    ravenpy.PolishCfg(
        ravenpy.AlignCfg(3, -5, -4),
        ravenpy.CudaCfg(0, 0, False),
        num_polishing_rounds
      )
  )

  ravenpy.graph_print_unitgs(graph, num_polishing_rounds)
