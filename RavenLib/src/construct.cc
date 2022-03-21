#include "raven/graph/construct.h"

#include <deque>

#include "biosoup/overlap.hpp"
#include "biosoup/timer.hpp"
#include "ram/minimizer_engine.hpp"
#include "raven/graph/overlap_utils.h"
#include "raven/graph/serialization/binary.h"
#include "edlib.h"

namespace raven {

void FindOverlapsAndCreatePiles(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    ram::MinimizerEngine& minimizer_engine,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    double freq, std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps) {
  constexpr std::size_t kMaxNumOverlaps = 32;

  piles.reserve(sequences.size());
  for (const auto& it : sequences) {
    piles.emplace_back(new Pile(it->id, it->inflated_len));
  }

  std::size_t bytes = 0;

  biosoup::Timer timer;

  for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
    bytes += sequences[i]->inflated_len;

    if (i != sequences.size() - 1ULL && bytes < (1ULL << 32ULL)) {
      continue;
    }

    bytes = 0;
    timer.Start();

    minimizer_engine.Minimize(sequences.begin() + j, sequences.begin() + i + 1,
                              true);
    minimizer_engine.Filter(freq);

    std::cerr << "[raven::Graph::Construct] minimized " << j << " - " << i + 1
              << " / " << sequences.size() << " " << std::fixed << timer.Stop()
              << "s" << std::endl;

    timer.Start();

    std::vector<std::uint32_t> num_overlaps(overlaps.size());
    for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
      num_overlaps[k] = overlaps[k].size();
    }

    std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;

    for (std::uint32_t k = 0; k < i + 1; ++k) {
      thread_futures.emplace_back(thread_pool->Submit(
          [&](std::uint32_t i) -> std::vector<biosoup::Overlap> {
            return minimizer_engine.Map(sequences[i], true, true, true);
          },
          k));

      bytes += sequences[k]->inflated_len;
      if (k != i && bytes < (1U << 30U)) {
        continue;
      }
      bytes = 0;

      for (auto& it : thread_futures) {
        for (const auto& jt : it.get()) {
          overlaps[jt.lhs_id].emplace_back(jt);
          overlaps[jt.rhs_id].emplace_back(OverlapReverse(jt));
        }
      }
      thread_futures.clear();

      std::vector<std::future<void>> void_futures;
      for (const auto& it : piles) {
        if (overlaps[it->id()].empty() ||
            overlaps[it->id()].size() == num_overlaps[it->id()]) {
          continue;
        }

        void_futures.emplace_back(thread_pool->Submit(
            [&](std::uint32_t i) -> void {
              piles[i]->AddLayers(overlaps[i].begin() + num_overlaps[i],
                                  overlaps[i].end());

              num_overlaps[i] = std::min(overlaps[i].size(), kMaxNumOverlaps);

              if (overlaps[i].size() < kMaxNumOverlaps) {
                return;
              }

              std::sort(overlaps[i].begin(), overlaps[i].end(),
                        [&](const biosoup::Overlap& lhs,
                            const biosoup::Overlap& rhs) -> bool {
                          return GetOverlapLength(lhs) > GetOverlapLength(rhs);
                        });

              std::vector<biosoup::Overlap> tmp;
              tmp.insert(tmp.end(), overlaps[i].begin(),
                         overlaps[i].begin() + kMaxNumOverlaps);  // NOLINT
              tmp.swap(overlaps[i]);
            },
            it->id()));
      }
      for (const auto& it : void_futures) {
        it.wait();
      }
    }

    std::cerr << "[raven::Graph::Construct] mapped sequences " << std::fixed
              << timer.Stop() << "s" << std::endl;

    j = i + 1;
  }
}

void TrimAndAnnotatePiles(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps) {
  biosoup::Timer timer;
  timer.Start();

  std::vector<std::future<void>> thread_futures;
  for (std::uint32_t i = 0; i < piles.size(); ++i) {
    thread_futures.emplace_back(thread_pool->Submit(
        [&](std::uint32_t i) -> void {
          piles[i]->FindValidRegion(4);
          if (piles[i]->is_invalid()) {
            std::vector<biosoup::Overlap>().swap(overlaps[i]);
          } else {
            piles[i]->FindMedian();
            piles[i]->FindChimericRegions();
          }
        },
        i));
  }

  for (const auto& it : thread_futures) {
    it.wait();
  }
  thread_futures.clear();

  std::cerr << "[raven::Graph::Construct] annotated piles " << std::fixed
            << timer.Stop() << "s" << std::endl;
}

void ResolveContainedReads(
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    double identity) {
  biosoup::Timer timer;

  timer.Start();

  std::vector<std::future<void>> futures;
  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    futures.emplace_back(thread_pool->Submit(
        [&] (std::uint32_t i) -> void {
          std::uint32_t k = 0;
          for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
            if (!OverlapUpdate(overlaps[i][j], piles)) {
              continue;
            }

            const auto& it = overlaps[i][j];

            auto lhs = sequences[it.lhs_id]->InflateData(
                it.lhs_begin,
                it.lhs_end - it.lhs_begin);

            auto rhs = sequences[it.rhs_id]->InflateData(
                it.rhs_begin,
                it.rhs_end - it.rhs_begin);
            if (!it.strand) {
              biosoup::NucleicAcid rhs_{"", rhs};
              rhs_.ReverseAndComplement();
              rhs = rhs_.InflateData();
            }

            auto result = edlibAlign(
                lhs.c_str(), lhs.size(),
                rhs.c_str(), rhs.size(),
                edlibDefaultAlignConfig());

            auto score = result.status == EDLIB_STATUS_OK ?
                1. - static_cast<double>(result.editDistance) / std::max(lhs.size(), rhs.size()) :  // NOLINT
                0.;

            edlibFreeAlignResult(result);

            if (score < identity) {
              continue;
            }
            overlaps[i][k++] = overlaps[i][j];
          }
          overlaps[i].resize(k);
        },
        i));
  }
  for (const auto& it : futures) {
    it.wait();
  }

  std::cerr << "[raven::Graph::Construct] filtered overlaps "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  timer.Start();

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    std::uint32_t k = 0;
    for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
      if (!OverlapUpdate(overlaps[i][j], piles)) {
        continue;
      }
      std::uint32_t type = GetOverlapType(overlaps[i][j], piles);
      if (type == 1 && !piles[overlaps[i][j].rhs_id]->is_maybe_chimeric()) {
        piles[i]->set_is_contained();
      } else if (type == 2 && !piles[i]->is_maybe_chimeric()) {
        piles[overlaps[i][j].rhs_id]->set_is_contained();
      } else {
        overlaps[i][k++] = overlaps[i][j];
      }
    }
    overlaps[i].resize(k);
  }
  for (std::uint32_t i = 0; i < piles.size(); ++i) {
    if (piles[i]->is_contained()) {
      piles[i]->set_is_invalid();

      std::vector<biosoup::Overlap>().swap(overlaps[i]);
    }
  }

  std::cerr << "[raven::Graph::Construct] removed contained sequences "
            << std::fixed << timer.Stop() << "s" << std::endl;
}

void ResolveChimericSequences(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  biosoup::Timer timer;

  timer.Start();

  std::vector<std::uint16_t> medians;
  for (const auto& it : piles) {
    if (it->median() != 0) {
      medians.emplace_back(it->median());
    }
  }
  std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                   medians.end());
  std::uint16_t median = medians[medians.size() / 2];

  std::vector<std::future<void>> thread_futures;
  for (const auto& it : piles) {
    if (it->is_invalid()) {
      continue;
    }
    thread_futures.emplace_back(thread_pool->Submit(
        [&](std::uint32_t i) -> void {
          piles[i]->ClearChimericRegions(median);
          if (piles[i]->is_invalid()) {
            std::vector<biosoup::Overlap>().swap(overlaps[i]);
          }
        },
        it->id()));
  }
  for (const auto& it : thread_futures) {
    it.wait();
  }

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    std::uint32_t k = 0;
    for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
      if (OverlapUpdate(overlaps[i][j], piles)) {
        overlaps[i][k++] = overlaps[i][j];
      }
    }
    overlaps[i].resize(k);
  }

  for (const auto& it : overlaps) {
    for (const auto& jt : it) {
      std::uint32_t type = GetOverlapType(jt, piles);
      if (type == 1) {
        piles[jt.lhs_id]->set_is_contained();
        piles[jt.lhs_id]->set_is_invalid();
      } else if (type == 2) {
        piles[jt.rhs_id]->set_is_contained();
        piles[jt.rhs_id]->set_is_invalid();
      }
    }
  }

  overlaps.clear();

  std::cerr << "[raven::Graph::Construct] removed chimeric sequences "
            << std::fixed << timer.Stop() << "s" << std::endl;
}

void FindOverlapsAndRepetetiveRegions(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    ram::MinimizerEngine& minimizerEngine, double freq, std::uint8_t kmer_len, double identity,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  biosoup::Timer timer;

  std::sort(
      sequences.begin(), sequences.end(),
      [&](const std::unique_ptr<biosoup::NucleicAcid>& lhs,
          const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
        return piles[lhs->id]->is_invalid() <
                   piles[rhs->id]->is_invalid() ||  // NOLINT
               (piles[lhs->id]->is_invalid() == piles[rhs->id]->is_invalid() &&
                lhs->id < rhs->id);  // NOLINT
      });

  std::uint32_t s = 0;
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    if (piles[sequences[i]->id]->is_invalid()) {
      s = i;
      break;
    }
  }

  // map valid reads to each other
  overlaps.resize(sequences.size() + 1);
  std::size_t bytes = 0;
  for (std::uint32_t i = 0, j = 0; i < s; ++i) {
    bytes += sequences[i]->inflated_len;
    if (i != s - 1 && bytes < (1U << 30)) {
      continue;
    }
    bytes = 0;

    timer.Start();

    minimizerEngine.Minimize(sequences.begin() + j, sequences.begin() + i + 1);

    std::cerr << "[raven::Graph::Construct] minimized " << j << " - " << i + 1
              << " / " << s << " " << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
    minimizerEngine.Filter(freq);
    for (std::uint32_t k = 0; k < i + 1; ++k) {
      thread_futures.emplace_back(thread_pool->Submit(
          [&](std::uint32_t i) -> std::vector<biosoup::Overlap> {
            std::vector<std::uint32_t> filtered;
            auto dst = minimizerEngine.Map(sequences[i],
                                           true,   // avoid equal
                                           true,   // avoid symmetric
                                           false,  // minhash
                                           &filtered);
            piles[sequences[i]->id]->AddKmers(filtered, kmer_len,
                                              sequences[i]);  // NOLINT
          
            std::uint32_t k = 0;
            for (std::uint32_t j = 0; j < dst.size(); ++j) {
              if (!OverlapUpdate(dst[j], piles)) {
                continue;
              }

              const auto& jt = dst[j];

              auto lhs = sequences[sequences_map[jt.lhs_id]]->InflateData(
                  jt.lhs_begin,
                  jt.lhs_end - jt.lhs_begin);

              auto rhs = sequences[sequences_map[jt.rhs_id]]->InflateData(
                  jt.rhs_begin,
                  jt.rhs_end - jt.rhs_begin);
              if (!jt.strand) {
                biosoup::NucleicAcid rhs_{"", rhs};
                rhs_.ReverseAndComplement();
                rhs = rhs_.InflateData();
              }

              auto result = edlibAlign(
                  lhs.c_str(), lhs.size(),
                  rhs.c_str(), rhs.size(),
                  edlibDefaultAlignConfig());

              auto score = result.status == EDLIB_STATUS_OK ?
                  1. - static_cast<double>(result.editDistance) / std::max(lhs.size(), rhs.size()) :  // NOLINT
                  0.;

              edlibFreeAlignResult(result);

              if (score < identity) {
                continue;
              }
              dst[k++] = dst[j];
            }
            dst.resize(k);
            
            return dst;
          },
          k));
    }
    for (auto& it : thread_futures) {
      for (auto& jt : it.get()) {
        if (!OverlapUpdate(jt, piles)) {
          continue;
        }
        std::uint32_t type = GetOverlapType(jt, piles);
        if (type == 0) {
          continue;
        } else if (type == 1) {
          piles[jt.lhs_id]->set_is_contained();
        } else if (type == 2) {
          piles[jt.rhs_id]->set_is_contained();
        } else {
          if (!overlaps.back().empty() &&
              overlaps.back().back().lhs_id == jt.lhs_id &&
              overlaps.back().back().rhs_id == jt.rhs_id) {
            if (GetOverlapLength(overlaps.back().back()) <
                GetOverlapLength(jt)) {
              overlaps.back().back() = jt;
            }
          } else {
            overlaps.back().emplace_back(jt);
          }
        }
      }
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] mapped valid sequences "
              << std::fixed << timer.Stop() << "s" << std::endl;

    j = i + 1;
  }

  timer.Start();

  std::vector<std::future<void>> thread_futures;
  for (const auto& pile : piles) {
    if (pile->is_contained()) {
      pile->set_is_invalid();
    }
  }

  {
    std::uint32_t k = 0;
    for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
      if (OverlapUpdate(overlaps.back()[i], piles)) {
        overlaps.back()[k++] = overlaps.back()[i];
      }
    }
    overlaps.back().resize(k);
  }

  std::cerr << "[raven::Graph::Construct] updated overlaps " << std::fixed
            << timer.Stop() << "s" << std::endl;

  std::sort(sequences.begin(), sequences.end(),
            [&](const std::unique_ptr<biosoup::NucleicAcid>& lhs,
                const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
              return lhs->id < rhs->id;
            });
}

void ResolveRepeatInducedOverlaps(
    const std::shared_ptr<thread_pool::ThreadPool>& thread_pool,
    const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  biosoup::Timer timer;

  timer.Start();

  while (true) {
    auto components = ConnectedComponents(overlaps, sequences, piles);
    for (const auto& it : components) {
      std::vector<std::uint16_t> medians;
      for (const auto& jt : it) {
        medians.emplace_back(piles[jt]->median());
      }
      std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                       medians.end());
      std::uint16_t median = medians[medians.size() / 2];

      std::vector<std::future<void>> futures;
      for (const auto& jt : it) {
        futures.emplace_back(thread_pool->Submit(
            [&](std::uint32_t i) -> void {
              piles[i]->FindRepetitiveRegions(median);
            },
            jt));
      }
      for (const auto& it : futures) {
        it.wait();
      }
    }

    for (const auto& it : overlaps.back()) {
      piles[it.lhs_id]->UpdateRepetitiveRegions(it);
      piles[it.rhs_id]->UpdateRepetitiveRegions(it);
    }

    bool is_changed = false;
    std::uint32_t j = 0;
    for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
      const auto& it = overlaps.back()[i];
      if (piles[it.lhs_id]->CheckRepetitiveRegions(it) ||
          piles[it.rhs_id]->CheckRepetitiveRegions(it)) {
        is_changed = true;
      } else {
        overlaps.back()[j++] = it;
      }
    }
    overlaps.back().resize(j);

    if (!is_changed) {
      break;
    }

    for (const auto& it : components) {
      for (const auto& jt : it) {
        piles[jt]->ClearRepetitiveRegions();
      }
    }
  }

  std::cerr << "[raven::Graph::Construct] removed false overlaps " << std::fixed
            << timer.Stop() << "s" << std::endl;

  timer.Start();
}

void ConstructAssemblyGraph(
    Graph& graph, const std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  biosoup::Timer timer;
  timer.Start();

  std::vector<std::int32_t> sequence_to_node(piles.size(), -1);

  std::vector<std::unique_ptr<Node>> nodes;
  std::vector<std::unique_ptr<Edge>> edges;

  auto emplace_node_through_factory = [&nodes,
                                       &graph](auto&&... args) -> Node* {
    return nodes
        .emplace_back(graph.node_factory.MakeUnique(
            std::forward<decltype(args)>(args)...))
        .get();
  };

  auto emplace_edge_through_factory = [&edges,
                                       &graph](auto&&... args) -> Edge* {
    return edges
        .emplace_back(graph.edge_factory.MakeUnique(
            std::forward<decltype(args)>(args)...))
        .get();
  };

  for (const auto& it : piles) {  // create nodes
    if (it->is_invalid()) {
      continue;
    }

    auto sequence = biosoup::NucleicAcid{
        sequences[it->id()]->name,
        sequences[it->id()]->InflateData(it->begin(),
                                         it->end() - it->begin())};  // NOLINT

    sequence_to_node[it->id()] = graph.node_factory.NextIndex();
    auto node = emplace_node_through_factory(sequence);

    sequence.ReverseAndComplement();
    auto rc_node = emplace_node_through_factory(sequence);

    node->pair = rc_node;
    rc_node->pair = node;
  }

  std::cerr << "[raven::Graph::Construct] stored " << nodes.size()
            << " nodes "  // NOLINT
            << std::fixed << timer.Stop() << "s" << std::endl;

  for (auto& it : overlaps.back()) {  // create edges
    if (!OverlapFinalize(it, piles)) {
      continue;
    }

    auto tail = nodes[sequence_to_node[it.lhs_id]].get();
    auto head = nodes[sequence_to_node[it.rhs_id] + 1 - it.strand].get();

    auto length = it.lhs_begin - it.rhs_begin;
    auto length_pair = (piles[it.rhs_id]->length() - it.rhs_end) -
                       (piles[it.lhs_id]->length() - it.lhs_end);

    if (it.score == 4) {
      std::swap(head, tail);
      length *= -1;
      length_pair *= -1;
    }

    auto edge = emplace_edge_through_factory(tail, head, length);
    auto rc_edge =
        emplace_edge_through_factory(head->pair, tail->pair, length_pair);

    edge->pair = rc_edge;
    rc_edge->pair = edge;
  }

  std::cerr << "[raven::Graph::Construct] stored " << edges.size() << " edges "
            << std::fixed << timer.Stop() << "s" << std::endl;

  graph.nodes = std::move(nodes);
  graph.edges = std::move(edges);
}

void ConstructGraph(
    Graph& graph, std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::shared_ptr<thread_pool::ThreadPool>& thread_pool, bool checkpoints,
    OverlapPhaseCfg const cfg) {
  if (sequences.empty() || graph.stage > -4) {
    return;
  }

  biosoup::Timer timer;
  timer.Start();

  ram::MinimizerEngine minimizerEngine{thread_pool, cfg.kmer_len,
                                       cfg.window_len};

  std::vector<std::vector<biosoup::Overlap>> overlaps;

  overlaps.resize(sequences.size());

  if (graph.stage == -5) {
    FindOverlapsAndCreatePiles(thread_pool, minimizerEngine, sequences,
                               cfg.freq, graph.piles, overlaps);
    TrimAndAnnotatePiles(thread_pool, graph.piles, overlaps);

    ResolveContainedReads(graph.piles, overlaps, sequences, thread_pool, cfg.identity);
    ResolveChimericSequences(thread_pool, graph.piles, overlaps, sequences);

    ++graph.stage;
    if (checkpoints) {
      biosoup::Timer localTimer;
      localTimer.Start();
      StoreGraphToFile(graph);
      std::cerr << "[raven::Graph::Construct] reached checkpoint " << std::fixed
                << localTimer.Stop() << "s" << std::endl;
    }
  }

  if (graph.stage == -4) {
    FindOverlapsAndRepetetiveRegions(thread_pool, minimizerEngine, cfg.freq,
                                     cfg.kmer_len, cfg.identity, graph.piles, overlaps,
                                     sequences);
    ResolveRepeatInducedOverlaps(thread_pool, graph.piles, overlaps, sequences);

    ConstructAssemblyGraph(graph, graph.piles, overlaps, sequences);

    ++graph.stage;
    if (checkpoints) {
      biosoup::Timer localTimer;
      localTimer.Start();

      StoreGraphToFile(graph);
      std::cerr << "[raven::Graph::Construct] reached checkpoint " << std::fixed
                << localTimer.Stop() << "s" << std::endl;
    }
  }

  std::cerr << "[raven::Graph::Construct] " << std::fixed << timer.Stop() << "s"
            << std::endl;
}

}  // namespace raven
