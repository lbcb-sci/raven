#include "raven/graph/polish.hpp"

#include "biosoup/timer.hpp"
#include "racon/polisher.hpp"
#include "raven/graph/common.h"
#include "raven/graph/serialization/binary.h"

namespace raven {

void Polish(std::shared_ptr<thread_pool::ThreadPool> thread_pool, Graph& graph,
            bool checkpoints,
            const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
            const PolishCfg& cfg) {
  if (sequences.empty() || cfg.num_rounds == 0) {
    return;
  }

  auto unitigs = GetUnitigs(graph);

  if (unitigs.empty()) {
    return;
  }

  graph.piles.clear();

  double avg_q = 0.;
  for (const auto& it : sequences) {
    if (it->block_quality.empty()) {
      continue;
    }
    double q =
        std::accumulate(it->block_quality.begin(), it->block_quality.end(), 0.);
    avg_q += q / it->block_quality.size();
  }
  if (avg_q == 0.) {  // when all values equal to '!'
    for (const auto& it : sequences) {
      it->block_quality.clear();
    }
  } else {
    avg_q /= sequences.size();
  }

  auto polisher = racon::Polisher::Create(
      thread_pool, avg_q, 0.3, 500, true, cfg.align_cfg.match,
      cfg.align_cfg.mismatch, cfg.align_cfg.gap,

      cfg.cuda_cfg.poa_batches, cfg.cuda_cfg.banded_alignment,
      cfg.cuda_cfg.alignment_batches);

  while (graph.stage < static_cast<std::int32_t>(cfg.num_rounds)) {
    auto polished = polisher->Polish(unitigs, sequences, false);
    unitigs.swap(polished);

    for (const auto& it : unitigs) {  // store unitigs
      const auto& node = graph.nodes[std::atoi(&it->name[3])];
      std::size_t tag;
      if ((tag = it->name.rfind(':')) != std::string::npos) {
        if (std::atof(&it->name[tag + 1]) > 0) {
          if (node->is_circular) {  // rotate
            auto s = it->InflateData();
            std::size_t b = 0.42 * s.size();
            s = s.substr(b) + s.substr(0, b);
            it->deflated_data = biosoup::NucleicAcid{"", s}.deflated_data;
          }

          node->is_polished = node->pair->is_polished = true;
          node->sequence.deflated_data = node->pair->sequence.deflated_data =
              it->deflated_data;  // NOLINT
          node->sequence.inflated_len = node->pair->sequence.inflated_len =
              it->inflated_len;  // NOLINT
        }
      }
    }

    ++graph.stage;
    if (checkpoints) {
      biosoup::Timer timer{};
      timer.Start();
      StoreGraphToFile(graph);

      std::cerr << "[raven::Graph::Polish] reached checkpoint " << std::fixed
                << timer.Stop() << "s" << std::endl;
    }
  }
}
}  // namespace raven
