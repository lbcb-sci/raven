#include "raven/graph/serialization/graph_repr.h"

namespace raven {

void PrintGfa(const Graph& graph, const std::string& path) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << "S\t" << it->sequence.name << "\t" << it->sequence.InflateData()
       << "\tLN:i:" << it->sequence.inflated_len << "\tRC:i:" << it->count
       << std::endl;
    if (it->is_circular) {
      os << "L\t" << it->sequence.name << "\t" << '+' << "\t"
         << it->sequence.name << "\t" << '+' << "\t0M" << std::endl;
    }
  }
  for (const auto& it : graph.edges) {
    if (it == nullptr || it->is_rc()) {
      continue;
    }
    os << "L\t" << it->tail->sequence.name << "\t"
       << (it->tail->is_rc() ? '-' : '+')  // NOLINT
       << "\t" << it->head->sequence.name << "\t"
       << (it->head->is_rc() ? '-' : '+')  // NOLINT
       << "\t" << it->tail->sequence.inflated_len - it->length << 'M'
       << std::endl;
  }
  os.close();
}

void PrintCsv(const Graph& graph, const std::string& path) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->sequence.inflated_len << " RC:i:" << it->count << ","
       << it->pair->id << " [" << it->pair->id / 2 << "]"
       << " LN:i:" << it->pair->sequence.inflated_len
       << " RC:i:" << it->pair->count << ",0,-" << std::endl;
  }
  for (const auto& it : graph.edges) {
    if (it == nullptr) {
      continue;
    }
    os << it->tail->id << " [" << it->tail->id / 2 << "]"
       << " LN:i:" << it->tail->sequence.inflated_len
       << " RC:i:" << it->tail->count << "," << it->head->id << " ["
       << it->head->id / 2 << "]"
       << " LN:i:" << it->head->sequence.inflated_len
       << " RC:i:" << it->head->count << ",1," << it->id << " " << it->length
       << " " << it->weight << std::endl;
  }
  for (const auto& it : graph.nodes) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->sequence.inflated_len << " RC:i:" << it->count << ","
       << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->sequence.inflated_len << " RC:i:" << it->count
       << ",1,-" << std::endl;
  }
  os.close();
}

void PrintJson(const raven::Graph& graph, const std::string& path) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  cereal::JSONOutputArchive archive(os);

  for (const auto& it : graph.piles) {
    if (it->is_invalid()) {
      continue;
    }

    archive(cereal::make_nvp(std::to_string(it->id()), *(it.get())));
  }
}

}  // namespace raven
