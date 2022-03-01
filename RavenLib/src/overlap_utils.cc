#include "raven/graph/overlap_utils.h"

namespace raven {

biosoup::Overlap OverlapReverse(const biosoup::Overlap& o) {
  return {o.rhs_id,    o.rhs_begin, o.rhs_end, o.lhs_id,
          o.lhs_begin, o.lhs_end,   o.score,   o.strand};
}

std::uint32_t GetOverlapLength(const biosoup::Overlap& o) {
  return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
}

bool OverlapUpdate(biosoup::Overlap& o,
                   const std::vector<std::unique_ptr<Pile>>& piles) {
  if (piles[o.lhs_id]->is_invalid() || piles[o.rhs_id]->is_invalid()) {
    return false;
  }

  if (o.lhs_begin >= piles[o.lhs_id]->end() ||
      o.lhs_end <= piles[o.lhs_id]->begin() ||
      o.rhs_begin >= piles[o.rhs_id]->end() ||
      o.rhs_end <= piles[o.rhs_id]->begin()) {
    return false;
  }

  std::uint32_t lhs_begin =
      o.lhs_begin + (o.strand ? (o.rhs_begin < piles[o.rhs_id]->begin()
                                     ? piles[o.rhs_id]->begin() - o.rhs_begin
                                     : 0)
                              : (o.rhs_end > piles[o.rhs_id]->end()
                                     ? o.rhs_end - piles[o.rhs_id]->end()
                                     : 0));
  std::uint32_t lhs_end =
      o.lhs_end - (o.strand ? (o.rhs_end > piles[o.rhs_id]->end()
                                   ? o.rhs_end - piles[o.rhs_id]->end()
                                   : 0)
                            : (o.rhs_begin < piles[o.rhs_id]->begin()
                                   ? piles[o.rhs_id]->begin() - o.rhs_begin
                                   : 0));

  std::uint32_t rhs_begin =
      o.rhs_begin + (o.strand ? (o.lhs_begin < piles[o.lhs_id]->begin()
                                     ? piles[o.lhs_id]->begin() - o.lhs_begin
                                     : 0)
                              : (o.lhs_end > piles[o.lhs_id]->end()
                                     ? o.lhs_end - piles[o.lhs_id]->end()
                                     : 0));
  std::uint32_t rhs_end =
      o.rhs_end - (o.strand ? (o.lhs_end > piles[o.lhs_id]->end()
                                   ? o.lhs_end - piles[o.lhs_id]->end()
                                   : 0)
                            : (o.lhs_begin < piles[o.lhs_id]->begin()
                                   ? piles[o.lhs_id]->begin() - o.lhs_begin
                                   : 0));

  if (lhs_begin >= piles[o.lhs_id]->end() ||
      lhs_end <= piles[o.lhs_id]->begin() ||
      rhs_begin >= piles[o.rhs_id]->end() ||
      rhs_end <= piles[o.rhs_id]->begin()) {
    return false;
  }

  lhs_begin = std::max(lhs_begin, piles[o.lhs_id]->begin());
  lhs_end = std::min(lhs_end, piles[o.lhs_id]->end());
  rhs_begin = std::max(rhs_begin, piles[o.rhs_id]->begin());
  rhs_end = std::min(rhs_end, piles[o.rhs_id]->end());

  if (lhs_begin >= lhs_end || lhs_end - lhs_begin < 84 ||
      rhs_begin >= rhs_end || rhs_end - rhs_begin < 84) {
    return false;
  }

  o.lhs_begin = lhs_begin;
  o.lhs_end = lhs_end;
  o.rhs_begin = rhs_begin;
  o.rhs_end = rhs_end;

  return true;
}

std::uint32_t GetOverlapType(const biosoup::Overlap& o,
                             const std::vector<std::unique_ptr<Pile>>& piles) {
  std::uint32_t lhs_length = piles[o.lhs_id]->end() - piles[o.lhs_id]->begin();
  std::uint32_t lhs_begin = o.lhs_begin - piles[o.lhs_id]->begin();
  std::uint32_t lhs_end = o.lhs_end - piles[o.lhs_id]->begin();

  std::uint32_t rhs_length = piles[o.rhs_id]->end() - piles[o.rhs_id]->begin();
  std::uint32_t rhs_begin =
      o.strand ? o.rhs_begin - piles[o.rhs_id]->begin()
               : rhs_length - (o.rhs_end - piles[o.rhs_id]->begin());
  std::uint32_t rhs_end =
      o.strand ? o.rhs_end - piles[o.rhs_id]->begin()
               : rhs_length - (o.rhs_begin - piles[o.rhs_id]->begin());

  std::uint32_t overhang = std::min(lhs_begin, rhs_begin) +
                           std::min(lhs_length - lhs_end, rhs_length - rhs_end);

  if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
      rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
    return 0;  // internal
  }
  if (lhs_begin <= rhs_begin && lhs_length - lhs_end <= rhs_length - rhs_end) {
    return 1;  // lhs contained
  }
  if (rhs_begin <= lhs_begin && rhs_length - rhs_end <= lhs_length - lhs_end) {
    return 2;  // rhs contained
  }
  if (lhs_begin > rhs_begin) {
    return 3;  // lhs -> rhs
  }
  return 4;  // rhs -> lhs
}

 bool OverlapFinalize(biosoup::Overlap& o,
                            const std::vector<std::unique_ptr<Pile>>& piles) {
  o.score = GetOverlapType(o, piles);
  if (o.score < 3) {
    return false;
  }

  o.lhs_begin -= piles[o.lhs_id]->begin();
  o.lhs_end -= piles[o.lhs_id]->begin();

  o.rhs_begin -= piles[o.rhs_id]->begin();
  o.rhs_end -= piles[o.rhs_id]->begin();
  if (!o.strand) {
    auto rhs_begin = o.rhs_begin;
    o.rhs_begin = piles[o.rhs_id]->length() - o.rhs_end;
    o.rhs_end = piles[o.rhs_id]->length() - rhs_begin;
  }
  return true;
}

 std::vector<std::vector<std::uint32_t>> ConnectedComponents(
    const std::vector<std::vector<biosoup::Overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::vector<std::unique_ptr<Pile>>& piles) {
  std::vector<std::vector<std::uint32_t>> connections(sequences.size());

  for (const auto& it : overlaps) {
    for (const auto& jt : it) {
      if (GetOverlapType(jt, piles) > 2) {
        connections[jt.lhs_id].emplace_back(jt.rhs_id);
        connections[jt.rhs_id].emplace_back(jt.lhs_id);
      }
    }
  }

  std::vector<std::vector<std::uint32_t>> dst;
  std::vector<bool> isVisited(sequences.size(), false);

  for (std::uint32_t i = 0; i < connections.size(); ++i) {
    if (piles[i]->is_invalid() || isVisited[i]) {
      continue;
    }

    dst.resize(dst.size() + 1);
    std::deque<std::uint32_t> que = {i};

    while (!que.empty()) {
      std::uint32_t j = que.front();
      que.pop_front();

      if (isVisited[j]) {
        continue;
      }
      isVisited[j] = true;
      dst.back().emplace_back(j);

      for (const auto& it : connections[j]) {
        que.emplace_back(it);
      }
    }
  }

  return dst;
}

}  // namespace raven
