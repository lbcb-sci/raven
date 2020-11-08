// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include <algorithm>
#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

#include "biosoup/timer.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"
#include "racon/polisher.hpp"

#include "edlib.h"  // NOLINT

namespace raven {

Graph::Node::Node(const biosoup::Sequence& sequence)
    : id(num_objects++),
      pid(sequence.id),
      name(sequence.name),
      data(sequence.data),
      count(1),
      is_circular(),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {}

Graph::Node::Node(Node* begin, Node* end)
    : id(num_objects++),
      pid(id % 2 ? end->pid : begin->pid),
      name(),
      data(),
      count(),
      is_circular(begin == end),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {

  auto it = begin;
  while (true) {
    data += it->outedges.front()->Label();
    count += it->count;
    if ((it = it->outedges.front()->head) == end) {
      break;
    }
  }
  if (begin != end) {
    data += end->data;
    count += end->count;
  }

  name = (is_unitig() ? "Utg" : "Ctg") + std::to_string(id & (~1UL));
}

Graph::Edge::Edge(Node* tail, Node* head, std::uint32_t length)
    : id(num_objects++),
      length(length),
      weight(0),
      tail(tail),
      head(head),
      pair() {
  tail->outedges.emplace_back(this);
  head->inedges.emplace_back(this);
}

std::atomic<std::uint32_t> Graph::Node::num_objects{0};
std::atomic<std::uint32_t> Graph::Edge::num_objects{0};

Graph::Graph(bool weaken, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
          thread_pool :
          std::make_shared<thread_pool::ThreadPool>(1)),
      minimizer_engine_(weaken ? 29 : 15, weaken ? 9 : 5, thread_pool_),
      stage_(-5),
      piles_(),
      nodes_(),
      edges_() {}

void Graph::Construct(std::vector<std::unique_ptr<biosoup::Sequence>>& sequences) {  // NOLINT
  if (sequences.empty() || stage_ > -4) {
    return;
  }

  std::vector<std::vector<biosoup::Overlap>> overlaps(sequences.size());

  // biosoup::Overlap helper functions
  auto overlap_reverse = [] (const biosoup::Overlap& o) -> biosoup::Overlap {
    return biosoup::Overlap(
        o.rhs_id, o.rhs_begin, o.rhs_end,
        o.lhs_id, o.lhs_begin, o.lhs_end,
        o.score,
        o.strand);
  };
  auto overlap_length = [] (const biosoup::Overlap& o) -> std::uint32_t {
    return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
  };
  auto overlap_update = [&] (biosoup::Overlap& o) -> bool {
    if (piles_[o.lhs_id]->is_invalid() ||
        piles_[o.rhs_id]->is_invalid()) {
      return false;
    }
    if (o.lhs_begin >= piles_[o.lhs_id]->end() ||
        o.lhs_end <= piles_[o.lhs_id]->begin() ||
        o.rhs_begin >= piles_[o.rhs_id]->end() ||
        o.rhs_end <= piles_[o.rhs_id]->begin()) {
      return false;
    }

    std::uint32_t lhs_begin = o.lhs_begin + (o.strand ?
        (o.rhs_begin < piles_[o.rhs_id]->begin() ?
            piles_[o.rhs_id]->begin() - o.rhs_begin : 0)
        :
        (o.rhs_end > piles_[o.rhs_id]->end() ?
            o.rhs_end - piles_[o.rhs_id]->end() : 0));
    std::uint32_t lhs_end = o.lhs_end - (o.strand ?
        (o.rhs_end > piles_[o.rhs_id]->end() ?
            o.rhs_end - piles_[o.rhs_id]->end() : 0)
        :
        (o.rhs_begin < piles_[o.rhs_id]->begin() ?
            piles_[o.rhs_id]->begin() - o.rhs_begin : 0));

    std::uint32_t rhs_begin = o.rhs_begin + (o.strand ?
        (o.lhs_begin < piles_[o.lhs_id]->begin() ?
            piles_[o.lhs_id]->begin() - o.lhs_begin : 0)
        :
        (o.lhs_end > piles_[o.lhs_id]->end() ?
            o.lhs_end - piles_[o.lhs_id]->end() : 0));
    std::uint32_t rhs_end = o.rhs_end - (o.strand ?
        (o.lhs_end > piles_[o.lhs_id]->end() ?
            o.lhs_end - piles_[o.lhs_id]->end() : 0)
        :
        (o.lhs_begin < piles_[o.lhs_id]->begin() ?
            piles_[o.lhs_id]->begin() - o.lhs_begin : 0));

    if (lhs_begin >= piles_[o.lhs_id]->end() ||
        lhs_end <= piles_[o.lhs_id]->begin() ||
        rhs_begin >= piles_[o.rhs_id]->end() ||
        rhs_end <= piles_[o.rhs_id]->begin()) {
      return false;
    }

    lhs_begin = std::max(lhs_begin, piles_[o.lhs_id]->begin());
    lhs_end = std::min(lhs_end, piles_[o.lhs_id]->end());
    rhs_begin = std::max(rhs_begin, piles_[o.rhs_id]->begin());
    rhs_end = std::min(rhs_end, piles_[o.rhs_id]->end());

    if (lhs_begin >= lhs_end || lhs_end - lhs_begin < 84 ||
        rhs_begin >= rhs_end || rhs_end - rhs_begin < 84) {
      return false;
    }

    o.lhs_begin = lhs_begin;
    o.lhs_end = lhs_end;
    o.rhs_begin = rhs_begin;
    o.rhs_end = rhs_end;

    return true;
  };
  auto overlap_type = [&] (const biosoup::Overlap& o) -> std::uint32_t {
    std::uint32_t lhs_length =
        piles_[o.lhs_id]->end() - piles_[o.lhs_id]->begin();
    std::uint32_t lhs_begin = o.lhs_begin - piles_[o.lhs_id]->begin();
    std::uint32_t lhs_end = o.lhs_end - piles_[o.lhs_id]->begin();

    std::uint32_t rhs_length =
        piles_[o.rhs_id]->end() - piles_[o.rhs_id]->begin();
    std::uint32_t rhs_begin = o.strand ?
        o.rhs_begin - piles_[o.rhs_id]->begin() :
        rhs_length - (o.rhs_end - piles_[o.rhs_id]->begin());
    std::uint32_t rhs_end = o.strand ?
        o.rhs_end - piles_[o.rhs_id]->begin():
        rhs_length - (o.rhs_begin - piles_[o.rhs_id]->begin());

    std::uint32_t overhang =
        std::min(lhs_begin, rhs_begin) +
        std::min(lhs_length - lhs_end, rhs_length - rhs_end);

    if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
        rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
      return 0;  // internal
    }
    if (lhs_begin <= rhs_begin &&
        lhs_length - lhs_end <= rhs_length - rhs_end) {
      return 1;  // lhs contained
    }
    if (rhs_begin <= lhs_begin &&
        rhs_length - rhs_end <= lhs_length - lhs_end) {
      return 2;  // rhs contained
    }
    if (lhs_begin > rhs_begin) {
      return 3;  // lhs -> rhs
    }
    return 4;  // rhs -> lhs
  };
  auto overlap_finalize = [&] (biosoup::Overlap& o) -> bool {
    o.score = overlap_type(o);
    if (o.score < 3) {
      return false;
    }

    o.lhs_begin -= piles_[o.lhs_id]->begin();
    o.lhs_end   -= piles_[o.lhs_id]->begin();

    o.rhs_begin -= piles_[o.rhs_id]->begin();
    o.rhs_end   -= piles_[o.rhs_id]->begin();
    if (!o.strand) {
       auto rhs_begin = o.rhs_begin;
       o.rhs_begin = piles_[o.rhs_id]->length() - o.rhs_end;
       o.rhs_end = piles_[o.rhs_id]->length() - rhs_begin;
    }
    return true;
  };

  auto connected_components = [&] () -> std::vector<std::vector<std::uint32_t>> {  // NOLINT
    std::vector<std::vector<std::uint32_t>> connections(sequences.size());
    for (const auto& it : overlaps) {
      for (const auto& jt : it) {
        if (overlap_type(jt) > 2) {
          connections[jt.lhs_id].emplace_back(jt.rhs_id);
          connections[jt.rhs_id].emplace_back(jt.lhs_id);
        }
      }
    }

    std::vector<std::vector<std::uint32_t>> dst;
    std::vector<char> is_visited(sequences.size(), false);
    for (std::uint32_t i = 0; i < connections.size(); ++i) {
      if (piles_[i]->is_invalid() || is_visited[i]) {
        continue;
      }

      dst.resize(dst.size() + 1);
      std::deque<std::uint32_t> que = { i };
      while (!que.empty()) {
        std::uint32_t j = que.front();
        que.pop_front();

        if (is_visited[j]) {
          continue;
        }
        is_visited[j] = true;
        dst.back().emplace_back(j);

        for (const auto& it : connections[j]) {
          que.emplace_back(it);
        }
      }
    }

    return dst;
  };
  // biosoup::Overlap helper functions

  if (stage_ == -5) {  // checkpoint test
    Store();
  }

  biosoup::Timer timer{};

  if (stage_ == -5) {  // find overlaps and create piles
    for (const auto& it : sequences) {
      piles_.emplace_back(new Pile(it->id, it->data.size()));
    }
    std::size_t bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
      bytes += sequences[i]->data.size();
      if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine_.Minimize(
          sequences.begin() + j,
          sequences.begin() + i + 1,
          true);
      minimizer_engine_.Filter(0.001);

      std::cerr << "[raven::Graph::Construct] minimized "
                << j << " - " << i + 1 << " / " << sequences.size() << " "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      std::vector<std::uint32_t> num_overlaps(overlaps.size());
      for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
        num_overlaps[k] = overlaps[k].size();
      }

      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;

      for (std::uint32_t k = 0; k < i + 1; ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&] (std::uint32_t i) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[i], true, true, true);
            },
            k));

        bytes += sequences[k]->data.size();
        if (k != i && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        for (auto& it : thread_futures) {
          for (const auto& jt : it.get()) {
            overlaps[jt.lhs_id].emplace_back(jt);
            overlaps[jt.rhs_id].emplace_back(overlap_reverse(jt));
          }
        }
        thread_futures.clear();

        std::vector<std::future<void>> void_futures;
        for (const auto& it : piles_) {
          if (overlaps[it->id()].empty() ||
              overlaps[it->id()].size() == num_overlaps[it->id()]) {
            continue;
          }

          void_futures.emplace_back(thread_pool_->Submit(
              [&] (std::uint32_t i) -> void {
                piles_[i]->AddLayers(
                    overlaps[i].begin() + num_overlaps[i],
                    overlaps[i].end());

                num_overlaps[i] = std::min(
                    overlaps[i].size(),
                    static_cast<std::size_t>(16));

                if (overlaps[i].size() < 16) {
                  return;
                }

                std::sort(overlaps[i].begin(), overlaps[i].end(),
                    [&] (const biosoup::Overlap& lhs,
                         const biosoup::Overlap& rhs) -> bool {
                      return overlap_length(lhs) > overlap_length(rhs);
                    });

                std::vector<biosoup::Overlap> tmp;
                tmp.insert(tmp.end(), overlaps[i].begin(), overlaps[i].begin() + 16);  // NOLINT
                tmp.swap(overlaps[i]);
              },
              it->id()));
        }
        for (const auto& it : void_futures) {
          it.wait();
        }
      }

      std::cerr << "[raven::Graph::Construct] mapped sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      j = i + 1;
    }
  }

  if (stage_ == -5) {  // trim and annotate piles
    timer.Start();

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      thread_futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> void {
            piles_[i]->FindValidRegion(4);
            if (piles_[i]->is_invalid()) {
              std::vector<biosoup::Overlap>().swap(overlaps[i]);
            } else {
              piles_[i]->FindMedian();
              piles_[i]->FindChimericRegions();
            }
          },
          i));
    }
    for (const auto& it : thread_futures) {
      it.wait();
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] annotated piles "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -5) {  // resolve contained reads
    timer.Start();

    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
      std::uint32_t k = 0;
      for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
        if (!overlap_update(overlaps[i][j])) {
          continue;
        }
        std::uint32_t type = overlap_type(overlaps[i][j]);
        if (type == 1 &&
            !piles_[overlaps[i][j].rhs_id]->is_maybe_chimeric()) {
          piles_[i]->set_is_contained();
        } else if (type == 2 &&
            !piles_[i]->is_maybe_chimeric()) {
          piles_[overlaps[i][j].rhs_id]->set_is_contained();
        } else {
          overlaps[i][k++] = overlaps[i][j];
        }
      }
      overlaps[i].resize(k);
    }
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      if (piles_[i]->is_contained()) {
        piles_[i]->set_is_invalid();
        std::vector<biosoup::Overlap>().swap(overlaps[i]);
      }
    }

    std::cerr << "[raven::Graph::Construct] removed contained sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -5) {  // resolve chimeric sequences
    timer.Start();

    while (true) {
      auto components = connected_components();
      for (const auto& it : components) {
        std::vector<std::uint32_t> medians;
        for (const auto& jt : it) {
          medians.emplace_back(piles_[jt]->median());
        }
        std::nth_element(
            medians.begin(),
            medians.begin() + medians.size() / 2,
            medians.end());
        std::uint32_t median = medians[medians.size() / 2];

        std::vector<std::future<void>> thread_futures;
        for (const auto& jt : it) {
          thread_futures.emplace_back(thread_pool_->Submit(
              [&] (std::uint32_t i) -> void {
                piles_[i]->ClearChimericRegions(median);
                if (piles_[i]->is_invalid()) {
                  std::vector<biosoup::Overlap>().swap(overlaps[i]);
                }
              },
              jt));
          }
        for (const auto& it : thread_futures) {
          it.wait();
        }
        thread_futures.clear();
      }

      bool is_changed = false;
      for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        std::uint32_t k = 0;
        for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
          if (overlap_update(overlaps[i][j])) {
            overlaps[i][k++] = overlaps[i][j];
          } else {
            is_changed = true;
          }
        }
        overlaps[i].resize(k);
      }

      if (!is_changed) {
        for (const auto& it : overlaps) {
          for (const auto& jt : it) {
            std::uint32_t type = overlap_type(jt);
            if (type == 1) {
              piles_[jt.lhs_id]->set_is_contained();
              piles_[jt.lhs_id]->set_is_invalid();
            } else if (type == 2) {
              piles_[jt.rhs_id]->set_is_contained();
              piles_[jt.rhs_id]->set_is_invalid();
            }
          }
        }
        overlaps.clear();
        break;
      }
    }

    std::cerr << "[raven::Graph::Construct] removed chimeric sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -5) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Construct] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
  }

  if (stage_ == -4) {  // clear piles for sensitive overlaps
    timer.Start();

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      if (piles_[i]->is_invalid()) {
        continue;
      }
      thread_futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> void {
            piles_[i]->ClearValidRegion();
          },
          i));
    }
    for (const auto& it : thread_futures) {
      it.wait();
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] cleared piles "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -4) {  // find overlaps and update piles with repetitive regions
    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<biosoup::Sequence>& lhs,
             const std::unique_ptr<biosoup::Sequence>& rhs) -> bool {
          return piles_[lhs->id]->is_invalid() <  piles_[rhs->id]->is_invalid() ||  // NOLINT
                (piles_[lhs->id]->is_invalid() == piles_[rhs->id]->is_invalid() && lhs->id < rhs->id);  // NOLINT
        });

    std::uint32_t s = 0;
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
      if (piles_[sequences[i]->id]->is_invalid()) {
        s = i;
        break;
      }
    }

    // map invalid reads to valid reads
    overlaps.resize(sequences.size() + 1);
    std::size_t bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < s; ++i) {
      bytes += sequences[i]->data.size();
      if (i != s - 1 && bytes < (1ULL << 32)) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine_.Minimize(
          sequences.begin() + j,
          sequences.begin() + i + 1,
          true);

      std::cerr << "[raven::Graph::Construct] minimized "
                << j << " - " << i + 1 << " / " << s << " "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      minimizer_engine_.Filter(0.00001);
      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
      for (std::uint32_t k = s; k < sequences.size(); ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&] (std::uint32_t i) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[i], true, false, true);
            },
            k));

        bytes += sequences[k]->data.size();
        if (k != sequences.size() - 1 && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        for (auto& it : thread_futures) {
          for (const auto& jt : it.get()) {
            overlaps[jt.rhs_id].emplace_back(jt);
          }
        }
        thread_futures.clear();

        std::vector<std::future<void>> void_futures;
        for (std::uint32_t k = j; k < i + 1; ++k) {
          if (overlaps[sequences[k]->id].empty()) {
            continue;
          }
          void_futures.emplace_back(thread_pool_->Submit(
              [&] (std::uint32_t i) -> void {
                piles_[i]->AddLayers(overlaps[i].begin(), overlaps[i].end());
                std::vector<biosoup::Overlap>().swap(overlaps[i]);
              },
              sequences[k]->id));
        }
        for (const auto& it : void_futures) {
          it.wait();
        }
      }

      std::cerr << "[raven::Graph::Construct] mapped invalid sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      j = i + 1;
    }

    // map valid reads to each other
    bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < s; ++i) {
      bytes += sequences[i]->data.size();
      if (i != s - 1 && bytes < (1U << 30)) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine_.Minimize(
          sequences.begin() + j,
          sequences.begin() + i + 1);

      std::cerr << "[raven::Graph::Construct] minimized "
                << j << " - " << i + 1 << " / " << s << " "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
      minimizer_engine_.Filter(0.001);
      for (std::uint32_t k = 0; k < i + 1; ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&] (std::uint32_t i) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[i], true, true);
            },
            k));
      }
      for (auto& it : thread_futures) {
        for (auto& jt : it.get()) {
          if (!overlap_update(jt)) {
            continue;
          }
          std::uint32_t type = overlap_type(jt);
          if (type == 0) {
            continue;
          } else if (type == 1) {
            piles_[jt.lhs_id]->set_is_contained();
          } else if (type == 2) {
            piles_[jt.rhs_id]->set_is_contained();
          } else {
            if (overlaps.back().size() &&
                overlaps.back().back().lhs_id == jt.lhs_id &&
                overlaps.back().back().rhs_id == jt.rhs_id) {
              if (overlap_length(overlaps.back().back()) < overlap_length(jt)) {
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
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      j = i + 1;
    }

    timer.Start();

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      if (piles_[i]->is_contained()) {
        piles_[i]->set_is_invalid();
        continue;
      }
      if (piles_[i]->is_invalid()) {
        continue;
      }
      thread_futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> void {
            piles_[i]->ClearInvalidRegion();
            piles_[i]->FindMedian();
          },
          i));
    }
    for (const auto& it : thread_futures) {
      it.wait();
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] updated piles "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    {
      std::uint32_t k = 0;
      for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
        if (overlap_update(overlaps.back()[i])) {
          overlaps.back()[k++] = overlaps.back()[i];
        }
      }
      overlaps.back().resize(k);
    }

    std::cerr << "[raven::Graph::Construct] updated overlaps "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<biosoup::Sequence>& lhs,
             const std::unique_ptr<biosoup::Sequence>& rhs) -> bool {
          return lhs->id < rhs->id;
        });
  }

  if (stage_ == -4) {  // resolve repeat induced overlaps
    timer.Start();

    while (true) {
      auto components = connected_components();
      for (const auto& it : components) {
        std::vector<std::uint32_t> medians;
        for (const auto& jt : it) {
          medians.emplace_back(piles_[jt]->median());
        }
        std::nth_element(
            medians.begin(),
            medians.begin() + medians.size() / 2,
            medians.end());
        std::uint32_t median = medians[medians.size() / 2];

        std::vector<std::future<void>> futures;
        for (const auto& jt : it) {
          futures.emplace_back(thread_pool_->Submit(
              [&] (std::uint32_t i) -> void {
                piles_[i]->FindRepetitiveRegions(median);
              },
              jt));
        }
        for (const auto& it : futures) {
          it.wait();
        }
      }

      for (const auto& it : overlaps.back()) {
        piles_[it.lhs_id]->UpdateRepetitiveRegions(it);
        piles_[it.rhs_id]->UpdateRepetitiveRegions(it);
      }

      bool is_changed = false;
      std::uint32_t j = 0;
      for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
        const auto& it = overlaps.back()[i];
        if (piles_[it.lhs_id]->CheckRepetitiveRegions(it) ||
            piles_[it.rhs_id]->CheckRepetitiveRegions(it)) {
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
          piles_[jt]->ClearRepetitiveRegions();
        }
      }
    }

    std::cerr << "[raven::Graph::Construct] removed false overlaps "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();
  }

  if (stage_ == -4) {  // construct assembly graph
    std::vector<std::int32_t> sequence_to_node(piles_.size(), -1);
    for (const auto& it : piles_) {  // create nodes
      if (it->is_invalid()) {
        continue;
      }

      auto sequence = biosoup::Sequence{
          sequences[it->id()]->name,
          sequences[it->id()]->data.substr(it->begin(), it->end() - it->begin())};  // NOLINT
      sequence.id = it->id();

      sequence_to_node[it->id()] = Node::num_objects;

      auto node = std::make_shared<Node>(sequence);
      sequence.ReverseAndComplement();
      nodes_.emplace_back(node);
      nodes_.emplace_back(std::make_shared<Node>(sequence));
      node->pair = nodes_.back().get();
      node->pair->pair = node.get();
    }

    std::cerr << "[raven::Graph::Construct] stored " << nodes_.size() << " nodes "  // NOLINT
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    for (auto& it : overlaps.back()) {  // create edges
      if (!overlap_finalize(it)) {
        continue;
      }

      auto tail = nodes_[sequence_to_node[it.lhs_id]].get();
      auto head = nodes_[sequence_to_node[it.rhs_id] + 1 - it.strand].get();

      auto length = it.lhs_begin - it.rhs_begin;
      auto length_pair =
          (piles_[it.rhs_id]->length() - it.rhs_end) -
          (piles_[it.lhs_id]->length() - it.lhs_end);

      if (it.score == 4) {
        std::swap(head, tail);
        length *= -1;
        length_pair *= -1;
      }

      auto edge = std::make_shared<Edge>(tail, head, length);
      edges_.emplace_back(edge);
      edges_.emplace_back(std::make_shared<Edge>(head->pair, tail->pair, length_pair));  // NOLINT
      edge->pair = edges_.back().get();
      edge->pair->pair = edge.get();
    }

    std::cerr << "[raven::Graph::Construct] stored " << edges_.size() << " edges "  // NOLINT
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -4) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Construct] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
  }

  std::cerr << "[raven::Graph::Construct] "
            << std::fixed <<  timer.elapsed_time() << "s"
            << std::endl;
}  // NOLINT

void Graph::Assemble() {
  if (stage_ < -3 || stage_ > -1) {
    return;
  }

  biosoup::Timer timer{};

  if (stage_ == -3) {  // remove transitive edges
    timer.Start();

    RemoveTransitiveEdges();

    std::cerr << "[raven::Graph::Assemble] removed transitive edges "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -3) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Assemble] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -2) {  // remove tips and bubbles
    timer.Start();

    while (true) {
      std::uint32_t num_changes = RemoveTips();
      num_changes += RemoveBubbles();
      if (num_changes == 0) {
        break;
      }
    }

    std::cerr << "[raven::Graph::Assemble] removed tips and bubbles "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -2) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Assemble] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -1) {  // remove long edges
    timer.Start();

    CreateUnitigs(42);  // speed up force directed layout
    RemoveLongEdges(16);

    std::cerr << "[raven::Graph::Assemble] removed long edges "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    CreateUnitigs(7);
    std::unordered_set<std::uint32_t> valid_piles;
    for (const auto& it : nodes_) {
      if (it) {
        valid_piles.emplace(it->pid);
      }
    }
    for (const auto& it : piles_) {
      if (valid_piles.count(it->id()) == 0) {
        it->set_is_invalid();
      }
    }
    PrintJSON("leftovers.json");
    PrintGFA("graph.gfa");

    std::ofstream os("neighbors.json");
    cereal::JSONOutputArchive archive(os);
    for (const auto& it : nodes_) {
      if (it && it->id % 2) {
        std::vector<std::uint32_t> neighbors;
        for (const auto& jt : it->outedges) {
          neighbors.emplace_back(jt->head->pid);
        }
        for (const auto& jt : it->inedges) {
          neighbors.emplace_back(jt->tail->pid);
        }
        archive(cereal::make_nvp(std::to_string(it->pid), neighbors));
      }
    }
  }

  if (stage_ == -1) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Assemble] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  timer.Start();

  while (true) {  // TODO(rvaser): check if necessary
     std::uint32_t num_changes = RemoveTips();
     num_changes += RemoveBubbles();
     if (num_changes == 0) {
       break;
     }
  }

  timer.Stop();
  std::cerr << "[raven::Graph::Assemble] "
            << std::fixed << timer.elapsed_time() << "s"
            << std::endl;
}

void Graph::RemoveMarked(const std::string& path) {
  if (path.empty()) {
    return;
  }

  biosoup::Timer timer{};
  timer.Start();

  std::ifstream is(path);
  if (is.is_open()) {
    std::unordered_set<std::uint32_t> indices;

    std::string tail, head;
    while (is >> tail >> head) {
      for (const auto& it : nodes_) {
        if (!it) {
          continue;
        }
        if (it->name == tail) {
          if (tail == head) {
            for (const auto& jt : it->inedges) {
              indices.emplace(jt->id);
              indices.emplace(jt->pair->id);
            }
            for (const auto& jt : it->outedges) {
              indices.emplace(jt->id);
              indices.emplace(jt->pair->id);
            }
          } else {
            for (const auto& jt : it->outedges) {
              if (jt->head->name == head) {
                indices.emplace(jt->id);
                indices.emplace(jt->pair->id);
                break;
              }
            }
            for (const auto& jt : it->inedges) {
              if (jt->tail->name == head) {
                indices.emplace(jt->id);
                indices.emplace(jt->pair->id);
                break;
              }
            }
          }
          break;
        }
      }
    }
    is.close();

    RemoveEdges(indices, true);
    PrintGFA("updated_graph.gfa");

    Store();
  }

  timer.Stop();
  std::cerr << "[raven::Graph::RemoveMarked] "
            << std::fixed << timer.elapsed_time() << "s"
            << std::endl;
}

std::uint32_t Graph::RemoveTransitiveEdges() {
  auto is_comparable = [] (double a, double b) -> bool {
    double eps = 0.12;
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
           (b >= a * (1 - eps) && b <= a * (1 + eps));
  };

  std::vector<Edge*> candidate(nodes_.size(), nullptr);
  std::unordered_set<std::uint32_t> marked_edges;
  for (const auto& it : nodes_) {
    if (it == nullptr) {
      continue;
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = jt;
    }
    for (auto jt : it->outedges) {
      for (auto kt : jt->head->outedges) {
        if (candidate[kt->head->id] &&
            is_comparable(jt->length + kt->length, candidate[kt->head->id]->length)) {  // NOLINT
          marked_edges.emplace(candidate[kt->head->id]->id);
          marked_edges.emplace(candidate[kt->head->id]->pair->id);
        }
      }
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = nullptr;
    }
  }

  for (auto i : marked_edges) {  // store for force directed layout
    if (i & 1) {
      auto lhs = edges_[i]->tail->id & ~1UL;
      auto rhs = edges_[i]->head->id & ~1UL;
      nodes_[lhs]->transitive.emplace(rhs);
      nodes_[rhs]->transitive.emplace(lhs);
    }
  }

  RemoveEdges(marked_edges);
  return marked_edges.size() / 2;
}

std::uint32_t Graph::RemoveTips() {
  std::uint32_t num_tips = 0;
  std::vector<char> is_visited(nodes_.size(), 0);

  for (const auto& it : nodes_) {
    if (it == nullptr || is_visited[it->id] || !it->is_tip()) {
      continue;
    }
    bool is_circular = false;
    std::uint32_t num_sequences = 0;

    auto end = it.get();
    while (!end->is_junction()) {
      num_sequences += end->count;
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 ||
          end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    if (is_circular || end->outdegree() == 0 || num_sequences > 5) {
      continue;
    }

    std::unordered_set<std::uint32_t> marked_edges;
    for (auto jt : end->outedges) {
      if (jt->head->indegree() > 1) {
        marked_edges.emplace(jt->id);
        marked_edges.emplace(jt->pair->id);
      }
    }
    if (marked_edges.size() / 2 == end->outedges.size()) {  // delete whole
      auto begin = it.get();
      while (begin != end) {
        marked_edges.emplace(begin->outedges.front()->id);
        marked_edges.emplace(begin->outedges.front()->pair->id);
        begin = begin->outedges.front()->head;
      }
      ++num_tips;
    }
    RemoveEdges(marked_edges, true);
  }

  return num_tips;
}

std::uint32_t Graph::RemoveBubbles() {
  std::vector<std::uint32_t> distance(nodes_.size(), 0);
  std::vector<Node*> predecessor(nodes_.size(), nullptr);

  // path helper functions
  auto path_extract = [&] (Node* begin, Node* end) -> std::vector<Node*> {
    std::vector<Node*> dst;
    while (end != begin) {
      dst.emplace_back(end);
      end = predecessor[end->id];
    }
    dst.emplace_back(begin);
    std::reverse(dst.begin(), dst.end());
    return dst;
  };
  auto path_type = [] (const std::vector<Node*>& path) -> bool {
    if (path.empty()) {
      return false;
    }
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
      if (path[i]->is_junction()) {
        return false;  // complex
      }
    }
    return true;  // without branches
  };
  auto bubble_type = [&] (
      const std::vector<Node*>& lhs,
      const std::vector<Node*>& rhs) -> bool {
    if (lhs.empty() || rhs.empty()) {
      return false;
    }
    std::unordered_set<Node*> intersection;
    for (auto it : lhs) {
      intersection.emplace(it);
    }
    for (auto it : rhs) {
      intersection.emplace(it);
    }
    if (lhs.size() + rhs.size() - 2 != intersection.size()) {
      return false;
    }
    for (auto it : lhs) {
      if (intersection.count(it->pair) != 0) {
        return false;
      }
    }
    if (path_type(lhs) && path_type(rhs)) {  // both without branches
      return true;
    }

    auto path_sequence = [] (const std::vector<Node*>& path) ->
        std::unique_ptr<biosoup::Sequence> {
      auto sequence =
          std::unique_ptr<biosoup::Sequence>(new biosoup::Sequence());
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        for (auto it : path[i]->outedges) {
          if (it->head == path[i + 1]) {
            sequence->data += it->Label();
            break;
          }
        }
      }
      sequence->data += path.back()->data;
      return sequence;
    };

    auto ls = path_sequence(lhs);
    auto rs = path_sequence(rhs);

    if (std::min(ls->data.size(), rs->data.size()) <
        std::max(ls->data.size(), rs->data.size()) * 0.8) {
      return false;
    }

    auto overlaps = minimizer_engine_.Map(ls, rs);

    std::uint32_t matches = 0;
    for (const auto& it : overlaps) {
      matches = std::max(matches, it.score);
    }
    return matches > 0.5 * std::min(ls->data.size(), rs->data.size());
  };
  // path helper functions

  std::uint32_t num_bubbles = 0;
  for (const auto& it : nodes_) {
    if (it == nullptr || it->outdegree() < 2) {
      continue;
    }

    // BFS
    Node* begin = it.get();
    Node* end = nullptr;
    Node* other_end = nullptr;
    std::deque<Node*> que{begin};
    std::vector<Node*> visited{1, begin};
    while (!que.empty() && !end) {
      auto jt = que.front();
      que.pop_front();

      for (auto kt : jt->outedges) {
        if (kt->head == begin) {  // cycle
          continue;
        }
        if (distance[jt->id] + kt->length > 500000) {  // out of reach
          continue;
        }
        distance[kt->head->id] = distance[jt->id] + kt->length;
        visited.emplace_back(kt->head);
        que.emplace_back(kt->head);

        if (predecessor[kt->head->id]) {  // found bubble
          end = kt->head;
          other_end = jt;
          break;
        }

        predecessor[kt->head->id] = jt;
      }
    }
    std::unordered_set<std::uint32_t> marked_edges;
    if (end) {
      auto lhs = path_extract(begin, end);
      auto rhs = path_extract(begin, other_end);
      rhs.emplace_back(end);

      if (bubble_type(lhs, rhs)) {
        std::uint32_t lhs_count = 0;
        for (auto jt : lhs) {
          lhs_count += jt->count;
        }
        std::uint32_t rhs_count = 0;
        for (auto jt : rhs) {
          rhs_count += jt->count;
        }
        marked_edges = FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
        if (marked_edges.empty()) {
          marked_edges = FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
        }
      }
    }

    for (auto jt : visited) {
      distance[jt->id] = 0;
      predecessor[jt->id] = nullptr;
    }

    RemoveEdges(marked_edges, true);
    num_bubbles += 1 - marked_edges.empty();
  }

  return num_bubbles;
}

std::uint32_t Graph::RemoveLongEdges(std::uint32_t num_rounds) {
  std::uint32_t num_long_edges = 0;

  for (std::uint32_t i = 0; i < num_rounds; ++i) {
    CreateForceDirectedLayout();

    std::unordered_set<std::uint32_t> marked_edges;
    for (const auto& it : nodes_) {
      if (it == nullptr || it->outdegree() < 2) {
        continue;
      }
      for (auto jt : it->outedges) {
        for (auto kt : it->outedges) {
          if (jt != kt && jt->weight * 2.0 < kt->weight) {  // TODO(rvaser)
            marked_edges.emplace(kt->id);
            marked_edges.emplace(kt->pair->id);
          }
        }
      }
    }
    RemoveEdges(marked_edges);
    num_long_edges += marked_edges.size() / 2;

    RemoveTips();
  }

  return num_long_edges;
}

void Graph::CreateForceDirectedLayout(const std::string& path) {
  std::ofstream os;
  bool is_first = true;
  if (!path.empty()) {
    os.open(path);
    os << "{" << std::endl;
  }

  std::vector<std::unordered_set<std::uint32_t>> components;
  std::vector<char> is_visited(nodes_.size(), 0);
  for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
    if (nodes_[i] == nullptr || is_visited[i]) {
      continue;
    }

    components.resize(components.size() + 1);

    std::deque<std::uint32_t> que = { i };
    while (!que.empty()) {
      std::uint32_t j = que.front();
      que.pop_front();

      if (is_visited[j]) {
        continue;
      }
      const auto& node = nodes_[j];
      is_visited[node->id] = 1;
      is_visited[node->pair->id] = 1;
      components.back().emplace((node->id >> 1) << 1);

      for (auto it : node->inedges) {
        que.emplace_back(it->tail->id);
      }
      for (auto it : node->outedges) {
        que.emplace_back(it->head->id);
      }
    }
  }
  std::vector<char>().swap(is_visited);

  std::sort(components.begin(), components.end(),
      [] (const std::unordered_set<std::uint32_t>& lhs,
          const std::unordered_set<std::uint32_t>& rhs) {
        return lhs.size() > rhs.size();
      });

  static std::uint64_t seed = 21;
  seed <<= 1;

  std::mt19937 generator(seed);
  std::uniform_real_distribution<> distribution(0., 1.);

  struct Point {
    Point() = default;
    Point(double x, double y)
        : x(x),
          y(y) {}

    bool operator==(const Point& other) const {
      return x == other.x && y == other.y;
    }
    Point operator+(const Point& other) const {
      return Point(x + other.x, y + other.y);
    }
    Point& operator+=(const Point& other) {
      x += other.x;
      y += other.y;
      return *this;
    }
    Point operator-(const Point& other) const {
      return Point(x - other.x, y - other.y);
    }
    Point operator*(double c) const {
      return Point(x * c, y * c);
    }
    Point& operator/=(double c) {
      x /= c;
      y /= c;
      return *this;
    }
    double norm() const {
      return sqrt(x * x + y * y);
    }

    double x;
    double y;
  };

  struct Quadtree {
    Quadtree(Point nucleus, double width)
        : nucleus(nucleus),
          width(width),
          center(0, 0),
          mass(0),
          subtrees() {
    }

    bool add(const Point& p) {
      if (nucleus.x - width > p.x || p.x > nucleus.x + width ||
          nucleus.y - width > p.y || p.y > nucleus.y + width) {
        return false;
      }
      ++mass;
      if (mass == 1) {
        center = p;
      } else if (subtrees.empty()) {
        if (center == p) {
          return true;
        }
        double w = width / 2;
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y - w), w);
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y - w), w);
        for (auto& it : subtrees) {
          if (it.add(center)) {
            break;
          }
        }
      }
      for (auto& it : subtrees) {
        if (it.add(p)) {
          break;
        }
      }
      return true;
    }

    void centre() {
      if (subtrees.empty()) {
        return;
      }
      center = Point(0, 0);
      for (auto& it : subtrees) {
        it.centre();
        center += it.center * it.mass;
      }
      center /= mass;
    }

    Point force(const Point& p, double k) const {
      auto delta = p - center;
      auto distance = delta.norm();
      if (width * 2 / distance < 1) {
        return delta * (mass * (k * k) / (distance * distance));
      }
      delta = Point(0, 0);
      for (const auto& it : subtrees) {
        delta += it.force(p, k);
      }
      return delta;
    }

    Point nucleus;
    double width;
    Point center;
    std::uint32_t mass;
    std::vector<Quadtree> subtrees;
  };

  std::uint32_t c = 0;
  for (const auto& component : components) {
    if (component.size() < 6) {
      continue;
    }

    bool has_junctions = false;
    for (const auto& it : component) {
      if (nodes_[it]->is_junction()) {
        has_junctions = true;
        break;
      }
    }
    if (has_junctions == false) {
      continue;
    }

    // update transitive edges
    for (const auto& n : component) {
      std::unordered_set<std::uint32_t> valid;
      for (const auto& m : nodes_[n]->transitive) {
        if (component.find(m) != component.end()) {
          valid.emplace(m);
        }
      }
      nodes_[n]->transitive.swap(valid);
    }

    std::uint32_t num_iterations = 100;
    double k = sqrt(1. / static_cast<double>(component.size()));
    double t = 0.1;
    double dt = t / static_cast<double>(num_iterations + 1);

    std::vector<Point> points(nodes_.size());
    for (const auto& it : component) {
      points[it].x = distribution(generator);
      points[it].y = distribution(generator);
    }

    for (std::uint32_t i = 0; i < num_iterations; ++i) {
      Point x = {0, 0}, y = {0, 0};
      for (const auto& n : component) {
        x.x = std::min(x.x, points[n].x);
        x.y = std::max(x.y, points[n].x);
        y.x = std::min(y.x, points[n].y);
        y.y = std::max(y.y, points[n].y);
      }
      double w = (x.y - x.x) / 2, h = (y.y - y.x) / 2;

      Quadtree tree(Point(x.x + w, y.x + h), std::max(w, h) + 0.01);
      for (const auto& n : component) {
        tree.add(points[n]);
      }
      tree.centre();

      std::vector<std::future<void>> thread_futures;
      std::vector<Point> displacements(nodes_.size(), Point(0, 0));

      auto thread_task = [&](std::uint32_t n) -> void {
        auto displacement = tree.force(points[n], k);
        for (auto e : nodes_[n]->inedges) {
          auto m = (e->tail->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (auto e : nodes_[n]->outedges) {
          auto m = (e->head->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (const auto& m : nodes_[n]->transitive) {
          auto delta = points[n] - points[m];
          auto distance = delta.norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        auto length = displacement.norm();
        if (length < 0.01) {
          length = 0.1;
        }
        displacements[n] = displacement * (t / length);
        return;
      };

      for (const auto& n : component) {
        thread_futures.emplace_back(thread_pool_->Submit(thread_task, n));
      }
      for (const auto& it : thread_futures) {
        it.wait();
      }
      for (const auto& n : component) {
        points[n] += displacements[n];
      }

      t -= dt;
    }

    for (const auto& it : edges_) {
      if (it == nullptr || it->id & 1) {
        continue;
      }
      auto n = (it->tail->id >> 1) << 1;
      auto m = (it->head->id >> 1) << 1;

      if (component.find(n) != component.end() &&
          component.find(m) != component.end()) {
        it->weight = (points[n] - points[m]).norm();
        it->pair->weight = it->weight;
      }
    }

    if (!path.empty()) {
      if (!is_first) {
        os << "," << std::endl;
      }
      is_first = false;

      os << "    \"component_" << c++ << "\": {" << std::endl;

      bool is_first_node = true;
      os << "      \"nodes\": {" << std::endl;
      for (const auto& it : component) {
        if (!is_first_node) {
          os << "," << std::endl;
        }
        is_first_node = false;
        os << "        \"" << it << "\": [";
        os << points[it].x << ", ";
        os << points[it].y << ", ";
        os << (nodes_[it]->is_junction() ? 1 : 0) << ", ";
        os << nodes_[it]->count << "]";
      }
      os << std::endl << "      }," << std::endl;

      bool is_first_edge = true;
      os << "      \"edges\": [" << std::endl;
      for (const auto& it : component) {
        for (auto e : nodes_[it]->inedges) {
          auto o = (e->tail->id >> 1) << 1;
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 0]";
        }
        for (auto e : nodes_[it]->outedges) {
          auto o = (e->head->id >> 1) << 1;
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 0]";
        }
        for (const auto& o : nodes_[it]->transitive) {
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 1]";
        }
      }
      os << std::endl << "      ]" << std::endl;
      os << "    }";
    }
  }

  if (!path.empty()) {
    os << std::endl << "}";
    os << std::endl;
    os.close();
  }
}

void Graph::Polish(
    const std::vector<std::unique_ptr<biosoup::Sequence>>& sequences,
    std::uint8_t match,
    std::uint8_t mismatch,
    std::uint8_t gap,
    std::uint32_t cuda_poa_batches,
    bool cuda_banded_alignment,
    std::uint32_t cuda_alignment_batches,
    std::uint32_t num_rounds) {
  if (sequences.empty() || num_rounds == 0) {
    return;
  }

  auto unitigs = GetUnitigs();
  if (unitigs.empty()) {
    return;
  }

  double q = 0.;
  for (const auto& it : sequences) {
    if (it->quality.empty()) {
      continue;
    }
    double p = 0.;
    for (const auto& jt : it->quality) {
      p += jt - 33;
    }
    q += p / it->quality.size();
  }
  if (q == 0.) {  // when all values equal to '!'
    for (const auto& it : sequences) {
      it->quality.clear();
    }
  } else {
    q /= sequences.size();
  }

  auto polisher = racon::Polisher::Create(
      q, 0.3, 500, true,
      match, mismatch, gap,
      thread_pool_,
      cuda_poa_batches,
      cuda_banded_alignment,
      cuda_alignment_batches);

  while (stage_ < static_cast<std::int32_t>(num_rounds)) {
    polisher->Initialize(unitigs, sequences);
    auto polished = polisher->Polish(false);
    unitigs.swap(polished);

    for (const auto& it : unitigs) {  // store unitigs
      const auto& node = nodes_[std::atoi(&it->name[3])];
      std::size_t tag;
      if ((tag = it->name.rfind(':')) != std::string::npos) {
        if (std::atof(&it->name[tag + 1]) > 0) {
          node->is_polished = true;
          node->data = it->data;
          it->ReverseAndComplement();
          node->pair->data = it->data;
          it->ReverseAndComplement();
        }
      }
    }

    biosoup::Timer timer{};
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Polish] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }
}

std::uint32_t Graph::CreateUnitigs(std::uint32_t epsilon) {
  std::unordered_set<std::uint32_t> marked_edges;
  std::vector<std::shared_ptr<Node>> unitigs;
  std::vector<std::shared_ptr<Edge>> unitig_edges;
  std::vector<std::uint32_t> node_updates(nodes_.size(), 0);
  std::vector<char> is_visited(nodes_.size(), 0);

  for (const auto& it : nodes_) {
    if (it == nullptr || is_visited[it->id] || it->is_junction()) {
      continue;
    }

    std::uint32_t extension = 1;

    bool is_circular = false;
    auto begin = it.get();
    while (!begin->is_junction()) {  // extend left
      is_visited[begin->id] = 1;
      is_visited[begin->pair->id] = 1;
      if (begin->indegree() == 0 ||
          begin->inedges.front()->tail->is_junction()) {
        break;
      }
      begin = begin->inedges.front()->tail;
      ++extension;
      if (begin == it.get()) {
        is_circular = true;
        break;
      }
    }

    auto end = it.get();
    while (!end->is_junction()) {  // extend right
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 ||
          end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      ++extension;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    if (!is_circular && begin == end) {
      continue;
    }
    if (!is_circular && extension < 2 * epsilon + 2) {
      continue;
    }

    if (begin != end) {  // remove nodes near junctions
      for (std::uint32_t i = 0; i < epsilon; ++i) {
        begin = begin->outedges.front()->head;
      }
      for (std::uint32_t i = 0; i < epsilon; ++i) {
        end = end->inedges.front()->tail;
      }
    }

    auto unitig = std::make_shared<Node>(begin, end);
    unitigs.emplace_back(unitig);
    unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair));
    unitig->pair = unitigs.back().get();
    unitig->pair->pair = unitig.get();

    if (begin != end) {  // connect unitig to graph
      if (begin->indegree()) {
        marked_edges.emplace(begin->inedges.front()->id);
        marked_edges.emplace(begin->inedges.front()->pair->id);

        auto edge = std::make_shared<Edge>(
            begin->inedges.front()->tail,
            unitig.get(),
            begin->inedges.front()->length);
        unitig_edges.emplace_back(edge);
        unitig_edges.emplace_back(std::make_shared<Edge>(
            unitig->pair,
            begin->inedges.front()->pair->head,
            begin->inedges.front()->pair->length + unitig->pair->data.size() - begin->pair->data.size()));  // NOLINT
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
      if (end->outdegree()) {
        marked_edges.emplace(end->outedges.front()->id);
        marked_edges.emplace(end->outedges.front()->pair->id);

        auto edge = std::make_shared<Edge>(
            unitig.get(),
            end->outedges.front()->head,
            end->outedges.front()->length + unitig->data.size() - end->data.size());  // NOLINT
        unitig_edges.emplace_back(edge);
        unitig_edges.emplace_back(std::make_shared<Edge>(
            end->outedges.front()->pair->tail,
            unitig->pair,
            end->outedges.front()->pair->length));
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
    }

    auto jt = begin;
    while (true) {
      marked_edges.emplace(jt->outedges.front()->id);
      marked_edges.emplace(jt->outedges.front()->pair->id);

      // update transitive edges
      node_updates[jt->id & ~1UL] = unitig->id;
      unitig->transitive.insert(
         nodes_[jt->id & ~1UL]->transitive.begin(),
         nodes_[jt->id & ~1UL]->transitive.end());

      if ((jt = jt->outedges.front()->head) == end) {
        break;
      }
    }
  }

  nodes_.insert(nodes_.end(), unitigs.begin(), unitigs.end());
  edges_.insert(edges_.end(), unitig_edges.begin(), unitig_edges.end());
  RemoveEdges(marked_edges, true);

  for (const auto& it : nodes_) {  // update transitive edges
    if (it) {
      std::unordered_set<std::uint32_t> valid;
      for (auto jt : it->transitive) {
        valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
      }
      it->transitive.swap(valid);
    }
  }

  return unitigs.size() / 2;
}

std::vector<std::unique_ptr<biosoup::Sequence>> Graph::GetUnitigs(
    bool drop_unpolished) {

  CreateUnitigs();

  biosoup::Sequence::num_objects = 0;

  std::vector<std::unique_ptr<biosoup::Sequence>> dst;
  for (const auto& it : nodes_) {
    if (it == nullptr || it->is_rc() || !it->is_unitig()) {
      continue;
    }
    if (drop_unpolished && !it->is_polished) {
      continue;
    }

    std::string name = it->name +
        " LN:i:" + std::to_string(it->data.size()) +
        " RC:i:" + std::to_string(it->count) +
        " XO:i:" + std::to_string(it->is_circular);

    dst.emplace_back(new biosoup::Sequence(name, it->data));
  }

  return dst;
}

void Graph::RemoveEdges(
    const std::unordered_set<std::uint32_t>& indices,
    bool remove_nodes) {

  auto erase_remove = [] (std::vector<Edge*>& edges, Edge* marked) -> void {
    edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
  };

  std::unordered_set<std::uint32_t> node_indices;
  for (auto i : indices) {
    if (remove_nodes) {
      node_indices.emplace(edges_[i]->tail->id);
      node_indices.emplace(edges_[i]->head->id);
    }
    erase_remove(edges_[i]->tail->outedges, edges_[i].get());
    erase_remove(edges_[i]->head->inedges, edges_[i].get());
  }
  if (remove_nodes) {
    for (auto i : node_indices) {
      if (nodes_[i]->outdegree() == 0 && nodes_[i]->indegree() == 0) {
        nodes_[i].reset();
      }
    }
  }
  for (auto i : indices) {
    edges_[i].reset();
  }
}

std::unordered_set<std::uint32_t> Graph::FindRemovableEdges(
    const std::vector<Node*>& path) {
  if (path.empty()) {
    return std::unordered_set<std::uint32_t>{};
  }

  auto find_edge = [] (Node* tail, Node* head) -> Edge* {
    for (auto it : tail->outedges) {
      if (it->head == head) {
        return it;
      }
    }
    return nullptr;
  };

  // find first node with multiple in edges
  std::int32_t pref = -1;
  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->indegree() > 1) {
      pref = i;
      break;
    }
  }

  // find last node with multiple out edges
  std::int32_t suff = -1;
  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->outdegree() > 1) {
      suff = i;
    }
  }

  std::unordered_set<std::uint32_t> dst;
  if (pref == -1 && suff == -1) {  // remove whole path
    for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
    return dst;
  }

  if (pref != -1 && path[pref]->outdegree() > 1) {  // complex path
    return dst;  // empty
  }
  if (suff != -1 && path[suff]->indegree() > 1) {  // complex path
    return dst;  // empty
  }

  if (pref == -1) {  // remove everything after last suffix node
    for (std::uint32_t i = suff; i < path.size() - 1; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  } else if (suff == -1) {  // remove everything before first prefix node
    for (std::int32_t i = 0; i < pref; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  } else if (suff < pref) {  // remove everything in between
    for (std::int32_t i = suff; i < pref; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  }
  return dst;  // empty
}

void Graph::PrintJSON(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  cereal::JSONOutputArchive archive(os);
  for (const auto& it : piles_) {
    if (it->is_invalid()) {
      continue;
    }
    archive(cereal::make_nvp(std::to_string(it->id()), *(it.get())));
  }
}

void Graph::PrintCSV(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : nodes_) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->data.size()
       << " RC:i:" << it->count
       << ","
       << it->pair->id << " [" << it->pair->id / 2 << "]"
       << " LN:i:" << it->pair->data.size()
       << " RC:i:" << it->pair->count
       << ",0,-"
       << std::endl;
  }
  for (const auto& it : edges_) {
    if (it == nullptr) {
      continue;
    }
    biosoup::Sequence lhs("lhs", it->tail->data.substr(it->length));
    biosoup::Sequence rhs("rhs", it->head->data.substr(0, lhs.data.size()));
    EdlibAlignResult result = edlibAlign(
        lhs.data.c_str(), lhs.data.size(),
        rhs.data.c_str(), rhs.data.size(),
        edlibDefaultAlignConfig());
    double score = 1 - result.editDistance / static_cast<double>(lhs.data.size());  // NOLINT

    os << it->tail->id << " [" << it->tail->id / 2 << "]"
       << " LN:i:" << it->tail->data.size()
       << " RC:i:" << it->tail->count
       << ","
       << it->head->id << " [" << it->head->id / 2 << "]"
       << " LN:i:" << it->head->data.size()
       << " RC:i:" << it->head->count
       << ",1,"
       << it->id << " " << it->length << " " << it->weight << " " << score
       << std::endl;
  }
  for (const auto& it : nodes_) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->data.size()
       << " RC:i:" << it->count
       << ","
       << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->data.size()
       << " RC:i:" << it->count
       << ",1,-"
       << std::endl;
  }
  os.close();
}

void Graph::PrintGFA(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : nodes_) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << "S\t" << it->name
       << "\t"  << it->data
       << "\tLN:i:" << it->data.size()
       << "\tRC:i:" << it->count
       << std::endl;
    if (it->is_circular) {
      os << "L\t" << it->name << "\t" << '+'
         << "\t"  << it->name << "\t" << '+'
         << "\t0M"
         << std::endl;
    }
  }
  for (const auto& it : edges_) {
    if (it == nullptr || it->is_rc()) {
      continue;
    }
    os << "L\t" << it->tail->name << "\t" << (it->tail->is_rc() ? '-' : '+')
       << "\t"  << it->head->name << "\t" << (it->head->is_rc() ? '-' : '+')
       << "\t"  << it->tail->data.size() - it->length << 'M'
       << std::endl;
  }
  os.close();
}

void Graph::Store() const {
  std::ofstream os("raven.cereal");
  try {
    cereal::BinaryOutputArchive archive(os);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[raven::Graph::Store] error: unable to store archive");
  }
}

void Graph::Load() {
  std::ifstream is("raven.cereal");
  try {
    cereal::BinaryInputArchive archive(is);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[raven::Graph::Load] error: unable to load archive");
  }
}

}  // namespace raven
