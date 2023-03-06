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
#include "edlib.h"  // NOLINT
#include "racon/polisher.hpp"


namespace raven {

Graph::Node::Node(const biosoup::NucleicAcid& sequence)
    : id(num_objects++),
      sequence(sequence),
      count(1),
      is_unitig(),
      is_circular(),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {}

Graph::Node::Node(Node* begin, Node* end)
    : id(num_objects++),
      sequence(),
      count(),
      is_unitig(),
      is_circular(begin == end),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {
  std::string data{};

  auto it = begin;
  while (true) {
    data += it->outedges.front()->Label();
    count += it->count;
    if ((it = it->outedges.front()->head) == end) {
      break;
    }
  }
  if (begin != end) {
    data += end->sequence.InflateData();
    count += end->count;
  }

  is_unitig = count > 5 && data.size() > 9999;

  sequence = biosoup::NucleicAcid(
    (is_unitig ? "Utg" : "Ctg") + std::to_string(id & (~1UL)),
    data);
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

Graph::Graph(
    bool weaken,
    bool checkpoints,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
          thread_pool :
          std::make_shared<thread_pool::ThreadPool>(1)),
      stage_(-5),
      checkpoints_(checkpoints),
      accurate_(weaken),
      disagreement_(),
      annotations_(),
      piles_(),
      nodes_(),
      edges_() {}

void Graph::Construct(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,  // NOLINT
    double disagreement,
    unsigned split) {
  disagreement_ = disagreement;
  if (sequences.empty() || stage_ > -4) {
    return;
  }

  std::vector<std::vector<biosoup::Overlap>> overlaps(sequences.size());

  // biosoup::Overlap helper functions
  
  auto edlib_alignment_reverse = [] (const std::string& s) -> std::string {
    std::string rs = "";

    for(int i = 0; i < s.length(); i++){
     
      switch(s[i]){
        case 'D':
          rs += 'I';
          break;
        case 'I':
          rs += 'D';
          break;
        default:
          rs += s[i];
      };
    }

    return rs;

  };

    auto overlap_reverse = [&edlib_alignment_reverse] (const biosoup::Overlap& o) -> biosoup::Overlap {
    return biosoup::Overlap(
        o.rhs_id, o.rhs_begin, o.rhs_end,
        o.lhs_id, o.lhs_begin, o.lhs_end,
        o.score, edlib_alignment_reverse(o.alignment),
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
  annotations_.resize(sequences.size());
  struct base_pile{
    std::uint32_t a;
    std::uint32_t c;
    std::uint32_t g;
    std::uint32_t t;
    std::uint32_t i;
    std::uint32_t d;
  };
 

  std::vector<std::vector<base_pile>> snp_base_pile(sequences.size());
  for (const auto& it : sequences) {
    snp_base_pile[it->id].resize(it->inflated_len);
  }


  // annotations_ helper functions

  auto edlib_wrapper = [&] (
      std::uint32_t i,
      const biosoup::Overlap& it,
      const std::string& lhs,
      const std::string& rhs) -> std::string {
    std::string cigar; 
    EdlibAlignResult result = edlibAlign(
        lhs.c_str(), lhs.size(),
        rhs.c_str(), rhs.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0)); // align lhs and rhs
    if (result.status == EDLIB_STATUS_OK) {
      std::uint32_t lhs_pos = it.lhs_begin;
      std::uint32_t rhs_pos = 0;
      for (int j = 0; j < result.alignmentLength; ++j) {
        switch (result.alignment[j]) {
          case 0:
          case 3: {
            switch (rhs[rhs_pos]) {
              case 'A': ++snp_base_pile[i][lhs_pos].a; break;
              case 'C': ++snp_base_pile[i][lhs_pos].c; break;
              case 'G': ++snp_base_pile[i][lhs_pos].g; break;
              case 'T': ++snp_base_pile[i][lhs_pos].t; break;
              default: break; // if they align
            }
            ++lhs_pos;
            ++rhs_pos;
            break;
          }
          case 1: {
            ++snp_base_pile[i][lhs_pos].i;
            ++lhs_pos;
            break; // insertion on the left hand side
          }
          case 2: {
            ++snp_base_pile[i][lhs_pos].d;
            ++rhs_pos;
            break; // deletion on the left hand side
          }
          default: break;
        }
      }
    }
    cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    edlibFreeAlignResult(result);
    return cigar;
  };

  auto call_snps = [&](std::uint32_t i)-> void{
    std::vector<base_pile> tmp = snp_base_pile[i];

    std::vector<std::uint32_t> cov;
    cov.reserve(tmp.size());

    for (const auto& jt : tmp) {
      cov.emplace_back(jt.a + jt.c + jt.g + jt.t);
      //cov.emplace_back(jt.a + jt.c + jt.g + jt.t + jt.d + jt.i);
    }

    std::nth_element(cov.begin(), cov.begin() + cov.size() / 2, cov.end());
    double m = cov[cov.size() / 2] * 2. / 3.;

    std::size_t j = 0;
    for (const auto& jt : tmp){
      std::vector<double> counts = {
          static_cast<double>(jt.a),
          static_cast<double>(jt.c),
          static_cast<double>(jt.g),
          static_cast<double>(jt.t),
          static_cast<double>(jt.d),
          static_cast<double>(jt.i)
      };

      double sum = std::accumulate(counts.begin(), counts.end(), 0);

      if(use_frequencies){
        for (auto& kt : counts){
          kt /= sum;
        };
      };

      if (sum > m){
        std::size_t variants = 0;
        for(const auto& it : counts){
          if(use_frequencies){
            if(freq_low_th < it && it < freq_high_th){
              ++variants;
            }
          } else{
            if(it > variant_call_th){
              ++variants;
            }
          }
        }
      if (variants > 1) annotations_[i].emplace(j);
      

      };
    ++j;
    };


  };



  auto annotation_extract = [&] (
      std::uint32_t i,
      std::uint32_t begin,
      std::uint32_t end,
      std::uint32_t len,
      bool strand) -> std::unordered_set<std::uint32_t> {
    std::unordered_set<std::uint32_t> dst;
    if (annotations_[i].empty()) {
      return dst;
    }
    for (const auto& it : annotations_[i]) {
      if (begin <= it && it <= end) {
        dst.emplace(strand ? it : len - 1 - it);
      }
    }
    return dst;
  };
  // annotations_ helper functions

  if (stage_ == -5 && checkpoints_) {  // checkpoint test
    Store();
  }

  biosoup::Timer timer{};

  std::uint32_t kmer_len = accurate_ ? 29 : 15;
  ram::MinimizerEngine minimizer_engine{
      thread_pool_,
      kmer_len,
      accurate_ ? 9U : 5U
  };

  if (stage_ == -5) {  // find overlaps and create piles
    for (const auto& it : sequences) {
      piles_.emplace_back(new Pile(it->id, it->inflated_len));
    }
    std::size_t bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
      bytes += sequences[i]->inflated_len;
      if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine.Minimize(
          sequences.begin() + j,
          sequences.begin() + i + 1,
          true);
      minimizer_engine.Filter(0.001);

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
            [&] (std::uint32_t i) -> std::vector<biosoup::Overlap> { // map sequences and fill out the potential snp list
              std::vector<biosoup::Overlap> ovlps = minimizer_engine.Map(sequences[i], true, false, true);
              std::vector<biosoup::Overlap> ovlps_final;

              for(const auto& ovlp : ovlps){
                auto lhs = sequences[i]->InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);
                biosoup::NucleicAcid rhs_{"", sequences[ovlp.rhs_id]->InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin)};

                if(!ovlp.strand) rhs_.ReverseAndComplement();

                auto rhs = rhs_.InflateData();

                biosoup::Overlap ovlp_tmp{ovlp.lhs_id, ovlp.lhs_begin, ovlp.lhs_end, ovlp.rhs_id, ovlp.rhs_begin, ovlp.rhs_end, ovlp.score, edlib_wrapper(i, ovlp, lhs, rhs), ovlp.strand};
                ovlps_final.emplace_back(ovlp_tmp);
              
              };

              call_snps(i);
              return ovlps_final;
            },
            k));

        bytes += sequences[k]->inflated_len;
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

    if(print_snp_data){
      std::ofstream outdata;
      outdata.open("snp_annotations.anno");

      for (std::uint32_t i = 0; i < annotations_.size(); ++i) {
        if (annotations_[i].empty()) {
          continue;
        }
        outdata << i;
        for (const auto& jt : annotations_[i]) {
          outdata << " " << jt;
        }
        outdata << std::endl;
      }

    };

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

    std::vector<std::future<void>> futures;
    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> void {
            std::uint32_t k = 0;
            for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
              if (!overlap_update(overlaps[i][j])) {
                continue;
              }

              const auto& it = overlaps[i][j];

              auto lhs_anno = annotation_extract(
                  it.lhs_id,
                  it.lhs_begin,
                  it.lhs_end,
                  sequences[it.lhs_id]->inflated_len,
                  true);

              auto rhs_anno = annotation_extract(
                  it.rhs_id,
                  it.rhs_begin,
                  it.rhs_end,
                  sequences[it.rhs_id]->inflated_len,
                  it.strand);

              if (!lhs_anno.empty() || !rhs_anno.empty()) {
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

                EdlibAlignResult result = edlibAlign(
                    lhs.c_str(), lhs.size(),
                    rhs.c_str(), rhs.size(),
                    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));  // NOLINT

                if (result.status == EDLIB_STATUS_OK) {
                  std::uint32_t lhs_pos = it.lhs_begin;
                  std::uint32_t rhs_pos = it.strand ?
                      it.rhs_begin :
                      sequences[it.rhs_id]->inflated_len - it.rhs_end;

                  std::uint32_t mismatches = 0;
                  std::uint32_t snps = 0;

                  for (int a = 0; a < result.alignmentLength; ++a) {
                    if (lhs_anno.find(lhs_pos) != lhs_anno.end() ||
                        rhs_anno.find(rhs_pos) != rhs_anno.end()) {
                      ++snps;
                      if (result.alignment[a] == 3) {
                        ++mismatches;
                      }
                    }
                    switch (result.alignment[a]) {
                      case 0:
                      case 3: {
                        ++lhs_pos;
                        ++rhs_pos;
                        break;
                      }
                      case 1: {
                        ++lhs_pos;
                        break;
                      }
                      case 2: {
                        ++rhs_pos;
                        break;
                      }
                      default: break;
                    }
                  }

                  edlibFreeAlignResult(result);

                  if (mismatches / static_cast<double>(snps) > disagreement_) {
                    continue;
                  }
                }
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
        std::vector<std::uint16_t> medians;
        for (const auto& jt : it) {
          medians.emplace_back(piles_[jt]->median());
        }
        std::nth_element(
            medians.begin(),
            medians.begin() + medians.size() / 2,
            medians.end());
        std::uint16_t median = medians[medians.size() / 2];

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
    ++stage_;
    if (checkpoints_) {
      timer.Start();
      Store();
      std::cerr << "[raven::Graph::Construct] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
  }

  if (stage_ == -4) {  // find overlaps and repetitive regions
    std::sort(sequences.begin(), sequences.end(),
        [&] (const std::unique_ptr<biosoup::NucleicAcid>& lhs,
             const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
          return piles_[lhs->id]->is_invalid() <  piles_[rhs->id]->is_invalid() ||  // NOLINT
                (piles_[lhs->id]->is_invalid() == piles_[rhs->id]->is_invalid() && lhs->id < rhs->id);  // NOLINT
        });

    std::vector<std::uint32_t> sequences_map(sequences.size());
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
      sequences_map[sequences[i]->id] = i;
    }

    std::uint32_t s = 0;
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
      if (piles_[sequences[i]->id]->is_invalid()) {
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

      minimizer_engine.Minimize(
          sequences.begin() + j,
          sequences.begin() + i + 1);

      std::cerr << "[raven::Graph::Construct] minimized "
                << j << " - " << i + 1 << " / " << s << " "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
      minimizer_engine.Filter(0.001);
      for (std::uint32_t k = 0; k < i + 1; ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&] (std::uint32_t i) -> std::vector<biosoup::Overlap> {
              std::vector<std::uint32_t> filtered;
              auto dst = minimizer_engine.Map(
                  sequences[i],
                  true,  // avoid equal
                  true,  // avoid symmetric
                  false,  // minhash
                  &filtered);
              piles_[sequences[i]->id]->AddKmers(filtered, kmer_len, sequences[i]); // NOLINT

              std::uint32_t k = 0;
              for (std::uint32_t j = 0; j < dst.size(); ++j) {
                if (!overlap_update(dst[j]) || overlap_length(dst[j]) < 1000) {
                  continue;
                }

                const auto& it = dst[j];

                auto lhs_anno = annotation_extract(
                    it.lhs_id,
                    it.lhs_begin,
                    it.lhs_end,
                    sequences[sequences_map[it.lhs_id]]->inflated_len,
                    true);

                auto rhs_anno = annotation_extract(
                    it.rhs_id,
                    it.rhs_begin,
                    it.rhs_end,
                    sequences[sequences_map[it.rhs_id]]->inflated_len,
                    it.strand);

                if (!lhs_anno.empty() || !rhs_anno.empty()) {
                  auto lhs = sequences[sequences_map[it.lhs_id]]->InflateData(
                      it.lhs_begin,
                      it.lhs_end - it.lhs_begin);

                  auto rhs = sequences[sequences_map[it.rhs_id]]->InflateData(
                      it.rhs_begin,
                      it.rhs_end - it.rhs_begin);
                  if (!it.strand) {
                    biosoup::NucleicAcid rhs_{"", rhs};
                    rhs_.ReverseAndComplement();
                    rhs = rhs_.InflateData();
                  }

                  EdlibAlignResult result = edlibAlign(
                      lhs.c_str(), lhs.size(),
                      rhs.c_str(), rhs.size(),
                      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));  // NOLINT

                  if (result.status == EDLIB_STATUS_OK) {
                    std::uint32_t lhs_pos = it.lhs_begin;
                    std::uint32_t rhs_pos = it.strand ?
                        it.rhs_begin :
                        sequences[sequences_map[it.rhs_id]]->inflated_len - it.rhs_end;  // NOLINT

                    std::uint32_t mismatches = 0;
                    std::uint32_t snps = 0;

                    for (int a = 0; a < result.alignmentLength; ++a) {
                      if (lhs_anno.find(lhs_pos) != lhs_anno.end() ||
                          rhs_anno.find(rhs_pos) != rhs_anno.end()) {
                        ++snps;
                        if (result.alignment[a] == 3) {
                          ++mismatches;
                        }
                      }
                      switch (result.alignment[a]) {
                        case 0:
                        case 3: {
                          ++lhs_pos;
                          ++rhs_pos;
                          break;
                        }
                        case 1: {
                          ++lhs_pos;
                          break;
                        }
                        case 2: {
                          ++rhs_pos;
                          break;
                        }
                        default: break;
                      }
                    }

                    edlibFreeAlignResult(result);

                    if (mismatches / static_cast<double>(snps) > disagreement_) {
                      continue;
                    }
                  }
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
      }
    }

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
        [&] (const std::unique_ptr<biosoup::NucleicAcid>& lhs,
             const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
          return lhs->id < rhs->id;
        });
  }

  if (stage_ == -4) {  // resolve repeat induced overlaps
    timer.Start();

    while (true) {
      auto components = connected_components();
      for (const auto& it : components) {
        std::vector<std::uint16_t> medians;
        for (const auto& jt : it) {
          medians.emplace_back(piles_[jt]->median());
        }
        std::nth_element(
            medians.begin(),
            medians.begin() + medians.size() / 2,
            medians.end());
        std::uint16_t median = medians[medians.size() / 2];

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
    Node::num_objects = 0;
    Edge::num_objects = 0;

    std::vector<std::int32_t> sequence_to_node(piles_.size(), -1);
    for (const auto& it : piles_) {  // create nodes
      if (it->is_invalid()) {
        continue;
      }

      std::unordered_set<std::uint32_t> annotations;
      for (const auto& jt : annotations_[it->id()]) {
        if (it->begin() <= jt && jt < it->end()) {
          annotations.emplace(jt - it->begin());
        }
      }
      annotations_[it->id()].swap(annotations);

      auto sequence = biosoup::NucleicAcid{
          sequences[it->id()]->name,
          sequences[it->id()]->InflateData(it->begin(), it->end() - it->begin())};  // NOLINT
      sequence.id = it->id();

      sequence_to_node[it->id()] = Node::num_objects;

      auto node = std::make_shared<Node>(sequence);
      sequence.ReverseAndComplement();
      nodes_.emplace_back(node);
      nodes_.emplace_back(std::make_shared<Node>(sequence));
      node->pair = nodes_.back().get();
      node->pair->pair = node.get();

      if (it->id() < split) {
        node->color = 1;
        node->pair->color = 1;
      }
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
    ++stage_;
    if (checkpoints_) {
      timer.Start();
      Store();
      std::cerr << "[raven::Graph::Construct] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
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

    PrintGfa("after_transitive.gfa");
  }

  if (stage_ == -3) {  // checkpoint
    ++stage_;
    if (checkpoints_) {
      timer.Start();
      Store();
      std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
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

    PrintGfa("after_bubble.gfa");
  }

  if (stage_ == -2) {  // checkpoint
    ++stage_;
    if (checkpoints_) {
      timer.Start();
      Store();
      std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
  }

  if (stage_ == -1) {  // remove long edges
    timer.Start();

    if (annotations_.empty()) {
      CreateUnitigs(42);  // speed up force directed layout
    }
    RemoveLongEdges(16);

    std::cerr << "[raven::Graph::Assemble] removed long edges "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    PrintGfa("after_force.gfa");

    timer.Start();

    while (true) {
      std::uint32_t num_changes = RemoveTips();
      if (annotations_.empty()) {
        num_changes += RemoveBubbles();
      }
      if (num_changes == 0) {
        break;
      }
    }

    if (annotations_.empty()) {
      SalvagePlasmids();
    } else {
      SalvageHaplotypes();
    }

    timer.Stop();
  }

  if (stage_ == -1) {  // checkpoint
    ++stage_;
    if (checkpoints_) {
      timer.Start();
      Store();
      std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
  }

  std::cerr << "[raven::Graph::Assemble] "
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

  // annotations_ helper functions
  auto annotation_extract = [&] (
      std::uint32_t i,
      std::uint32_t begin,
      std::uint32_t end,
      std::uint32_t len,
      bool strand) -> std::unordered_set<std::uint32_t> {
    std::unordered_set<std::uint32_t> dst;
    if (annotations_[i].empty()) {
      return dst;
    }
    for (const auto& it : annotations_[i]) {
      if (begin <= it && it <= end) {
        dst.emplace(strand ? it : len - 1 - it);
      }
    }
    return dst;
  };
  // annotations_ helper functions

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
  auto path_is_simple = [] (const std::vector<Node*>& path) -> bool {
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
  auto path_sequence = [] (const std::vector<Node*>& path) -> std::string {
    std::string data{};
    for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
      for (auto it : path[i]->outedges) {
        if (it->head == path[i + 1]) {
          data += it->Label();
          break;
        }
      }
    }
    data += path.back()->sequence.InflateData();
    return data;
  };
  auto path_annotation = [&] (const std::vector<Node*>& path)
      -> std::unordered_set<std::uint32_t> {
    std::unordered_set<std::uint32_t> dst;
    std::uint32_t offset = 0;
    for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
      for (auto it : path[i]->outedges) {
        if (it->head == path[i + 1]) {
          const auto& annotations = annotation_extract(
              it->tail->sequence.id,
              0,
              it->length,
              it->tail->sequence.inflated_len,
              !it->tail->is_rc());
          for (const auto& jt : annotations) {
            dst.emplace(offset + jt);
          }
          offset += it->length;
        }
      }
    }
    const auto& annotations = annotation_extract(
        path.back()->sequence.id,
        0,
        path.back()->sequence.inflated_len,
        path.back()->sequence.inflated_len,
        !path.back()->is_rc());
    for (const auto& jt : annotations) {
      dst.emplace(offset + jt);
    }
    return dst;
  };
  auto bubble_pop = [&] (
      const std::vector<Node*>& lhs,
      const std::vector<Node*>& rhs) -> std::unordered_set<std::uint32_t> {
    if (lhs.empty() || rhs.empty()) {
      return std::unordered_set<std::uint32_t>{};
    }

    // check BFS result
    std::unordered_set<Node*> bubble;
    bubble.insert(lhs.begin(), lhs.end());
    bubble.insert(rhs.begin(), rhs.end());
    if (lhs.size() + rhs.size() - 2 != bubble.size()) {
      return std::unordered_set<std::uint32_t>{};
    }
    for (const auto& it : lhs) {
      if (bubble.count(it->pair) != 0) {
        return std::unordered_set<std::uint32_t>{};
      }
    }

    if (!path_is_simple(lhs) || !path_is_simple(rhs)) {  // complex path(s)
      // check poppability
      if (FindRemovableEdges(lhs).empty() && FindRemovableEdges(rhs).empty()) {
        return std::unordered_set<std::uint32_t>{};
      }

      // check similarity
      auto l = path_sequence(lhs);
      auto r = path_sequence(rhs);
      if (std::min(l.size(), r.size()) < std::max(l.size(), r.size()) * 0.8) {
        return std::unordered_set<std::uint32_t>{};
      }

      EdlibAlignResult result = edlibAlign(
          l.c_str(), l.size(),
          r.c_str(), r.size(),
          edlibDefaultAlignConfig());
      double score = 0;
      if (result.status == EDLIB_STATUS_OK) {
        score = 1 - result.editDistance /
            static_cast<double>(std::max(l.size(), r.size()));
        edlibFreeAlignResult(result);
      }
      if (score < 0.8) {
        return std::unordered_set<std::uint32_t>{};
      }
      if (!annotations_.empty()) {
        std::unordered_set<std::uint32_t> marked_edges;
        if (!path_is_simple(lhs)) {
          marked_edges = FindRemovableEdges(lhs);
          if (marked_edges.size() > 2) {
            marked_edges.clear();
          }
        }
        if (marked_edges.empty() && !path_is_simple(rhs)) {
          marked_edges = FindRemovableEdges(rhs);
          if (marked_edges.size() > 2) {
            marked_edges.clear();
          }
        }
        return marked_edges;
      }
    }

    if (!annotations_.empty()) {
      auto la = path_annotation(lhs);
      auto ra = path_annotation(rhs);

      if (!la.empty() && !ra.empty()) {
        auto l = path_sequence(lhs);
        auto r = path_sequence(rhs);

        EdlibAlignResult result = edlibAlign(
            l.c_str(), l.size(),
            r.c_str(), r.size(),
            edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

        if (result.status == EDLIB_STATUS_OK) {
          std::uint32_t lhs_pos = 0;
          std::uint32_t rhs_pos = 0;

          std::uint32_t mismatches = 0;
          std::uint32_t snps = 0;

          for (int a = 0; a < result.alignmentLength; ++a) {
            if (la.find(lhs_pos) != la.end() ||
                ra.find(rhs_pos) != ra.end()) {
              ++snps;
              if (result.alignment[a] == 3) {
                ++mismatches;
              }
            }
            switch (result.alignment[a]) {
              case 0:
              case 3: {
                ++lhs_pos;
                ++rhs_pos;
                break;
              }
              case 1: {
                ++lhs_pos;
                break;
              }
              case 2: {
                ++rhs_pos;
                break;
              }
              default: break;
            }
          }

          edlibFreeAlignResult(result);

          if (mismatches / static_cast<double>(snps) > 0.1) {  // disagreement_) {
            return std::unordered_set<std::uint32_t>{};
          }
        }
      }
    }

    std::uint32_t lhs_count = 0;
    for (auto jt : lhs) {
      lhs_count += jt->count;
    }
    std::uint32_t rhs_count = 0;
    for (auto jt : rhs) {
      rhs_count += jt->count;
    }
    auto marked_edges = FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
    if (marked_edges.empty()) {
      marked_edges = FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
    }
    return marked_edges;
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
        if (distance[jt->id] + kt->length > 5000000) {  // out of reach
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
      marked_edges = bubble_pop(lhs, rhs);
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
    double Norm() const {
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

    bool Add(const Point& p) {
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
          if (it.Add(center)) {
            break;
          }
        }
      }
      for (auto& it : subtrees) {
        if (it.Add(p)) {
          break;
        }
      }
      return true;
    }

    void Centre() {
      if (subtrees.empty()) {
        return;
      }
      center = Point(0, 0);
      for (auto& it : subtrees) {
        it.Centre();
        center += it.center * it.mass;
      }
      center /= mass;
    }

    Point Force(const Point& p, double k) const {
      auto delta = p - center;
      auto distance = delta.Norm();
      if (width * 2 / distance < 1) {
        return delta * (mass * (k * k) / (distance * distance));
      }
      delta = Point(0, 0);
      for (const auto& it : subtrees) {
        delta += it.Force(p, k);
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
        tree.Add(points[n]);
      }
      tree.Centre();

      std::vector<std::future<void>> thread_futures;
      std::vector<Point> displacements(nodes_.size(), Point(0, 0));

      auto thread_task = [&](std::uint32_t n) -> void {
        auto displacement = tree.Force(points[n], k);
        for (auto e : nodes_[n]->inedges) {
          auto m = (e->tail->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (auto e : nodes_[n]->outedges) {
          auto m = (e->head->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (const auto& m : nodes_[n]->transitive) {
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        auto length = displacement.Norm();
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
        it->weight = (points[n] - points[m]).Norm();
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

std::uint32_t Graph::SalvagePlasmids() {
  CreateUnitigs();

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> plasmids;
  for (const auto& it : nodes_) {
    if (!it || it->is_rc() || it->is_unitig || !it->is_circular) {
      continue;
    }
    plasmids.emplace_back(new biosoup::NucleicAcid(it->sequence));
  }
  if (plasmids.empty()) {
    return 0;
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
  for (const auto& it : nodes_) {
    if (!it || it->is_rc() || !it->is_unitig) {
      continue;
    }
    unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
  }
  if (unitigs.empty()) {
    return 0;
  }

  // remove duplicates within plasmids
  std::sort(plasmids.begin(), plasmids.end(),
      [] (const std::unique_ptr<biosoup::NucleicAcid>& lhs,
          const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
        return lhs->inflated_len < rhs->inflated_len;
      });
  for (std::uint32_t i = 0; i < plasmids.size(); ++i) {
    plasmids[i]->id = i;
  }

  ram::MinimizerEngine minimizer_engine{thread_pool_};
  minimizer_engine.Minimize(plasmids.begin(), plasmids.end());
  minimizer_engine.Filter(0.001);
  for (auto& it : plasmids) {
    if (!minimizer_engine.Map(it, true, true).empty()) {
      it.reset();
    }
  }
  plasmids.erase(
      std::remove(plasmids.begin(), plasmids.end(), nullptr),
      plasmids.end());

  // remove duplicates between plasmids and unitigs
  minimizer_engine.Minimize(unitigs.begin(), unitigs.end(), true);
  minimizer_engine.Filter(0.001);
  for (auto& it : plasmids) {
    if (!minimizer_engine.Map(it, false, false).empty()) {
      it.reset();
    }
  }
  plasmids.erase(
      std::remove(plasmids.begin(), plasmids.end(), nullptr),
      plasmids.end());

  // update nodes
  for (const auto& it : plasmids) {
    const auto& node = nodes_[std::atoi(&it->name[3])];
    node->is_unitig = node->pair->is_unitig = true;
    node->sequence.name[0] = node->pair->sequence.name[0] = 'U';
  }

  return plasmids.size();
}

void Graph::SalvageHaplotypes() {
  ram::MinimizerEngine minimizer_engine{thread_pool_};

  while (true) {
    auto num_nodes = nodes_.size();

    CreateUnitigs();
    if (num_nodes == nodes_.size()) {
      break;
    }

    // extend primary
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
    for (const auto& it : nodes_) {
      if (!it || it->is_rc() || !it->is_unitig) {
        continue;
      }
      unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
      unitigs.back()->id = unitigs.size() - 1;
    }
    if (unitigs.empty()) {
      return;
    }

    minimizer_engine.Minimize(unitigs.begin(), unitigs.end());
    minimizer_engine.Filter(0.001);

    std::vector<std::vector<biosoup::Overlap>> overlaps(unitigs.size());
    for (const auto& it : unitigs) {
      for (const auto& jt : minimizer_engine.Map(it, true, true)) {
        std::uint32_t lhs_len = unitigs[jt.lhs_id]->inflated_len;
        std::uint32_t lhs_begin = jt.lhs_begin;
        std::uint32_t lhs_end = jt.lhs_end;

        std::uint32_t rhs_len = unitigs[jt.rhs_id]->inflated_len;
        std::uint32_t rhs_begin = jt.strand ? jt.rhs_begin : rhs_len - jt.rhs_end;  // NOLINT
        std::uint32_t rhs_end = jt.strand ? jt.rhs_end : rhs_len - jt.rhs_begin;

        std::uint32_t overhang =
            std::min(lhs_begin, rhs_begin) +
            std::min(lhs_len - lhs_end, rhs_len - rhs_end);

        if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
            rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
          continue;
        }
        if (lhs_end - lhs_begin < lhs_len * 0.9 &&
            rhs_end - rhs_begin < rhs_len * 0.9) {
          continue;
        }
        if ((lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) ||  // NOLINT
            (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end)) {  // NOLINT
          continue;
        }
        if (lhs_len > rhs_len) {
          overlaps[jt.lhs_id].emplace_back(jt);
        } else {
          overlaps[jt.rhs_id].emplace_back(
              jt.rhs_id, jt.rhs_begin, jt.rhs_end,
              jt.lhs_id, jt.lhs_begin, jt.lhs_end,
              jt.score,
              jt.strand);
        }
      }
    }
    for (auto& it : overlaps) {
      if (it.empty() || it.size() > 2) {
        continue;
      }

      const auto& unitig = unitigs[it.front().lhs_id];
      auto data = unitig->InflateData();

      for (auto& jt : it) {
        if (!jt.strand) {
          unitigs[jt.rhs_id]->ReverseAndComplement();
          auto tmp = jt.rhs_begin;
          jt.rhs_begin = unitigs[jt.rhs_id]->inflated_len - jt.rhs_end;
          jt.rhs_end = unitigs[jt.rhs_id]->inflated_len - tmp;
        }

        if (jt.lhs_begin > jt.rhs_begin) {
          data += unitigs[jt.rhs_id]->InflateData(jt.rhs_end);
        } else {
          data = unitigs[jt.rhs_id]->InflateData(0, jt.rhs_begin) + data;
        }

        if (!jt.strand) {
          unitigs[jt.rhs_id]->ReverseAndComplement();
        }
      }

      const auto& node = nodes_[std::atoi(&unitig->name[3])];

      auto na = biosoup::NucleicAcid(node->sequence.name, data);
      node->sequence = na;

      na.name = node->pair->sequence.name;
      na.ReverseAndComplement();
      node->pair->sequence = na;
    }

    // reconstruct alternative
    overlaps.clear();
    unitigs.clear();
    for (const auto& it : nodes_) {
      if (!it || it->is_rc() || !it->is_unitig) {
        continue;
      }
      unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
      unitigs.back()->id = unitigs.size() - 1;
    }
    overlaps.resize(unitigs.size());

    minimizer_engine.Minimize(unitigs.begin(), unitigs.end());
    minimizer_engine.Filter(0.001);

    for (const auto& it : unitigs) {
      for (const auto& jt : minimizer_engine.Map(it, true, true)) {
        std::uint32_t lhs_len = unitigs[jt.lhs_id]->inflated_len;
        std::uint32_t lhs_begin = jt.lhs_begin;
        std::uint32_t lhs_end = jt.lhs_end;

        std::uint32_t rhs_len = unitigs[jt.rhs_id]->inflated_len;
        std::uint32_t rhs_begin = jt.strand ? jt.rhs_begin : rhs_len - jt.rhs_end;  // NOLINT
        std::uint32_t rhs_end = jt.strand ? jt.rhs_end : rhs_len - jt.rhs_begin;

        std::uint32_t overhang =
            std::min(lhs_begin, rhs_begin) +
            std::min(lhs_len - lhs_end, rhs_len - rhs_end);

        if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
            rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
          continue;
        }
        if (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end) {
          overlaps[jt.lhs_id].emplace_back(jt);
        } else if (lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) {  // NOLINT
          overlaps[jt.rhs_id].emplace_back(
              jt.rhs_id, jt.rhs_begin, jt.rhs_end,
              jt.lhs_id, jt.lhs_begin, jt.lhs_end,
              jt.score,
              jt.strand);
        }
      }
    }
    for (auto& it : overlaps) {
      if (it.empty()) {
        continue;
      }

      const auto& unitig = unitigs[it.front().lhs_id];
      if (nodes_[std::atoi(&unitig->name[3])] == nullptr) {
        continue;
      }

      std::sort(it.begin(), it.end(),
          [] (const biosoup::Overlap& lhs,
              const biosoup::Overlap& rhs) -> bool {
            return lhs.lhs_begin < rhs.lhs_begin;
          });
      it.emplace_back(unitig->id, -1, -1, unitig->id, -1, -1, 0);  // dummy

      std::string data = unitig->InflateData(0, it.front().lhs_begin);

      std::unordered_set<std::uint32_t> marked_edges;
      for (std::uint32_t j = 0; j < it.size() - 1; ++j) {
        auto& jt = it[j];

        const auto& n = nodes_[std::atoi(&unitigs[jt.rhs_id]->name[3])];
        if (n == nullptr) {
          continue;
        }

        if (!jt.strand) {
          unitigs[jt.rhs_id]->ReverseAndComplement();
          auto tmp = jt.rhs_begin;
          jt.rhs_begin = unitigs[jt.rhs_id]->inflated_len - jt.rhs_end;
          jt.rhs_end = unitigs[jt.rhs_id]->inflated_len - tmp;
        }

        data += unitigs[jt.rhs_id]->InflateData(
            jt.rhs_begin,
            jt.rhs_end - jt.rhs_begin);
        if (jt.lhs_end < it[j + 1].lhs_begin) {
          data += unitig->InflateData(
              jt.lhs_end,
              it[j + 1].lhs_begin - jt.lhs_end);
        }

        for (const auto& kt : n->inedges) {
          marked_edges.emplace(kt->id);
          marked_edges.emplace(kt->pair->id);
        }
        for (const auto& kt : n->outedges) {
          marked_edges.emplace(kt->id);
          marked_edges.emplace(kt->pair->id);
        }

        if (!jt.strand) {
          unitigs[jt.rhs_id]->ReverseAndComplement();
        }
      }
      RemoveEdges(marked_edges);

      it.pop_back();
      for (const auto& jt : it) {
        auto id = std::atoi(&unitigs[jt.rhs_id]->name[3]);
        nodes_[id ^ 1].reset();
        nodes_[id].reset();
      }

      auto na = biosoup::NucleicAcid("", data);

      auto node = std::make_shared<Node>(na);
      node->sequence.name = "Utg" + std::to_string(node->id & (~1UL));
      node->is_unitig = true;
      nodes_.emplace_back(node);

      na.ReverseAndComplement();
      nodes_.emplace_back(std::make_shared<Node>(na));

      node->pair = nodes_.back().get();
      node->pair->pair = node.get();
      node->pair->sequence.name = node->sequence.name;
      node->pair->is_unitig = true;
    }
  }
}

void Graph::Polish(
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
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

  piles_.clear();  // save memory

  double avg_q = 0.;
  for (const auto& it : sequences) {
    if (it->block_quality.empty()) {
      continue;
    }
    double q = std::accumulate(
        it->block_quality.begin(),
        it->block_quality.end(),
        0.);
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
      thread_pool_,
      avg_q,
      0.3,
      500,
      true,
      match, mismatch, gap,
      cuda_poa_batches,
      cuda_banded_alignment,
      cuda_alignment_batches);

  while (stage_ < static_cast<std::int32_t>(num_rounds)) {
    auto polished = polisher->Polish(unitigs, sequences, false);
    unitigs.swap(polished);

    for (const auto& it : unitigs) {  // store unitigs
      const auto& node = nodes_[std::atoi(&it->name[3])];
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
          node->sequence.deflated_data = node->pair->sequence.deflated_data = it->deflated_data;  // NOLINT
          node->sequence.inflated_len = node->pair->sequence.inflated_len = it->inflated_len;  // NOLINT
        }
      }
    }

    ++stage_;
    if (checkpoints_) {
      biosoup::Timer timer{};
      timer.Start();
      Store();
      std::cerr << "[raven::Graph::Polish] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
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
            begin->inedges.front()->pair->length + unitig->pair->sequence.inflated_len - begin->pair->sequence.inflated_len));  // NOLINT
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
      if (end->outdegree()) {
        marked_edges.emplace(end->outedges.front()->id);
        marked_edges.emplace(end->outedges.front()->pair->id);

        auto edge = std::make_shared<Edge>(
            unitig.get(),
            end->outedges.front()->head,
            end->outedges.front()->length + unitig->sequence.inflated_len - end->sequence.inflated_len);  // NOLINT
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

std::vector<std::unique_ptr<biosoup::NucleicAcid>> Graph::GetUnitigs(
    bool drop_unpolished) {

  CreateUnitigs();

  biosoup::NucleicAcid::num_objects = 0;

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
  for (const auto& it : nodes_) {
    if (it == nullptr || it->is_rc() || !it->is_unitig) {
      continue;
    }
    if (drop_unpolished && !it->is_polished) {
      continue;
    }

    std::string name = it->sequence.name +
        " LN:i:" + std::to_string(it->sequence.inflated_len) +
        " RC:i:" + std::to_string(it->count) +
        " XO:i:" + std::to_string(it->is_circular);

    dst.emplace_back(new biosoup::NucleicAcid(
        name,
        it->sequence.InflateData()));
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

void Graph::PrintJson(const std::string& path) const {
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

void Graph::PrintCsv(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : nodes_) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << it->id << " [" << it->sequence.id << "]"
       << " LN:i:" << it->sequence.inflated_len
       << " RC:i:" << it->count
       << ","
       << it->pair->id << " [" << it->pair->sequence.id << "]"
       << " LN:i:" << it->pair->sequence.inflated_len
       << " RC:i:" << it->pair->count
       << ",0,-"
       << std::endl;
  }
  for (const auto& it : edges_) {
    if (it == nullptr) {
      continue;
    }
    os << it->tail->id << " [" << it->tail->sequence.id << "]"
       << " LN:i:" << it->tail->sequence.inflated_len
       << " RC:i:" << it->tail->count
       << ","
       << it->head->id << " [" << it->head->sequence.id << "]"
       << " LN:i:" << it->head->sequence.inflated_len
       << " RC:i:" << it->head->count
       << ",1,"
       << it->id << " " << it->length << " " << it->weight
       << std::endl;
  }
  for (const auto& it : nodes_) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }
    os << it->id << " [" << it->sequence.id << "]"
       << " LN:i:" << it->sequence.inflated_len
       << " RC:i:" << it->count
       << ","
       << it->id << " [" << it->sequence.id << "]"
       << " LN:i:" << it->sequence.inflated_len
       << " RC:i:" << it->count
       << ",1,-"
       << std::endl;
  }
  os.close();
}

void Graph::PrintGfa(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : nodes_) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << "S\t" << it->sequence.name
       << "\t"  << "*"  // it->sequence.InflateData()
       << "\tLN:i:" << it->sequence.inflated_len
       << "\tRC:i:" << it->count
       << "\tCL:z:" << (it->color ? "blue" : "orange")
       << std::endl;
    if (it->is_circular) {
      os << "L\t" << it->sequence.name << "\t" << '+'
         << "\t"  << it->sequence.name << "\t" << '+'
         << "\t0M"
         << std::endl;
    }
  }
  for (const auto& it : edges_) {
    if (it == nullptr || it->is_rc()) {
      continue;
    }
    os << "L\t" << it->tail->sequence.name << "\t" << (it->tail->is_rc() ? '-' : '+')  // NOLINT
       << "\t"  << it->head->sequence.name << "\t" << (it->head->is_rc() ? '-' : '+')  // NOLINT
       << "\t"  << it->tail->sequence.inflated_len - it->length << 'M'
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
