// Copyright (c) 2020 Robert Vaser

#ifndef RAVEN_PILE_HPP_
#define RAVEN_PILE_HPP_

#include <cstdint>
#include <vector>
#include <utility>

#include "biosoup/overlap.hpp"
#include "cereal/cereal.hpp"
#include "cereal/access.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"

namespace raven {

constexpr std::uint32_t kPSS = 4;  // shrink 2 ^ kPSS times

class Pile {
 public:
  Pile(std::uint32_t id, std::uint32_t len);

  Pile(const Pile&) = delete;
  Pile& operator=(const Pile&) = delete;

  Pile(Pile&&) = default;
  Pile& operator=(Pile&&) = default;

  ~Pile() = default;

  std::uint32_t id() const {
    return id_;
  }

  std::uint32_t begin() const {
    return begin_ << kPSS;
  }

  std::uint32_t end() const {
    return end_ << kPSS;
  }

  std::uint32_t length() const {
    return end() - begin();
  }

  std::uint32_t median() const {
    return median_;
  }

  bool is_invalid() const {
    return is_invalid_;
  }

  void set_is_invalid() {
    is_invalid_ = true;
  }

  bool is_contained() const {
    return is_contained_;
  }

  void set_is_contained() {
    is_contained_ = true;
  }

  bool is_chimeric() const {
    return is_chimeric_;
  }

  bool is_maybe_chimeric() const {
    return !chimeric_regions_.empty();
  }

  void set_is_chimeric() {
    is_chimeric_ = true;
  }

  bool is_repetitive() const {
    return is_repetitive_;
  }

  void set_is_repetitive() {
    is_repetitive_ = true;
  }

  // add coverage
  void AddLayers(
      std::vector<biosoup::Overlap>::const_iterator begin,
      std::vector<biosoup::Overlap>::const_iterator end);

  // store longest region with values greater or equal than given coverage
  void FindValidRegion(std::uint32_t coverage);

  // fill valid region with zeroes
  void ClearValidRegion();

  // fill everything outside valid region with zeroes
  void ClearInvalidRegion();

  // store median of valid region
  void FindMedian();

  // store coverage drops
  void FindChimericRegions();

  void RemoveChimericNotations() {
    chimeric_regions_.clear();
  }

  // update valid region to longest non-chimeric given the component median
  void ClearChimericRegions(std::uint32_t median);

  // store coverage spikes given component median
  void FindRepetitiveRegions(std::uint32_t median);

  // increase confidence in repetitive regions given an overlap
  void UpdateRepetitiveRegions(const biosoup::Overlap& o);

  // define relationship between repetitive regions and given overlap
  bool CheckRepetitiveRegions(const biosoup::Overlap& o);

  // remove all repetitive regions
  void ClearRepetitiveRegions();

 private:
  Pile() = default;

  template<class Archive>
  void serialize(Archive& archive) {  // NOLINT
    archive(
        CEREAL_NVP(id_),
        CEREAL_NVP(begin_),
        CEREAL_NVP(end_),
        CEREAL_NVP(median_),
        CEREAL_NVP(is_invalid_),
        CEREAL_NVP(is_contained_),
        CEREAL_NVP(is_chimeric_),
        CEREAL_NVP(is_repetitive_),
        CEREAL_NVP(data_),
        CEREAL_NVP(chimeric_regions_),
        CEREAL_NVP(repetitive_regions_));
  }

  friend cereal::access;

  using Region = std::pair<std::uint32_t, std::uint32_t>;

  // clear invalid region after update
  void UpdateValidRegion(std::uint32_t begin, std::uint32_t end);

  // merge overlapping regions
  static std::vector<Region> MergeRegions(const std::vector<Region>& regions);

  // find drop and spike regions
  std::vector<Region> FindSlopes(double q);

  std::uint32_t id_;
  std::uint32_t begin_;
  std::uint32_t end_;
  std::uint32_t median_;
  bool is_invalid_;
  bool is_contained_;
  bool is_chimeric_;
  bool is_repetitive_;
  std::vector<std::uint32_t> data_;
  std::vector<Region> chimeric_regions_;
  std::vector<Region> repetitive_regions_;
};

}  // namespace raven

#endif  // RAVEN_PILE_HPP_
