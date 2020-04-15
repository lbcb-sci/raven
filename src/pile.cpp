/*!
 * @file pile.cpp
 *
 * @brief Pile class source file
 */

#include <algorithm>
#include <sstream>
#include <exception>
#include <deque>

#include <iostream>

#include "pile.hpp"

namespace raven {

constexpr std::uint32_t kS = 4; // reduce memory 2^kS times

void mergeRegions(std::vector<std::pair<std::uint32_t, std::uint32_t>>& regions) {

    std::vector<std::pair<std::uint32_t, std::uint32_t>> dst;
    std::vector<char> is_merged(regions.size(), 0);
    for (std::uint32_t i = 0; i < regions.size(); ++i) {
        if (is_merged[i]) {
            continue;
        }
        for (std::uint32_t j = 0; j < regions.size(); ++j) {
            if (i != j && !is_merged[j] &&
                regions[i].first < regions[j].second &&
                regions[i].second > regions[j].first) {

                is_merged[j] = true;
                regions[i].first = std::min(regions[i].first, regions[j].first);
                regions[i].second = std::max(regions[i].second, regions[j].second);
            }
        }
        dst.emplace_back(regions[i].first, regions[i].second);
    }
    regions.swap(dst);
}

std::unique_ptr<Pile> createPile(std::uint32_t id, std::uint32_t size) {
    return std::unique_ptr<Pile>(new Pile(id, size >> kS));
}

Pile::Pile(std::uint32_t id, std::uint32_t size)
        : id_(id), data_(size, 0), begin_(0), end_(size), median_(0), state_(0),
        chimeric_regions_(), repetitive_regions_() {
}

std::uint32_t Pile::begin() const {
    return begin_ << kS;
}

std::uint32_t Pile::end() const {
    return end_ << kS;
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> Pile::find_slopes(double q) {

    using Subpile = std::deque<std::pair<std::int32_t, std::int32_t>>;

    auto subpile_add = [] (Subpile& s, std::int32_t value, std::int32_t position) -> void {
        while (!s.empty() && s.back().second <= value) {
            s.pop_back();
        }
        s.emplace_back(position, value);
    };
    auto subpile_update = [] (Subpile& s, std::int32_t position) {
        while (!s.empty() && s.front().first <= position) {
            s.pop_front();
        }
    };

    // find slopes
    std::vector<std::pair<std::uint32_t, std::uint32_t>> dst;

    std::int32_t w = 847 >> kS;
    std::int32_t data_size = data_.size();

    Subpile left_subpile;
    std::uint32_t first_down = 0, last_down = 0;
    bool found_down = false;

    Subpile right_subpile;
    std::uint32_t first_up = 0, last_up = 0;
    bool found_up = false;

    // find slope regions
    for (std::int32_t i = 0; i < w; ++i) {
        subpile_add(right_subpile, data_[i], i);
    }
    for (std::int32_t i = 0; i < data_size; ++i) {
        if (i > 0) {
            subpile_add(left_subpile, data_[i - 1], i - 1);
        }
        subpile_update(left_subpile, i - 1 - w);

        if (i < data_size - w) {
            subpile_add(right_subpile, data_[i + w], i + w);
        }
        subpile_update(right_subpile, i);

        std::int32_t c = data_[i] * q;
        if (i != 0 && left_subpile.front().second > c) {
            if (found_down) {
                if (i - last_down > 1) {
                    dst.emplace_back(first_down << 1 | 0, last_down);
                    first_down = i;
                }
            } else {
                found_down = true;
                first_down = i;
            }
            last_down = i;
        }
        if (i != (data_size - 1) && right_subpile.front().second > c) {
            if (found_up) {
                if (i - last_up > 1) {
                    dst.emplace_back(first_up << 1 | 1, last_up);
                    first_up = i;
                }
            } else {
                found_up = true;
                first_up = i;
            }
            last_up = i;
        }
    }
    if (found_down) {
        dst.emplace_back(first_down << 1 | 0, last_down);
    }
    if (found_up) {
        dst.emplace_back(first_up << 1 | 1, last_up);
    }

    if (dst.empty()) {
        return dst;
    }

    // separate overlaping slopes
    while (true) {
        std::sort(dst.begin(), dst.end());

        bool is_changed = false;
        for (std::uint32_t i = 0; i < dst.size() - 1; ++i) {
            if (dst[i].second < (dst[i + 1].first >> 1)) {
                continue;
            }

            if (dst[i].first & 1) {
                right_subpile.clear();
                found_up = false;

                std::uint32_t subpile_begin = dst[i].first >> 1;
                std::uint32_t subpile_end = std::min(dst[i].second, dst[i + 1].second);

                for (std::uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
                    subpile_add(right_subpile, data_[j], j);
                }
                for (std::uint32_t j = subpile_begin; j < subpile_end; ++j) {
                    subpile_update(right_subpile, j);
                    if (data_[j] * q < right_subpile.front().second) {
                        if (found_up) {
                            if (j - last_up > 1) {
                                dst.emplace_back(first_up << 1 | 1, last_up);
                                first_up = j;
                            }
                        } else {
                            found_up = true;
                            first_up = j;
                        }
                        last_up = j;
                    }
                }
                if (found_up) {
                    dst.emplace_back(first_up << 1 | 1, last_up);
                }
                dst[i].first = subpile_end << 1 | 1;

            } else {
                if (dst[i].second == (dst[i + 1].first >> 1)) {
                    continue;
                }

                left_subpile.clear();
                found_down = false;

                std::uint32_t subpile_begin = std::max(dst[i].first >> 1,
                    dst[i + 1].first >> 1);
                std::uint32_t subpile_end = dst[i].second;

                for (std::uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
                    if (!left_subpile.empty() && data_[j] * q < left_subpile.front().second) {
                        if (found_down) {
                            if (j - last_down > 1) {
                                dst.emplace_back(first_down << 1, last_down);
                                first_down = j;
                            }
                        } else {
                            found_down = true;
                            first_down = j;
                        }
                        last_down = j;
                    }
                    subpile_add(left_subpile, data_[j], j);
                }
                if (found_down) {
                    dst.emplace_back(first_down << 1, last_down);
                }
                dst[i].second = subpile_begin;
            }

            is_changed = true;
            break;
        }

        if (!is_changed) {
            break;
        }
    }

    // narrow slopes
    for (std::uint32_t i = 0; i < dst.size() - 1; ++i) {
        if ((dst[i].first & 1) && !(dst[i + 1].first & 1)) {

            std::uint32_t subpile_begin = dst[i].second;
            std::uint32_t subpile_end = dst[i + 1].first >> 1;

            if (subpile_end - subpile_begin > static_cast<std::uint32_t>(w)) {
                continue;
            }

            std::uint32_t max_coverage = 0;
            for (std::uint32_t j = subpile_begin + 1; j < subpile_end; ++j) {
                max_coverage = std::max(max_coverage, data_[j]);
            }

            std::uint32_t valid_point = dst[i].first >> 1;
            for (std::uint32_t j = dst[i].first >> 1; j <= subpile_begin; ++j) {
                if (max_coverage > data_[j] * q) {
                    valid_point = j;
                }
            }
            dst[i].second = valid_point;

            valid_point = dst[i + 1].second;
            for (uint32_t j = subpile_end; j <= dst[i + 1].second; ++j) {
                if (max_coverage > data_[j] * q) {
                    valid_point = j;
                    break;
                }
            }
            dst[i + 1].first = valid_point << 1 | 0;
        }
    }

    return dst;
}

void Pile::find_valid_region(std::uint32_t coverage) {

    std::uint32_t begin = 0, end = 0;
    for (std::uint32_t i = begin_; i < end_; ++i) {
        if (data_[i] < coverage) {
            continue;
        }
        for (std::uint32_t j = i + 1; j < end_; ++j) {
            if (data_[j] >= coverage) {
                continue;
            }
            if (end - begin < j - i) {
                begin = i;
                end = j;
            }
            i = j;
            break;
        }
    }

    set_valid_region(begin, end);
}

void Pile::set_valid_region(std::uint32_t begin, std::uint32_t end) {

    if (begin >= end || end - begin < 1260 >> kS) {
        set_invalid();
        return;
    }

    for (std::uint32_t i = begin_; i < begin; ++i) {
        data_[i] = 0;
    }
    for (std::uint32_t i = end; i < end_; ++i) {
        data_[i] = 0;
    }

    begin_ = begin;
    end_ = end;
}

void Pile::clear_valid_region() {
    std::fill(data_.begin() + begin_, data_.begin() + end_, 0);
}

void Pile::clear_invalid_region() {
    std::fill(data_.begin(), data_.begin() + begin_, 0);
    std::fill(data_.begin() + end_, data_.end(), 0);
}

void Pile::find_median() {

    std::vector<std::uint32_t> data(data_.begin() + begin_, data_.begin() + end_);
    std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
    median_ = data[data.size() / 2];
}

void Pile::add_layers(
    std::vector<biosoup::Overlap>::const_iterator begin,
    std::vector<biosoup::Overlap>::const_iterator end) {

    if (begin >= end) {
        return;
    }

    std::vector<std::uint32_t> boundaries;
    for (auto it = begin; it != end; ++it) {
        if (it->lhs_id == id_) {
            boundaries.emplace_back(((it->lhs_begin >> kS) + 1) << 1);
            boundaries.emplace_back(((it->lhs_end >> kS) - 1) << 1 | 1);
        } else if (it->rhs_id == id_) {
            boundaries.emplace_back(((it->rhs_begin >> kS) + 1) << 1);
            boundaries.emplace_back(((it->rhs_end >> kS) - 1) << 1 | 1);
        }
    }
    std::sort(boundaries.begin(), boundaries.end());

    std::uint32_t coverage = 0;
    std::uint32_t last_boundary = 0;
    for (const auto& it: boundaries) {
        if (coverage > 0) {
            for (std::uint32_t i = last_boundary; i < (it >> 1); ++i) {
                data_[i] += coverage;
            }
        }
        last_boundary = it >> 1;
        coverage += it & 1 ? -1 : 1;
    }
}

void Pile::find_chimeric_regions() {

    auto slopes = find_slopes(1.82);
    if (slopes.empty()) {
        return;
    }

    for (std::uint32_t i = 0; i < slopes.size() - 1; ++i) {
        if (!(slopes[i].first & 1) && (slopes[i + 1].first & 1)) {
            chimeric_regions_.emplace_back(slopes[i].first >> 1,
                slopes[i + 1].second);
        }
    }

    mergeRegions(chimeric_regions_);
}

void Pile::resolve_chimeric_regions(std::uint32_t median) {

    auto is_chimeric_region = [&] (const std::pair<std::uint32_t, std::uint32_t>& r) -> bool {
        for (std::uint32_t i = r.first; i <= r.second; ++i) {
            if (data_[i] * 1.84 <= median) {
                return true;
            }
        }
        return false;
    };

    std::uint32_t begin = 0, end = 0, last = begin_;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> unresolved_regions;

    for (const auto& it: chimeric_regions_) {
        if (begin_ > it.first || end_ < it.second) {
            continue;
        }
        if (is_chimeric_region(it)) {
            if (it.first - last > end - begin) {
                begin = last;
                end = it.first;
            }
            last = it.second;
        } else {
            unresolved_regions.emplace_back(it);
        }
    }
    if (end_ - last > end - begin) {
        begin = last;
        end = end_;
    }

    if (begin != begin_ || end != end_) {
        state_ |= (1U << 2);
    }

    chimeric_regions_.swap(unresolved_regions);
    set_valid_region(begin, end);
}

void Pile::find_repetitive_regions(std::uint32_t median) {

    auto slopes = find_slopes(1.42);
    if (slopes.empty()) {
        return;
    }

    auto is_repetitive_region = [&] (
        const std::pair<std::uint32_t, std::uint32_t>& begin,
        const std::pair<std::uint32_t, std::uint32_t>& end) -> bool {

        if (((end.first >> 1) + end.second) / 2 -
            ((begin.first >> 1) + begin.second) / 2 > 0.84 * (end_ - begin_)) {
            return false;
        }
        bool found_peak = false;
        std::uint32_t peak_value = 1.42 * std::max(data_[begin.second], data_[end.first >> 1]);
        std::uint32_t min_value = median * 1.42;
        std::uint32_t num_valid = 0;

        for (std::uint32_t i = begin.second + 1; i < (end.first >> 1); ++i) {
            if (data_[i] > min_value) {
                ++num_valid;
            }
            if (data_[i] > peak_value) {
                found_peak = true;
            }
        }

        if (!found_peak || num_valid < 0.9 * ((end.first >> 1) - begin.second)) {
            return false;
        }
        return true;
    };

    for (std::uint32_t i = 0; i < slopes.size() - 1; ++i) {
        if (!(slopes[i].first & 1)) {
            continue;
        }
        for (std::uint32_t j = i + 1; j < slopes.size(); ++j) {
            if (slopes[j].first & 1) {
                continue;
            }

            if (is_repetitive_region(slopes[i], slopes[j])) {
                repetitive_regions_.emplace_back(
                    slopes[i].second - 0.336 * (slopes[i].second - (slopes[i].first >> 1)),
                    (slopes[j].first >> 1) + 0.336 * (slopes[j].second - (slopes[j].first >> 1)));
            }
        }
    }

    mergeRegions(repetitive_regions_);
    for (auto& it: repetitive_regions_) {
        it.first = std::max(begin_, it.first) << 1;
        it.second = std::min(end_, it.second);
    }
}

void Pile::resolve_repetitive_regions(const biosoup::Overlap& o) {

    if (repetitive_regions_.empty()) {
        return;
    }
    if (id_ != o.lhs_id && id_ != o.rhs_id) {
        return;
    }

    std::uint32_t begin = (id_ == o.lhs_id ? o.lhs_begin : o.rhs_begin) >> kS;
    std::uint32_t end = (id_ == o.lhs_id ? o.lhs_end : o.rhs_end) >> kS;
    std::uint32_t fuzz = 420 >> kS;
    std::uint32_t offset = 0.1 * (end_ - begin_);

    for (auto& it: repetitive_regions_) {
        if (begin < it.second && (it.first >> 1) < end) {
            if ((it.first >> 1) < begin_ + offset && begin - begin_ < end_ - end) {
                if (end >= it.second + fuzz) {
                    it.first |= 1;
                }
            } else if (it.second > end_ - offset && begin - begin_ > end_ - end) {
                if (begin + fuzz <= (it.first >> 1)) {
                    it.first |= 1;
                }
            }
        }
    }
}

bool Pile::is_false_overlap(const biosoup::Overlap& o) {

    if (repetitive_regions_.empty()) {
        return false;
    }
    if (id_ != o.lhs_id && id_ != o.rhs_id) {
        return false;
    }

    std::uint32_t begin = (id_ == o.lhs_id ? o.lhs_begin : o.rhs_begin) >> kS;
    std::uint32_t end = (id_ == o.lhs_id ? o.lhs_end : o.rhs_end) >> kS;
    std::uint32_t fuzz = 420 >> kS;
    std::uint32_t offset = 0.1 * (end_ - begin_);

    for (const auto& it: repetitive_regions_) {
        if (begin < it.second && (it.first >> 1) < end) {
            if ((it.first >> 1) < begin_ + offset) {
                if (end < it.second + fuzz && (it.first & 1)) {
                    return true;
                }
            } else if (it.second > end_ - offset) {
                if (begin + fuzz > (it.first >> 1) && (it.first & 1)) {
                    return true;
                }
            }
        }
    }

    return false;
}

std::string Pile::to_json() const {

    std::stringstream ss;
    ss << "\"" << id_ << "\":{";

    ss << "\"data\":[";
    for (std::uint32_t i = 0; i < data_.size() - 1; ++i) {
        ss << data_[i] << ",";
    }
    ss << data_.back() << "],";

    ss << "\"begin\":" << begin_ << ",";
    ss << "\"end\":" << end_ << ",";

    ss << "\"chimeric\":[";
    for (std::uint32_t i = 0; i < chimeric_regions_.size(); ++i) {
        ss << chimeric_regions_[i].first << "," << chimeric_regions_[i].second;
        if (i < chimeric_regions_.size() - 1) {
            ss << ",";
        }
    }
    ss << "],";

    ss << "\"repetitive\":[";
    for (std::uint32_t i = 0; i < repetitive_regions_.size(); ++i) {
        ss << (repetitive_regions_[i].first >> 1) << "," << repetitive_regions_[i].second;
        if (i < repetitive_regions_.size() - 1) {
            ss << ",";
        }
    }
    ss << "],";

    ss << "\"median\":" << median_;
    ss << "}";

    return ss.str();
}

}
