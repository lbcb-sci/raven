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

#include "ram/overlap.hpp"

#include "pile.hpp"

namespace raven {

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
    return std::unique_ptr<Pile>(new Pile(id, size));
}

Pile::Pile(std::uint32_t id, std::uint32_t size)
        : id_(id), data_(size, 0), begin_(0), end_(size), median_(0), state_(0),
        chimeric_regions_() {
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

    std::int32_t w = 847;
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
                                dst.emplace_back(first_up, last_up);
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
                    dst.emplace_back(first_up, last_up);
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
                                dst.emplace_back(first_down, last_down);
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
                    dst.emplace_back(first_down, last_down);
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

    if (begin >= end || end - begin < 1260) {
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

void Pile::add_layers(std::vector<ram::Overlap>::const_iterator begin,
    std::vector<ram::Overlap>::const_iterator end) {

    if (begin >= end) {
        return;
    }

    std::vector<std::uint32_t> boundaries;
    for (auto it = begin; it != end; ++it) {
        if (it->q_id == id_) {
            boundaries.emplace_back((it->q_begin + 1) << 1);
            boundaries.emplace_back((it->q_end - 1) << 1 | 1);
        } else if (it->t_id == id_) {
            boundaries.emplace_back((it->t_begin + 1) << 1);
            boundaries.emplace_back((it->t_end - 1) << 1 | 1);
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

    ss << "\"median\":" << median_;
    ss << "}";

    return ss.str();
}

}
