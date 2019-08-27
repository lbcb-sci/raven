/*!
 * @file pile.cpp
 *
 * @brief Pile class source file
 */

#include <algorithm>
#include <sstream>
#include <exception>
#include <deque>

#include "ram/overlap.hpp"

#include "pile.hpp"

namespace raven {

std::unique_ptr<Pile> createPile(std::uint32_t id, std::uint32_t size) {
    return std::unique_ptr<Pile>(new Pile(id, size));
}

Pile::Pile(std::uint32_t id, std::uint32_t size)
        : id_(id), data_(size, 0), begin_(0), end_(size), median_(0), state_(0),
        chimeric_regions_() {
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> Pile::find_slopes(double q) {

    // find slopes

    // separate overlaping slopes

    // narrow slopes
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
        state_ |= 1U;
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
        if (it->q_id != id_) {
            continue;
        }
        boundaries.emplace_back((it->q_begin + 1) << 1);
        boundaries.emplace_back((it->q_end - 1) << 1 | 1);
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

    std::vector<std::uint32_t>().swap(boundaries);
}

void Pile::find_chimeric_regions() {

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

    ss << "\"median\":" << median_;
    ss << "}";

    return ss.str();
}

}
