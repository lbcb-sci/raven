/*!
 * @file pile.cpp
 *
 * @brief Pile class source file
 */

#include <algorithm>
#include <sstream>
#include <exception>

#include "ram/overlap.hpp"

#include "pile.hpp"

namespace raven {

std::unique_ptr<Pile> createPile(std::uint32_t id, std::uint32_t size) {
    return std::unique_ptr<Pile>(new Pile(id, size));
}

Pile::Pile(std::uint32_t id, std::uint32_t size)
        : id_(id), data_(size, 0), begin_(0), end_(size), median_(0) {
}

void Pile::add_layers(std::vector<ram::Overlap>::const_iterator begin,
    std::vector<ram::Overlap>::const_iterator end) {

    if (begin >= end) {
        return;
    }

    std::vector<std::uint32_t> boundaries;
    for (auto it = begin; it != end; ++it) {
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
