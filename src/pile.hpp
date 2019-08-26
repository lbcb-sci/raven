/*!
 * @file pile.hpp
 *
 * @brief Pile class header file
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <memory>

namespace ram {
    struct Overlap;
}

namespace raven {

class Pile;
std::unique_ptr<Pile> createPile(std::uint32_t id, std::uint32_t size);

class Pile {
public:
    ~Pile() = default;

    void add_layers(const std::vector<ram::Overlap>& overlaps);

    std::string to_json() const;

    friend std::unique_ptr<Pile> createPile(std::uint32_t id, std::uint32_t size);
private:
    Pile(std::uint32_t id, std::uint32_t size);
    Pile(const Pile&) = delete;
    const Pile& operator=(const Pile&) = delete;

    std::uint32_t id_;
    std::vector<std::uint32_t> data_;
    std::uint32_t begin_;
    std::uint32_t end_;
    std::uint32_t median_;
};

}
