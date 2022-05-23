#include "raven/io.h"

#include <array>
#include <fstream>
#include <functional>
#include <limits>
#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

namespace raven {

namespace detail {

static constexpr auto kFastaSuffxies =
    std::array<char const*, 4>{".fasta", "fasta.gz", ".fa", ".fa.gz"};

static constexpr auto kFastqSuffixes =
    std::array<char const*, 4>{".fastq", ".fastq.gz", ".fq", ".fa.gz"};

static auto IsSuffixFor(std::string_view const suffix,
                        std::string_view const query) -> bool {
  return suffix.length() <= query.length()
             ? suffix == query.substr(query.length() - suffix.length())
             : false;
}

static auto CreateParser(std::filesystem::path const& path)
    -> std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> {
  using namespace std::placeholders;
  if (std::filesystem::exists(path)) {
    if (std::any_of(kFastaSuffxies.cbegin(), kFastaSuffxies.cend(),
                    std::bind(IsSuffixFor, _1, path.c_str()))) {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastaParser>(path.c_str());
    } else if (std::any_of(kFastqSuffixes.cbegin(), kFastqSuffixes.cend(),
                           std::bind(IsSuffixFor, _1, path.c_str()))) {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastqParser>(path.c_str());
    }
  }

  throw std::invalid_argument(
      "[raven::detail::CreateParser] invalid file path: " + path.string());
}

}  // namespace detail

auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto sparser = detail::CreateParser(path);

  return sparser->Parse(std::numeric_limits<std::uint64_t>::max());
}

auto LoadSequences(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                   std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  {
    auto parse_futures = std::vector<
        std::future<std::vector<std::unique_ptr<biosoup::NucleicAcid>>>>();

    auto parse_async = [&thread_pool](std::filesystem::path const& path)
        -> std::future<std::vector<std::unique_ptr<biosoup::NucleicAcid>>> {
      return thread_pool->Submit(
          [](std::filesystem::path const& path)
              -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
            return LoadSequences(path);
          },
          std::cref(path));
    };

    std::transform(paths.cbegin(), paths.cend(),
                   std::back_inserter(parse_futures), parse_async);

    for (auto& pf : parse_futures) {
      auto sequences = pf.get();
      dst.insert(dst.end(), std::make_move_iterator(sequences.begin()),
                 std::make_move_iterator(sequences.end()));
    }

    std::sort(dst.begin(), dst.end(),
              [](std::unique_ptr<biosoup::NucleicAcid> const& lhs,
                 std::unique_ptr<biosoup::NucleicAcid> const& rhs) -> bool {
                return lhs->id < rhs->id;
              });
  }

  decltype(dst)(std::make_move_iterator(dst.begin()),
                std::make_move_iterator(dst.end()))
      .swap(dst);

  return dst;
}

}  // namespace raven
