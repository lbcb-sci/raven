#ifndef RAVEN_IO_H_
#define RAVEN_IO_H_

#include <filesystem>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "raven/export.h"
#include "thread_pool/thread_pool.hpp"

namespace raven {

RAVEN_EXPORT
auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

RAVEN_EXPORT
auto LoadSequences(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                   std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;


}  // namespace raven

#endif  // RAVEN_IO_H_
