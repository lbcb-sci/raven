#ifndef RAVEN_IO_H_
#define RAVEN_IO_H_

#include <filesystem>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "raven/export.h"
#include "thread_pool/thread_pool.hpp"

namespace raven {

RAVEN_EXPORT
std::vector<std::unique_ptr<biosoup::NucleicAcid>> LoadSequences(
    const std::filesystem::path& path);

RAVEN_EXPORT
std::vector<std::unique_ptr<biosoup::NucleicAcid>> LoadSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    const std::vector<std::filesystem::path>& paths);

}  // namespace raven

#endif  // RAVEN_IO_H_
