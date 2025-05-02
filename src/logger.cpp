#include <mtetcol/logger.h>

#include <spdlog/sinks/stdout_color_sinks.h>

namespace mtetcol {

spdlog::logger& logger()
{
    static auto default_logger = spdlog::stdout_color_mt("mtetcol");
    return *default_logger;
}

} // namespace mtetcol
