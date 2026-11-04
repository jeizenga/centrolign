#ifndef centrolign_logging_hpp
#define centrolign_logging_hpp

#include <string>
#include <mutex>
#include <vector>
#include <utility>
#include <cstdint>

namespace centrolign {

/*
 * Struct namespace for logging
 */
struct logging {
    
    // priority level for logging messages
    enum LoggingLevel {
        Silent = 0,
        Minimal = 1,
        Basic = 2,
        Verbose = 3,
        Debug = 4
    };

    // the level of logging we are currently using
    static LoggingLevel level;

    static void log(LoggingLevel priority, const std::string& msg);
    
    // the most recent message logged by each thread that is still running, as pairs of
    // (thread ID, message), excluding the calling thread
    static std::vector<std::pair<uint64_t, std::string>> recent_activity();
};

}

#endif /* centrolign_logging_hpp */
