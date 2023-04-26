#ifndef centrolign_logging_hpp
#define centrolign_logging_hpp

#include <string>

namespace centrolign {

namespace logging {

// priority level for logging messages
enum LoggingLevel {
    Silent = 0,
    Minimal = 1,
    Basic = 2,
    Verbose = 3,
    Debug = 4
};

// the level of logging we are currently using
static LoggingLevel level = Basic;

void log(LoggingLevel priority, const std::string& msg);

}

}

#endif /* centrolign_logging_hpp */
