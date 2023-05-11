#include "centrolign/logging.hpp"

#include <mutex>
#include <iostream>
#include <ctime>
#include <sstream>

namespace centrolign {

using namespace std;

namespace logging {

static recursive_mutex monitor;

struct LoggingStart {
    
    time_t start_time;
    clock_t start_clock;
    
    LoggingStart() {
        time(&start_time);
        start_clock = clock();
    }
    ~LoggingStart() = default;
};

static LoggingStart start;

string format_seconds(double secs) {
    stringstream strm;
    strm << fixed;
    strm.precision(0);
    if (secs <= 60.0) {
        strm << secs << " s";
    }
    else {
        strm.precision(1);
        double mins = secs / 60.0;
        if (mins <= 60.0) {
            strm << mins << " m";
        }
        else {
            double hrs = mins / 60.0;
            if (hrs <= 24.0) {
                strm << hrs << " h";
            }
            else {
                strm << (hrs / 24.0) << " d";
            }
        }
    }
    return strm.str();
}

void log(LoggingLevel priority, const std::string& msg) {
    if (priority <= level) {
        lock_guard<recursive_mutex> lock(monitor);
        
        time_t time_now;
        time(&time_now);
        clock_t clock_now = clock();
        
        double wall_secs = difftime(time_now, start.start_time);
        double cpu_secs = (clock_now - start.start_clock) / CLOCKS_PER_SEC;
        
        cerr << "[elapsed time: " << format_seconds(wall_secs) << " wall / " << format_seconds(cpu_secs) << " cpu] " << msg;
        if (msg.back() != '\n') {
            cerr << '\n';
        }
        
    }
}

}

}
