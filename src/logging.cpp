#include "centrolign/logging.hpp"

#include <mutex>
#include <iostream>
#include <ctime>
#include <sstream>
#include <atomic>
#include <cstdint>
#include <thread>
#include <map>
#include <vector>

namespace centrolign {

using namespace std;

static recursive_mutex monitor;
static atomic<uint64_t> next_thread_id{0};

uint64_t thread_id()
{
    thread_local uint64_t id = next_thread_id.fetch_add(1);
    return id;
}

// the most recent message from each running thread (allocated and never destroyed so that
// threads that exit during static destruction can still access it)
static map<uint64_t, string>& recent_messages() {
    static map<uint64_t, string>* messages = new map<uint64_t, string>();
    return *messages;
}

// clears a thread's recent message when the thread exits
struct RecentMessageGuard {
    RecentMessageGuard() = default;
    ~RecentMessageGuard() {
        lock_guard<recursive_mutex> lock(monitor);
        recent_messages().erase(thread_id());
    }
};

static void record_recent_message(const string& msg) {
    // skip memory reports, which don't indicate what the thread is doing
    static const vector<string> skip_prefixes{"Current memory", "Peak memory", "Allocator free memory"};
    for (const auto& prefix : skip_prefixes) {
        if (msg.compare(0, prefix.size(), prefix) == 0) {
            return;
        }
    }
    thread_local RecentMessageGuard guard;
    recent_messages()[thread_id()] = msg;
}

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

logging::LoggingLevel logging::level = logging::Basic;

void logging::log(logging::LoggingLevel priority, const std::string& msg) {
    if (priority <= level) {
        lock_guard<recursive_mutex> lock(monitor);
        
        time_t time_now;
        time(&time_now);
        clock_t clock_now = clock();
        
        tm* local = localtime(&time_now);
        static char buffer[100];
        strftime(buffer, 100, "%Y-%m-%d %H:%M", local);
        
        double wall_secs = difftime(time_now, start.start_time);
        double cpu_secs = (clock_now - start.start_clock) / CLOCKS_PER_SEC;
        
        cerr << "[" << buffer << " (thread " << thread_id() << "), elapsed: " << format_seconds(wall_secs) << " wall / " << format_seconds(cpu_secs) << " cpu] " << msg;
        if (msg.back() != '\n') {
            cerr << '\n';
        }
        
        record_recent_message(msg);
    }
}

vector<pair<uint64_t, string>> logging::recent_activity() {
    lock_guard<recursive_mutex> lock(monitor);
    vector<pair<uint64_t, string>> activity;
    uint64_t self = thread_id();
    for (const auto& rec : recent_messages()) {
        if (rec.first != self) {
            activity.emplace_back(rec.first, rec.second);
        }
    }
    return activity;
}

}
