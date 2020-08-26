#pragma once
#ifdef PROFILE
#include <callgrind.h>
#endif

struct Profiler {
#ifdef PROFILE
    Profiler() {
        CALLGRIND_START_INSTRUMENTATION;
    }
    ~Profiler() {
        CALLGRIND_STOP_INSTRUMENTATION;
    }
#endif
};
