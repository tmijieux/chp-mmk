#ifndef CHP_TIMER_H
#define CHP_TIMER_H

#include "perf/perf.h"
#include "proc.hpp"

namespace chp {

class timer {
    perf_t _p1, _p2;

public:
    void start();
    void stop();
    void print(proc const& p) const;
    
    virtual ~timer(){}
};

};

#endif // CHP_TIMER_H
