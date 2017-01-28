#ifndef CHP_TIMER_H
#define CHP_TIMER_H

typedef struct chp_timer_ chp_timer;

#include "perf/perf.h"
#include "proc.h"

struct chp_timer_ {
    perf_t p1, p2;
};

void chp_timer_start(chp_timer *T);
void chp_timer_stop(chp_timer *T);
void chp_timer_print(chp_timer *T, chp_proc *P);

#endif // CHP_TIMER_H
