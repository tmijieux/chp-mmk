#include <stdio.h>
#include <mpi.h>
#include "timer.h"

void chp_timer_start(chp_timer *T)
{
    perf(&T->p1);
}

void chp_timer_stop(chp_timer *T)
{
    perf(&T->p1);
    perf_diff(&T->p1, &T->p2);
}

void chp_timer_print(chp_timer *T, chp_proc *P)
{
    uint64_t micro, max_t;

    micro = perf_get_micro(&T->p2);
    MPI_Reduce(&micro, &max_t, 1, MPI_UNSIGNED_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (!P->rank)
        fprintf(stderr, "# %lu.%06lu s\n", max_t/1000000UL, max_t%1000000UL);
}
