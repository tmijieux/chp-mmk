#include <cstdio>
#include <mpi.h>

#include "timer.hpp"

using namespace chp;

void timer::start()
{
    perf(&_p1);
}

void timer::stop()
{
    perf(&_p2);
    perf_diff(&_p1, &_p2);
}

void timer::print(proc const& P) const
{
    uint64_t micro, max_t;

    micro = perf_get_micro(&_p2);
    MPI_Reduce(&micro, &max_t, 1, MPI_UNSIGNED_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (!P.rank())
        printf("# %lu.%06lu s\n", max_t/1000000UL, max_t%1000000UL);
}
