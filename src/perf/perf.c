#include "perf.h"
#include <stdlib.h>
#include <stdio.h>

void perf(perf_t * p)
{
    gettimeofday(p, NULL);
}

void perf_diff(const perf_t * begin, perf_t * end)
{
    end->tv_sec = end->tv_sec - begin->tv_sec;
    end->tv_usec = end->tv_usec - begin->tv_usec;
    if (end->tv_usec < 0) {
        (end->tv_sec)--;
        end->tv_usec += 1000000;
    }
}

void perf_printh(const perf_t * p)
{
    long m = p->tv_sec / 60;
    long s = p->tv_sec - m*60;
    long ms = p->tv_usec / 1000;
    long micros = p->tv_usec - ms*1000;

    //  printf("%ld sec %ld usec\n", p->tv_sec, p->tv_usec);
    printf("%ld:%ld:%ld:%ld\n",m,s,ms,micros);
}

void perf_printmicro(const perf_t * p)
{
    printf("%ld\n",p->tv_usec + ( p->tv_sec * 1000000) );
}

uint64_t perf_get_micro(const perf_t *p)
{
    return (uint64_t) p->tv_usec + (p->tv_sec * 1000000UL);
}

double perf_mflops(const perf_t * p, const uint64_t nb_op)
{
    return (double)nb_op / (p->tv_sec * 1000000 + p->tv_usec);
}

