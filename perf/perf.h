#ifndef PERF_H
#define PERF_H

#include <sys/time.h>
#include <stdint.h>

typedef struct timeval perf_t;

void
perf(perf_t * p);

void
perf_diff(const perf_t * begin, perf_t * end);

void
perf_printh(const perf_t * p);

void
perf_printmicro(const perf_t * p);

double
perf_mflops(const perf_t * p, const uint64_t nb_op);


uint64_t
perf_get_micro(const perf_t *p);


#endif
