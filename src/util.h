#ifndef TDP_UTIL_H
#define TDP_UTIL_H

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <immintrin.h> //AVX
#include <assert.h>
#include <omp.h>

#include "error.h"

#ifdef DEQUAL
#undef DEQUAL
#endif
#define DEQUAL(X_, Y_, H_) (fabs((X_) - (Y_)) < (H_))

#ifdef __FMA__
#define MM256_FMADD_PD(a, b, c) _mm256_fmadd_pd(a, b, c)
#else
#define MM256_FMADD_PD(a, b, c) _mm256_add_pd(_mm256_mul_pd(a, b), c)
#endif

#define SQUARE(x)  ((x)*(x))
#define SQUARE_UL(X) SQUARE(((uint64_t) (X)))

#define MALLOC_ARRAY(var, size) ((var) = malloc((size)*sizeof(*(var))))
#define CALLOC_ARRAY(var, size) ((var) = calloc((size), sizeof(*(var))))
#define ASSERT_MSG(msg, cond) assert(  ((void)(msg), (cond)) )


double *tdp_matrix_new(int m/*rows*/, int n/*columns*/);
double *tdp_avx256_aligned_matrix_new(int m/*rows*/, int n/*columns*/);
void tdp_matrix_zero(int m/*rows*/, int n/*columns*/, double *mat);
void tdp_matrix_one(int m/*row*/, int n/*column*/,
                    double value, double *mat, int lda/*leading dimension*/);
void tdp_matrix_fill(int m/*row*/, int n/*column*/,
                     double value, double *mat, int lda/*leading dimension*/);

void tdp_matrix_print(int m/*row*/, int n/*column*/,
                      double *mat, int lda/*leading dimension*/,
                      FILE *outstream);
void tdp_matrix_rand(int m/*rows*/, int n/*columns*/,
                     double *mat, double min, double max);

double *tdp_vector_new(int m);
void tdp_vector_rand(int m, double min, double max, double *v);
void tdp_vector_one(int m, double value, double *v);
void tdp_vector_zero(int m, double *v);
void tdp_vector_print(int m, double *v, FILE *out);

void tdp_print_cache_size(void);
uint64_t tdp_get_cache_size(int id);
double *tdp_cache_garbage(void);

#endif // TDP_UTIL_H
