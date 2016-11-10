#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <malloc.h>

#include "util.h"
#include "cblas.h"

/**
 * Pretty print the matrix 'mat'
 */
void tdp_matrix_print(int m/*rows*/, int n/*columns*/,
                      double *mat, int lda/*leading dimension*/,
                      FILE *out)
{

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            fprintf(out, "%3g ", mat[j*lda+i]);
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

/**
 * Return a new zero'd (m x n) matrix
 */
double *tdp_matrix_new(int m/*rows*/, int n/*columns*/)
{
    double *d;
    d = calloc(m*n, sizeof*d);
    return d;
}

/**
 * Set the matrix elements to random values taken in interval [min, max]
 */
void tdp_matrix_rand(int m/*rows*/, int n/*columns*/,
                     double *mat, double min, double max)
{
    for (int i = 0 ; i < m*n; ++i)
        mat[i] = min + ((double)rand() / RAND_MAX)*(max-min);
}

/**
 * Set matrix elements to 0.0
 */
void tdp_matrix_zero(int m/*rows*/, int n/*columns*/, double *mat)
{
    memset(mat, 0, sizeof*mat*m*n);
}

/**
 * Set the matrix' main diagonal elements to 'value', and other to 0.0
 */
void tdp_matrix_one(int m/*rows*/, int n/*columns*/,
                    double value, double *mat, int lda/*leading dimension*/)
{
    tdp_matrix_zero(m, n, mat);
    int M = min(m, n);
    for (int j = 0; j < M; ++j)
        mat[j*lda+j] = value;
}


/**
 * Set the matrix' main elements to 'value'
 */
void tdp_matrix_fill(int m/*rows*/, int n/*columns*/,
                     double value, double *mat, int lda/*leading dimension*/)
{
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            mat[j*lda+i] = value;
}

/**
 * Set the matrix' main diagonal elements to 'value', and other to 0.0
 */
void tdp_matrix_3one(int m/*rows*/, int n/*columns*/,
                     double v1, double v2,
                     double *mat, int lda/*leading dimension*/)
{
    tdp_matrix_zero(m, n, mat);
    int M = min(m, n);

    mat[0] = v1;
    mat[1] = v2;
    for (int j = 1; j < M-1; ++j) {
        mat[j*lda+j-1] = v2;
        mat[j*lda+j] = v1;
        mat[j*lda+j+1] = v2;
    }
    mat[(M-1)*lda+(M-2)] = v2;
    mat[(M-1)*lda+(M-1)] = v1;
}

/**
 * Return new zero'd vector
 */
double *tdp_vector_new(int m)
{
    double *d;
    d = calloc(m, sizeof*d);
    return d;
}

/**
 * Set vector elements with random values taken in interval [min, max]
 */
void tdp_vector_rand(int m, double min, double max, double *v)
{
    for (int i = 0; i < m; ++i)
        v[i] = min + ((double)rand() / RAND_MAX) * (max-min);
}

/**
 * Set all vector elements to 'value'
 */
void tdp_vector_one(int m, double value, double *v)
{
    for (int i = 0; i < m; ++i)
        v[i] = value;
}

/**
 * Set vector elements to 0.0
 */
void tdp_vector_zero(int m, double *v)
{
    memset(v, 0, m*sizeof v[0]);
}

/**
 * Pretty print the vector
 */
void tdp_vector_print(int m, double *v, FILE *out)
{
    for (int i = 0; i < m; ++i)
        fprintf(out, "%g\n", v[i]);
}

double *tdp_cache_garbage(void)
{
    uint64_t S = tdp_get_cache_size(3)*10;
    double *a = memalign(64, S);

    uint64_t s = S * (((uint64_t)log(S))+1);
    while (s > 0) {
        int i = rand() % (S/sizeof *a);
        int k = rand() % (S/sizeof *a);
        a[i] = a[k];
        s -= sizeof *a;
    }
    return a;
}

#define b(val, base, end)                                       \
    ((val << (__WORDSIZE-end-1)) >> (__WORDSIZE-end+base-1))

#define GET_CACHE_DETAILS(level_, ways_, partitions_, line_size_, sets_) \
    do {                                                                \
        uint64_t eax_, ebx_, ecx_, edx_;                                \
        __asm__( "cpuid"                                                \
                 : "=a"(eax_), "=b"(ebx_), "=c"(ecx_), "=d"(edx_)       \
                 : "a"(4), "b"(0), "c"(level_), "d"(0));                \
        ways_ = b(ebx_, 22, 31) + 1;                                    \
        partitions_ = b(ebx_, 12, 21) + 1;                              \
        line_size_ = b(ebx_, 0, 11) + 1;                                \
        sets_ = ecx_ + 1;                                               \
    }while(0)

#define GETCACHESIZE(level_)                                            \
    ({                                                                  \
        uint64_t ways, partitions, line_size, sets;                     \
        GET_CACHE_DETAILS(level_, ways, partitions, line_size, sets);   \
        (ways * partitions * line_size * sets) / 1024;                  \
    })

#define PRINT_CACHE_DETAILS(level_)                                     \
    do {                                                                \
        uint64_t ways, partitions, line_size, sets;                     \
        GET_CACHE_DETAILS(level_, ways, partitions, line_size, sets);   \
        printf("line size: %ld\n", line_size);                          \
        printf("ways: %ld\n", ways);                                    \
        printf("sets: %ld\n", sets);                                    \
        printf("partitions: %ld\n\n", partitions);                      \
    }while(0)

void tdp_print_cache_size(void)
{
    printf("L1 cache_size %ldK\n", GETCACHESIZE(1));
    PRINT_CACHE_DETAILS(1);

    printf("L2 cache_size %ldK\n", GETCACHESIZE(2));
    PRINT_CACHE_DETAILS(2);

    printf("L3 cache_size %ldK\n", GETCACHESIZE(3));
    PRINT_CACHE_DETAILS(3);
}

uint64_t tdp_get_cache_size(int id)
{
    assert(id >= 1 && id <= 3);
    return GETCACHESIZE(id);
}
