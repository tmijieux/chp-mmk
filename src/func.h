#ifndef CHP_FUNC_H
#define CHP_FUNC_H

#include <stdlib.h>
#include "proc.h"

typedef enum chp_func_type_ chp_func_type;

enum chp_func_type_ {
    CHP_STATIONARY,
    CHP_UNSTATIONARY,
};

typedef struct chp_func_ chp_func;

struct chp_func_ {
    chp_func_type type;
    const char *name;

    void (*rhs)(int const Nx,int const Ny,double const*X, double const*Y,double*F);
    void (*rhs_unsta)(int const Nx,int const Ny,
                      double const*X,double const*Y, double *F,
                      double const Lx, double const Ly, double const t);

    void (*bottom)(int const N,double const*A,double*B,double value);
    void (*top)(int const N,double const*A,double*B,double value);
    void (*right)(int const N,double const*A,double*B,double value);
    void (*left)(int const N,double const*A,double*B,double value);

    void (*U)(int const Nx, int const Ny, double const *X, double const *Y,
              double const Lx, double const Ly, double *U);

    chp_func *next;
};

extern unsigned func_list_length;
extern chp_func *func_list;

#define PASTE2_(x,y) x##y
#define PASTE2(x,y) PASTE2_(x, y)

#define REGISTER_FUNCTION(func_name, type_, bottom_, top_,      \
                          right_, left_, rhs_, rhs_unsta_, U_)  \
    static void __attribute__((constructor))                    \
    PASTE2(_register_func_,__COUNTER__)(void)                   \
    {                                                           \
        static chp_func meth = {                         \
            .type = CHP_##type_,                                \
            .name = #func_name   "  ("#type_")",                \
            .bottom = bottom_,                                  \
            .top = top_,                                        \
            .right = right_,                                    \
            .left = left_,                                      \
            .rhs = rhs_,                                        \
            .rhs_unsta = rhs_unsta_,                            \
            .U = U_,                                            \
            .next = NULL,                                       \
        };                                                      \
        meth.next = func_list;                                  \
        func_list = &meth;                                      \
        ++ func_list_length;                                    \
    }

void chp_func_specialize_rank(chp_func *func, int rank, int group_size);
void chp_func_init(chp_func *func, unsigned idx, chp_proc *P);
void chp_print_func_list(void);

#endif // CHP_FUNC_H
