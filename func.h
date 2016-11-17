#ifndef CHP_FUNC_H
#define CHP_FUNC_H

#include <stdlib.h>

enum chp_func_type {
    CHP_STATIONARY,
    CHP_UNSTATIONARY,
};

struct chp_func
{
    enum chp_func_type type;
    const char *name;

    void (*rhs)(int const Nx,int const Ny,double const*X, double const*Y,double*F);
    void (*rhs_unsta)(int const Nx,int const Ny,double const*X,double const*Y,
                      double *F, double const Lx, double const Ly,
                      double const t);

    void (*bottom)(int const N,double const*A,double*B,double value);
    void (*top)(int const N,double const*A,double*B,double value);
    void (*right)(int const N,double const*A,double*B,double value);
    void (*left)(int const N,double const*A,double*B,double value);

    void (*U)(int const Nx, int const Ny, double const *X, double const *Y,
              double const Lx, double const Ly, double *U);

    struct chp_func *next;
};

extern unsigned func_list_length;
extern struct chp_func *func_list;

#define PASTE2_(x,y) x##y
#define PASTE2(x,y) PASTE2_(x, y)

#define REGISTER_FUNCTION(func_name, type_, bottom_, top_,      \
                          right_, left_, rhs_, rhs_unsta_, U_)  \
    static void __attribute__((constructor))                    \
    PASTE2(_register_func_,__COUNTER__)(void)                   \
    {                                                           \
        static struct chp_func meth = {                         \
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

void chp_func_specialize_rank(struct chp_func *func, int rank, int group_size);
struct chp_func *chp_get_func_by_name(char const *name);
struct chp_func *chp_get_func_by_id(unsigned idx);

#endif // CHP_FUNC_H
