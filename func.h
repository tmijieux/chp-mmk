#ifndef PROJCHP_FUNC_H
#define PROJCHP_FUNC_H

#include <stdlib.h>

enum projchp_method_type {
    PROJCHP_STATIONARY,
    PROJCHP_UNSTATIONARY,
};

struct projchp_method {
    enum projchp_method_type type;
    const char *name;
    void (*g)(int const N, double const *A, double *B, double param);
    void (*h)(int const N, double const *A, double *B,
              double Lx_min, double Lx_max);
    
    void (*f)(int const Nx, int const Ny,
              double const *X, double const *Y,
              double *F);
    void (*fu)(int const Nx, int const Ny,
               double const *X, double const *Y,
               double *F, double const Lx_min, double Lx_max,
               double const Ly, double const t);
    void (*mU)(int const Nx, int const Ny,
               double const *X, double const *Y,
               double const Lx_min, double const Lx_max,
               double const Ly, double *U);
    struct projchp_method *next;
};

extern unsigned method_list_length;
extern struct projchp_method *method_list;

#define PASTE2_(x,y) x##y
#define PASTE2(x,y) PASTE2_(x, y)

#define REGISTER_FUNCTION(method_name, type_, g_, h_, f_, fu_, mU_)     \
    static void __attribute__((constructor))                            \
    PASTE2(_register_method_,__COUNTER__)(void)                         \
    {                                                                   \
        static struct projchp_method meth = {                           \
            .type = PROJCHP_##type_,                                    \
            .name = #method_name   "  ("#type_")",                      \
            .g = g_,                                                    \
            .h = h_,                                                    \
            .f = f_,                                                    \
            .fu = fu_,                                                  \
            .mU = mU_,                                                  \
            .next = NULL,                                               \
        };                                                              \
        meth.next = method_list;                                        \
        method_list = &meth;                                            \
        ++ method_list_length;                                          \
    }

#endif // PROJCHP_FUNC_H
