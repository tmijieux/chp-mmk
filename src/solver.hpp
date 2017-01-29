#ifndef CHP_SOLVER_H
#define CHP_SOLVER_H

#include <vector>
#include <string>
#include "equation.hpp"

namespace chp {

class solver {

protected:
    static constexpr const double EPSILON  = 1e-8;
    static const int MAX_NB_ITER = 20000;

    int Nx, Ny, N;
    double B, Cx, Cy;
    std::vector<double> _tmp[4];

    void gemv_5diag(const double *__restrict__ X, double *__restrict__ AX);
    virtual int _solve(const double *__restrict__ rhs, double *__restrict__ X) = 0;
    solver(equation const& eq);

public:
    static solver *make_solver(equation const& eq, std::string const &solver_arg);
    int solve(const double *rhs, const double *X0, double *X);
    virtual ~solver() {}
};


class gauss_seidel_solver : public solver {
private:
    void impl(const double * __restrict__ rhs, double *__restrict__ X_new,
              const double *__restrict__ X_old);
    virtual int _solve(const double *__restrict__ rhs, double *__restrict__ X) override;

public:
    gauss_seidel_solver(equation const& eq) : solver(eq) {}
};


class jacobi_solver : public solver {
private:
    void impl(const double * __restrict__ rhs, double *__restrict__ X_new,
              const double *__restrict__ X_old);
    virtual int _solve(const double *__restrict__ rhs, double *__restrict__ X) override;

public:
    jacobi_solver(equation const& eq) : solver(eq) {}
};

class CG_solver : public solver {
private:
    void impl(const double * __restrict__ rhs, double *__restrict__ X_new,
              const double *__restrict__ X_old);
    virtual int _solve(const double *__restrict__ rhs, double *__restrict__ X) override;

public:
    CG_solver(equation const& eq) : solver(eq) {}
};



};

#endif // CHP_SOLVER_H
