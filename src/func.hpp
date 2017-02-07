#ifndef CHP_FUNC_H
#define CHP_FUNC_H

#include <cstdlib>
#include <vector>
#include <string>

#include "proc.hpp"


namespace chp {

class func {
private:
    static std::vector<func*> func_list;
    const func& specialize_rank(proc const &);

public:

    enum type {
        STATIONARY,
        UNSTATIONARY,
    };

    typedef void (*rhs_t)(int const,int const,double const*,double const*,double*);
    typedef void (*rhs_unsta_t)(int const,int const,
                                double const*,double const*,double*,
                                double const, double const,double const);
    typedef void (*border_t)(int const, double const*, double*, double);
    typedef void (*U_t)(int const, int const, double const*, double const*,
                        double const, double const, double*);

    type m_type;
    std::string m_name;

    rhs_t m_rhs;
    rhs_unsta_t m_rhs_unsta;
    border_t m_bottom;
    border_t m_top;
    border_t m_right;
    border_t m_left;
    U_t m_U;

    func(std::string const&, func::type,
         rhs_t, rhs_unsta_t, border_t, border_t, border_t, border_t, U_t);

    static func get_func(unsigned idx, proc const &P);
    static void print_list();
};

#define PASTE2_(x,y) x##y
#define PASTE2(x,y) PASTE2_(x, y)

#define REGISTER_FUNCTION(func_name, type_, bottom_, top_,      \
                          right_, left_, rhs_, rhs_unsta_, U_)  \
static func PASTE2(_meth_,__COUNTER__)(                         \
    (#func_name   "  ("#type_")"),                              \
    (func::type_),                                              \
    rhs_,                                                       \
    rhs_unsta_,                                                 \
    bottom_,                                                    \
    top_,                                                       \
    right_,                                                     \
    left_,                                                      \
    U_                                                          \
);                                                              \

}; //namespace chp

#endif // CHP_FUNC_H
