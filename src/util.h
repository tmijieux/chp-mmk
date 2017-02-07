#ifndef CHP_UTIL_H
#define CHP_UTIL_H

#define SQUARE(x) ((x)*(x))

#include <vector>
#include <cassert>
#include "error.h"

#define ASSERT_MSG(msg, cond) assert(  ((void)(msg), (cond)) )

namespace chp {

struct vec : public std::vector<double>  {
    inline operator double* () { return &(*this)[0]; }
    inline operator const double* () const { return &(*this)[0]; }
};

};

#endif // CHP_UTIL_H
