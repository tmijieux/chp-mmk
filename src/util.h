#ifndef TDP_UTIL_H
#define TDP_UTIL_H

#define SQUARE(x) ((x)*(x))

#include <vector>
#include <cassert>
#include "error.h"

#define ASSERT_MSG(msg, cond) assert(  ((void)(msg), (cond)) )

namespace chp {

struct vec : public std::vector<double>  {
    operator double* ()  { return &(*this)[0]; }
};

};

#endif // TDP_UTIL_H
