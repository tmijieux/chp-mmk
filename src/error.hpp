#ifndef CHP_ERROR_H
#define CHP_ERROR_H

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <errno.h>

#include <exception>
#include <string>
#include <boost/format.hpp>
#include <sstream>

namespace chp {

class exception : public std::exception {
private:
    boost::format _f;

public:
    template<typename  T>
    exception &operator%(T a)
    {
        _f = _f % a;
        return *this;
    }

    exception(const std::string& s): _f(s)  { }
    const char *what()
    {
        std::stringstream ss;
        ss << _f;
        return ss.str().c_str();
    }
};

};

#define CHP_ERROR_LIST(ERROR)                   \
    ERROR(SUCCESS)                              \


#define CHP_ERROR_TO_ENUM(error) CHP_ERROR_##error,

enum chp_error_code {
    CHP_ERROR_LIST(CHP_ERROR_TO_ENUM)
};

char const *chp_errmsg(int errcode);

#ifndef __GNUC__
#define __PRETTY_FUNCTION__    __FUNCDNAME__
#endif


#define __FILENAME__ (strrchr(__FILE__, '/') ?                  \
                      strrchr(__FILE__, '/') + 1 : __FILE__)

#define chp_error(format_, ...)                                 \
    do {                                                        \
        fprintf(stderr, "ERROR: %s:%d|%s: ",  __FILENAME__ ,    \
                __LINE__, __PRETTY_FUNCTION__);                 \
        fprintf(stderr, (format_), ##__VA_ARGS__);              \
    } while(0)

#define chp_fatal(format_, ...)                                         \
    do {                                                                \
        fprintf(stderr, "FATAL ERROR: %s:%d|%s: ",  __FILENAME__ ,      \
                __LINE__, __PRETTY_FUNCTION__);                         \
        fprintf(stderr, (format_), ##__VA_ARGS__);                      \
        exit(EXIT_FAILURE);                                             \
    } while(0)

#define chp_warning(format_, ...)                               \
    do {                                                        \
        fprintf(stderr, "WARNING: %s:%d|%s: ",  __FILENAME__ ,  \
                __LINE__, __PRETTY_FUNCTION__);                 \
        fprintf(stderr, (format_), ##__VA_ARGS__);              \
    } while(0)

#ifdef DEBUG
#define chp_debug(format_, ...)                                         \
    do {                                                                \
        fprintf(stderr, "DEBUG: %s:%d|%s: " format_, __FILENAME__ ,     \
                __LINE__, __PRETTY_FUNCTION__, ##__VA_ARGS__);          \
    } while(0)

#else // DEBUG
#define chp_debug(format_, ...) ((void) (format_))
#endif // DEBUG

#endif // CHP_ERROR_H
