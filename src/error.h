#ifndef PROJCHP_ERROR_H
#define PROJCHP_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>


#define PROJCHP_ERROR_LIST(ERROR)                  \
    ERROR(SUCCESS)                              \


#define PROJCHP_ERROR_TO_ENUM(error) PROJCHP_ERROR_##error,

enum projchp_error_code {
    PROJCHP_ERROR_LIST(PROJCHP_ERROR_TO_ENUM)
};

char const *projchp_errmsg(int errcode);

#ifndef __GNUC__
#define __PRETTY_FUNCTION__    __FUNCDNAME__
#endif


#define __FILENAME__ (strrchr(__FILE__, '/') ?                  \
                      strrchr(__FILE__, '/') + 1 : __FILE__)

#define projchp_error(format_, ...)                                \
    do {                                                        \
        fprintf(stderr, "ERROR: %s:%d|%s: ",  __FILENAME__ ,    \
                __LINE__, __PRETTY_FUNCTION__);                 \
        fprintf(stderr, (format_), ##__VA_ARGS__);              \
    } while(0)

#define projchp_fatal(format_, ...)                                        \
    do {                                                                \
        fprintf(stderr, "FATAL ERROR: %s:%d|%s: ",  __FILENAME__ ,      \
                __LINE__, __PRETTY_FUNCTION__);                         \
        fprintf(stderr, (format_), ##__VA_ARGS__);                      \
        exit(EXIT_FAILURE);                                             \
    } while(0)

#define projchp_warning(format_, ...)                              \
    do {                                                        \
        fprintf(stderr, "WARNING: %s:%d|%s: ",  __FILENAME__ ,  \
                __LINE__, __PRETTY_FUNCTION__);                 \
        fprintf(stderr, (format_), ##__VA_ARGS__);              \
    } while(0)

#ifdef DEBUG
#define projchp_debug(format_, ...)                                        \
    do {                                                                \
        fprintf(stderr, "DEBUG: %s:%d|%s: " format_, __FILENAME__ ,     \
                __LINE__, __PRETTY_FUNCTION__, ##__VA_ARGS__);          \
    } while(0)

#else // DEBUG
#define projchp_debug(format_, ...) ((void) (format_))
#endif // DEBUG

#endif // PROJCHP_ERROR_H
