#include "error.h"

#define PROJCHP_ERROR_TO_STRING(error) [PROJCHP_ERROR_##error] = #error,

static char const *projchp_errmsg_[] = {
    PROJCHP_ERROR_LIST(PROJCHP_ERROR_TO_STRING)
};

char const *projchp_errmsg(int error_code)
{
    return projchp_errmsg_[error_code > 0 ? error_code : -error_code];
}
