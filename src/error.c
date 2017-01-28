#include "error.h"

#define CHP_ERROR_TO_STRING(error) [CHP_ERROR_##error] = #error,

static char const *chp_errmsg_[] = {
    CHP_ERROR_LIST(CHP_ERROR_TO_STRING)
};

char const *chp_errmsg(int error_code)
{
    return chp_errmsg_[error_code > 0 ? error_code : -error_code];
}
