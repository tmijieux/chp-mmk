#ifndef CHP_ERROR_H
#define CHP_ERROR_H

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

#endif // CHP_ERROR_H
