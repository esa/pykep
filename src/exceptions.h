#ifndef KEPLERIAN_TOOLBOX_EXCEPTIONS_H
#define KEPLERIAN_TOOLBOX_EXCEPTIONS_H

#include <string>

#ifdef USE_PAGMO
inline void throw_value_error(std::string s) {
    pagmo_throw(value_error, s);
}
#else

class kep_toolbox_error : public std::exception {
public:
    kep_toolbox_error(std::string _message)
	: message(_message)
    {}

    virtual ~kep_toolbox_error() throw()
    {}
    
    virtual const char* what() const throw() {
	return message.c_str();
    }
private:
    std::string message;
};

inline void throw_value_error(std::string s) {
    throw kep_toolbox_error(s);
}
#endif // USE_PAGMO

#endif // KEPLERIAN_TOOLBOX_EXCEPTIONS_H
