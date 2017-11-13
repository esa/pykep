/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef KEP_TOOLBOX_EXCEPTIONS_H
#define KEP_TOOLBOX_EXCEPTIONS_H

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

#endif // KEP_TOOLBOX_EXCEPTIONS_H
