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

#ifndef KEP_TOOLBOX_IO_HPP
#define KEP_TOOLBOX_IO_HPP

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define KEP_TOOLBOX_MAX_OUTPUT_LENGTH 20u

namespace kep_toolbox
{

// Forward declaration
template <typename... Args>
inline void stream(std::ostream &, const Args &...);

namespace detail
{

template <typename T>
inline void stream_impl(std::ostream &os, const T &x)
{
    os << x;
}

inline void stream_impl(std::ostream &os, const bool &b)
{
    if (b) {
        os << "true";
    } else {
        os << "false";
    }
}

template <typename T>
inline void stream_impl(std::ostream &os, const std::vector<T> &v)
{
    auto len = v.size();
    if (len <= KEP_TOOLBOX_MAX_OUTPUT_LENGTH) {
        os << '[';
        for (decltype(v.size()) i = 0u; i < v.size(); ++i) {
            stream(os, v[i]);
            if (i != v.size() - 1u) {
                os << ", ";
            }
        }
        os << ']';
    } else {
        os << '[';
        for (decltype(v.size()) i = 0u; i < KEP_TOOLBOX_MAX_OUTPUT_LENGTH; ++i) {
            stream(os, v[i], ", ");
        }
        os << "... ]";
    }
}

template <typename T, typename U>
inline void stream_impl(std::ostream &os, const std::pair<T, U> &p)
{
    stream(os, '(', p.first, ',', p.second, ')');
}

template <typename T, typename U>
inline void stream_impl(std::ostream &os, const std::map<T, U> &m)
{
    unsigned counter = 0;
    stream(os, '{');
    for (auto it = m.begin(); it != m.end(); ++counter) {
        if (counter == KEP_TOOLBOX_MAX_OUTPUT_LENGTH) {
            stream(os, "...");
            break;
        }
        stream(os, it->first, " : ", it->second);
        ++it;
        if (it != m.end()) {
            stream(os, ",  ");
        }
    }
    stream(os, '}');
}

template <typename T, typename... Args>
inline void stream_impl(std::ostream &os, const T &x, const Args &...args)
{
    stream_impl(os, x);
    stream_impl(os, args...);
}

// A small helper function that transforms x to string, using internally kep_toolbox::stream.
template <typename T>
inline std::string to_string(const T &x)
{
    std::ostringstream oss;
    stream(oss, x);
    return oss.str();
}

} // end of namespace detail

/// The kep_toolbox streaming function.
/**
 * This function will direct to the output stream \p os the input arguments \p args.
 *
 * @param os the target stream.
 * @param args the objects that will be directed to to \p os.
 */
template <typename... Args>
inline void stream(std::ostream &os, const Args &...args)
{
    detail::stream_impl(os, args...);
}

/// The kep_toolbox print function.
/**
 * This function is equivalent to calling kep_toolbox::stream with \p std::cout as first argument.
 *
 * @param args the objects that will be printed to screen.
 */
template <typename... Args>
inline void print(const Args &...args)
{
    stream(std::cout, args...);
}

} // end of namespace kep_toolbox

#undef KEP_TOOLBOX_MAX_OUTPUT_LENGTH

#endif
