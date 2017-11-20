/***************************************************************************
 *   Copyright (C) 2009 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef BOOST_PYTHON_P_EXCEPTIONS_H
#define BOOST_PYTHON_P_EXCEPTIONS_H

#include <boost/cast.hpp>
#include <boost/python/exception_translator.hpp>
#include <exception>

inline void ie_translator(const index_error &ie)
{
    PyErr_SetString(PyExc_IndexError, ie.what());
}

inline void ve_translator(const value_error &ve)
{
    PyErr_SetString(PyExc_ValueError, ve.what());
}

inline void te_translator(const type_error &te)
{
    PyErr_SetString(PyExc_TypeError, te.what());
}

inline void ae_translator(const assertion_error &ae)
{
    PyErr_SetString(PyExc_AssertionError, ae.what());
}

inline void nie_translator(const not_implemented_error &nie)
{
    PyErr_SetString(PyExc_NotImplementedError, nie.what());
}

inline void me_translator(const memory_error &me)
{
    PyErr_SetString(PyExc_MemoryError, me.what());
}

inline void zde_translator(const zero_division_error &zde)
{
    PyErr_SetString(PyExc_ZeroDivisionError, zde.what());
}

// Translators for some standard and Boost exceptions.

inline void ba_translator(const std::bad_alloc &ba)
{
    PyErr_SetString(PyExc_MemoryError, ba.what());
}

inline void oe_translator(const std::overflow_error &oe)
{
    PyErr_SetString(PyExc_OverflowError, oe.what());
}

inline void bnc_translator(const boost::numeric::bad_numeric_cast &bnc)
{
    PyErr_SetString(PyExc_OverflowError, bnc.what());
}

// Translate our C++ exceptions into Python exceptions.
inline void translate_p_exceptions()
{
    boost::python::register_exception_translator<index_error>(ie_translator);
    boost::python::register_exception_translator<value_error>(ve_translator);
    boost::python::register_exception_translator<type_error>(te_translator);
    boost::python::register_exception_translator<assertion_error>(ae_translator);
    boost::python::register_exception_translator<not_implemented_error>(nie_translator);
    boost::python::register_exception_translator<memory_error>(me_translator);
    boost::python::register_exception_translator<zero_division_error>(zde_translator);
    boost::python::register_exception_translator<std::bad_alloc>(ba_translator);
    boost::python::register_exception_translator<std::overflow_error>(oe_translator);
    boost::python::register_exception_translator<boost::numeric::bad_numeric_cast>(bnc_translator);
}

#endif
