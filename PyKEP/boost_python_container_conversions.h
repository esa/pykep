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

// Slightly adapted from:
// http://cctbx.svn.sourceforge.net/viewvc/cctbx/trunk/scitbx/boost_python/container_conversions.h

// Original license agreement follows.

/*

*** License agreement ***

cctbx Copyright (c) 2006, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).  All
rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes,
patches, or upgrades to the features, functionality or performance of
the source code ("Enhancements") to anyone; however, if you choose to
make your Enhancements available either publicly, or directly to
Lawrence Berkeley National Laboratory, without imposing a separate
written license agreement for such Enhancements, then you hereby grant
the following license: a  non-exclusive, royalty-free perpetual license
to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

*/

#ifndef BOOST_PYTHON_CONTAINER_CONVERSIONS_H
#define BOOST_PYTHON_CONTAINER_CONVERSIONS_H

#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/tuple.hpp>

template <typename ContainerType> struct to_tuple {
  static PyObject *convert(ContainerType const &a) {
    boost::python::list result;
    typedef typename ContainerType::const_iterator const_iter;
    for (const_iter p = a.begin(); p != a.end(); p++) {
      result.append(boost::python::object(*p));
    }
    return boost::python::incref(boost::python::tuple(result).ptr());
  }

  static const PyTypeObject *get_pytype() { return &PyTuple_Type; }
};

struct default_policy {
  static bool check_convertibility_per_element() { return false; }

  template <typename ContainerType>
  static bool check_size(boost::type<ContainerType>, std::size_t /*sz*/) {
    return true;
  }

  template <typename ContainerType>
  static void assert_size(boost::type<ContainerType>, std::size_t /*sz*/) {}

  template <typename ContainerType>
  static void reserve(ContainerType &a, std::size_t sz) {}
};

struct fixed_size_policy {
  static bool check_convertibility_per_element() { return true; }

  template <typename ContainerType>
  static bool check_size(boost::type<ContainerType>, std::size_t sz) {
    return std::tuple_size<ContainerType>::value == sz;
  }

  template <typename ContainerType>
  static void assert_size(boost::type<ContainerType>, std::size_t sz) {
    if (!check_size(boost::type<ContainerType>(), sz)) {
      PyErr_SetString(PyExc_RuntimeError,
                      "Insufficient elements for fixed-size array.");
      boost::python::throw_error_already_set();
    }
  }

  template <typename ContainerType>
  static void reserve(ContainerType & /*a*/, std::size_t sz) {
    if (sz > std::tuple_size<ContainerType>::value) {
      PyErr_SetString(PyExc_RuntimeError,
                      "Too many elements for fixed-size array.");
      boost::python::throw_error_already_set();
    }
  }

  template <typename ContainerType, typename ValueType>
  static void set_value(ContainerType &a, std::size_t i, ValueType const &v) {
    reserve(a, i + 1);
    a[i] = v;
  }
};

struct variable_capacity_policy : default_policy {
  template <typename ContainerType>
  static void reserve(ContainerType &a, std::size_t sz) {
    a.reserve(sz);
  }

  template <typename ContainerType, typename ValueType>
  static void set_value(ContainerType &a, std::size_t
#if !defined(NDEBUG)
                                              i
#endif
                        ,
                        ValueType const &v) {
    assert(a.size() == i);
    a.push_back(v);
  }
};

struct fixed_capacity_policy : variable_capacity_policy {
  template <typename ContainerType>
  static bool check_size(boost::type<ContainerType>, std::size_t sz) {
    return ContainerType::max_size() >= sz;
  }
};

struct linked_list_policy : default_policy {
  template <typename ContainerType, typename ValueType>
  static void set_value(ContainerType &a, std::size_t /*i*/,
                        ValueType const &v) {
    a.push_back(v);
  }
};

struct set_policy : default_policy {
  template <typename ContainerType, typename ValueType>
  static void set_value(ContainerType &a, std::size_t /*i*/,
                        ValueType const &v) {
    a.insert(v);
  }
};

template <typename ContainerType, typename ConversionPolicy>
struct from_python_sequence {
  typedef typename ContainerType::value_type container_element_type;

  from_python_sequence() {
    boost::python::converter::registry::push_back(
        &convertible, &construct, boost::python::type_id<ContainerType>());
  }

  static void *convertible(PyObject *obj_ptr) {
    if (!(PyList_Check(obj_ptr) || PyTuple_Check(obj_ptr) ||
          PyIter_Check(obj_ptr) || PyRange_Check(obj_ptr) ||
          (!PyBytes_Check(obj_ptr) && !PyUnicode_Check(obj_ptr) &&
           (obj_ptr->ob_type == 0 || Py_TYPE(obj_ptr->ob_type) == 0 ||
            Py_TYPE(obj_ptr->ob_type)->tp_name == 0 ||
            std::strcmp(Py_TYPE(obj_ptr->ob_type)->tp_name,
                        "Boost.Python.class") != 0) &&
           PyObject_HasAttrString(obj_ptr, "__len__") &&
           PyObject_HasAttrString(obj_ptr, "__getitem__"))))
      return 0;
    boost::python::handle<> obj_iter(
        boost::python::allow_null(PyObject_GetIter(obj_ptr)));
    if (!obj_iter.get()) { // must be convertible to an iterator
      PyErr_Clear();
      return 0;
    }
    if (ConversionPolicy::check_convertibility_per_element()) {
      int obj_size = PyObject_Length(obj_ptr);
      if (obj_size < 0) { // must be a measurable sequence
        PyErr_Clear();
        return 0;
      }
      if (!ConversionPolicy::check_size(boost::type<ContainerType>(), obj_size))
        return 0;
      bool is_range = PyRange_Check(obj_ptr);
      std::size_t i = 0;
      if (!all_elements_convertible(obj_iter, is_range, i))
        return 0;
      if (!is_range)
        assert(i == (size_t)obj_size);
    }
    return obj_ptr;
  }

  // This loop factored out by Achim Domma to avoid Visual C++
  // Internal Compiler Error.
  static bool all_elements_convertible(boost::python::handle<> &obj_iter,
                                       bool is_range, std::size_t &i) {
    for (;; i++) {
      boost::python::handle<> py_elem_hdl(
          boost::python::allow_null(PyIter_Next(obj_iter.get())));
      if (PyErr_Occurred()) {
        PyErr_Clear();
        return false;
      }
      if (!py_elem_hdl.get())
        break; // end of iteration
      boost::python::object py_elem_obj(py_elem_hdl);
      boost::python::extract<container_element_type> elem_proxy(py_elem_obj);
      if (!elem_proxy.check())
        return false;
      if (is_range)
        break; // in a range all elements are of the same type
    }
    return true;
  }

  static void
  construct(PyObject *obj_ptr,
            boost::python::converter::rvalue_from_python_stage1_data *data) {
    boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
    void *storage =
        ((boost::python::converter::rvalue_from_python_storage<ContainerType> *)
             data)
            ->storage.bytes;
    new (storage) ContainerType();
    data->convertible = storage;
    ContainerType &result = *((ContainerType *)storage);
    std::size_t i = 0;
    for (;; i++) {
      boost::python::handle<> py_elem_hdl(
          boost::python::allow_null(PyIter_Next(obj_iter.get())));
      if (PyErr_Occurred())
        boost::python::throw_error_already_set();
      if (!py_elem_hdl.get())
        break; // end of iteration
      boost::python::object py_elem_obj(py_elem_hdl);
      boost::python::extract<container_element_type> elem_proxy(py_elem_obj);
      ConversionPolicy::set_value(result, i, elem_proxy());
    }
    ConversionPolicy::assert_size(boost::type<ContainerType>(), i);
  }
};

template <typename ContainerType> struct to_tuple_mapping {
  to_tuple_mapping() {
    boost::python::to_python_converter<ContainerType, to_tuple<ContainerType>
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                                       ,
                                       true
#endif
                                       >();
  }
};

template <typename ContainerType, typename ConversionPolicy>
struct tuple_mapping : to_tuple_mapping<ContainerType> {
  tuple_mapping() { from_python_sequence<ContainerType, ConversionPolicy>(); }
};

template <typename ContainerType> struct tuple_mapping_fixed_size {
  tuple_mapping_fixed_size() {
    tuple_mapping<ContainerType, fixed_size_policy>();
  }
};

template <typename ContainerType> struct tuple_mapping_fixed_capacity {
  tuple_mapping_fixed_capacity() {
    tuple_mapping<ContainerType, fixed_capacity_policy>();
  }
};

template <typename ContainerType> struct tuple_mapping_variable_capacity {
  tuple_mapping_variable_capacity() {
    tuple_mapping<ContainerType, variable_capacity_policy>();
  }
};

template <typename ContainerType> struct tuple_mapping_set {
  tuple_mapping_set() { tuple_mapping<ContainerType, set_policy>(); }
};

#endif
