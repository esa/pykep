#ifndef PYKEP_UTILS_H
#define PYKEP_UTILS_H

#include <Python.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/tuple.hpp>
#include <boost/serialization/serialization.hpp>
#include <sstream>
#include <string>

template <class T>
inline T Py_copy_from_ctor(const T &x)
{
    return T(x);
}

template <class T>
inline T Py_deepcopy_from_ctor(const T &x, boost::python::dict)
{
    return T(x);
}

// Generic pickle suite for C++ classes with default constructor extensible from Python.
// We need to take care of handling the derived class' dict.
// And we need to define the serialize member for the boost::python::wrapper class

// Serialization for python wrapper.
namespace boost
{
namespace serialization
{
template <class Archive, class T>
void serialize(Archive &, boost::python::wrapper<T> &, const unsigned int)
{
}
}
}

template <class T>
struct python_class_pickle_suite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const T &)
    {
        return boost::python::make_tuple();
    }
    static boost::python::tuple getstate(boost::python::object obj)
    {
        T const &x = boost::python::extract<T const &>(obj)();
        std::stringstream ss;
        boost::archive::text_oarchive oa(ss);
        oa << x;
        return boost::python::make_tuple(obj.attr("__dict__"), ss.str());
    }
    static void setstate(boost::python::object obj, boost::python::tuple state)
    {
        using namespace boost::python;
        T &x = extract<T &>(obj)();
        if (len(state) != 2) {
            PyErr_SetObject(PyExc_ValueError, ("expected 2-item tuple in call to __setstate__; got %s" % state).ptr());
            throw_error_already_set();
        }
        // Restore the object's __dict__.
        dict d = extract<dict>(obj.attr("__dict__"))();
        d.update(state[0]);
        // Restore the internal state of the C++ object.
        const std::string str = extract<std::string>(state[1]);
        std::stringstream ss(str);
        boost::archive::text_iarchive ia(ss);
        ia >> x;
    }
    static bool getstate_manages_dict()
    {
        return true;
    }
};

template <class T>
inline void py_cpp_loads(T &x, const std::string &s)
{
    std::stringstream ss(s);
    boost::archive::text_iarchive ia(ss);
    ia >> x;
}

template <class T>
inline std::string py_cpp_dumps(const T &x)
{
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << x;
    return ss.str();
}

#endif // PYKEP_UTILS_H
