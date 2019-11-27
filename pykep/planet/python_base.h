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

#ifndef KEP_TOOLBOX_PYTHON_PLANET_BASE_H
#define KEP_TOOLBOX_PYTHON_PLANET_BASE_H

#include <boost/python/wrapper.hpp>
#include <string>

#include <keplerian_toolbox/detail/visibility.hpp>
#include <keplerian_toolbox/exceptions.hpp>
#include <keplerian_toolbox/planet/base.hpp>
#include <keplerian_toolbox/serialization.hpp>

// The fact we have multiple libarries core.pyd, planet.pyd etc. makes it
// necessary to have this visibility rule.
#if defined(_WIN32) || defined(__CYGWIN__)

#define PYKEP_DLL_PUBLIC __declspec(dllexport)

#endif

namespace kep_toolbox
{
namespace planet
{

// Wrapper for exporting the planet::base into python
class PYKEP_DLL_PUBLIC python_base : public base, public boost::python::wrapper<base>
{
public:
    /// Same constructor as plates::base
    python_base(double mu_central_body = 0.1, double mu_self = 0.1, double radius = 0.1, double safe_radius = 0.1,
                const std::string &name = "Unknown")
        : base(mu_central_body, mu_self, radius, safe_radius, name), boost::python::wrapper<base>()
    {
    }

    /// Clone pure virtual
    planet_ptr clone() const
    {
        planet_ptr retval = this->get_override("__get_deepcopy__")();
        if (!retval) {
            throw_value_error(
                "algorithms's __get_deepcopy__() method returns a NULL pointer, please check the implementation");
        }
        return retval;
    }

    /// Human readable virtual
    std::string human_readable_extra() const
    {
        if (boost::python::override f = this->get_override("human_readable_extra")) {
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1700)
            return boost::python::call<std::string>(this->get_override("human_readable_extra").ptr());
#else
            return f();
#endif
        }
        return base::human_readable_extra();
    }
    std::string default_human_readable_extra() const
    {
        return this->base::human_readable_extra();
    }

protected:
    /// Pure virtual ephemerides
    void eph_impl(double mjd2000, array3D &r, array3D &v) const
    {
        if (boost::python::override f = this->get_override("eph_impl")) {
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1700)
            boost::python::call<void>(f.ptr());
#else
            f(mjd2000, r, v);
            return;
#endif
        }
        throw_value_error("ephemerides have not been implemented!!");
    }

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &boost::serialization::base_object<base>(*this);
        ar &boost::serialization::base_object<boost::python::wrapper<base>>(*this);
    }
};
}
} // namespaces

BOOST_CLASS_EXPORT(kep_toolbox::planet::python_base)

#endif
