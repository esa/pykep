/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/self.hpp>
#include <boost/python/operators.hpp>

#include "python_base.h"
#include "../../src/planets/base.h"
#include "../../src/planets/keplerian.h"
#include "../../src/planets/jpl_low_precision.h"
#include "../../src/planets/mpcorb.h"
#include "../utils.h"

using namespace boost::python;
using namespace kep_toolbox;

// Wrappers for the ephemerides computations
static inline tuple eph_wrapper1(const kep_toolbox::planets::base &p, const kep_toolbox::epoch &when)
{
	kep_toolbox::array3D r, v;
	p.eph(when.mjd2000(),r,v);
	return boost::python::make_tuple(r,v);
}

static inline tuple eph_wrapper2(const kep_toolbox::planets::base &p, double when)
{
	kep_toolbox::array3D r, v;
	p.eph(when,r,v);
	return boost::python::make_tuple(r,v);
}

// Wrapper to expose planets deriving from base
template <class Planet>
static inline class_<Planet, bases<planets::base> > planet_wrapper(const char *name, const char *descr)
{
	class_<Planet, bases<planets::base> > retval(name,descr,init<const Planet &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Planet>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Planet>);
	retval.def_pickle(generic_pickle_suite<Planet>());
	return retval;
}

// Wrapper to expose planets deriving from keplerian
template <class Planet>
static inline class_<Planet, bases<planets::keplerian> > planet_kep_wrapper(const char *name, const char *descr)
{
	class_<Planet, bases<planets::keplerian> > retval(name,descr,init<const Planet &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Planet>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Planet>);
	retval.def_pickle(generic_pickle_suite<Planet>());
	return retval;
}

BOOST_PYTHON_MODULE(_planets) {
	// Disable docstring c++ signature to allow sphinx autodoc to work properly
	docstring_options doc_options;
	doc_options.disable_signatures();

	// Base planet class. This must be python_base as to allow the virtual methods handled in python
	class_<planets::python_base, boost::noncopyable>("_base", "All planets inherit from this class", init<>())
		.add_property("safe_radius", &planets::python_base::get_safe_radius, &planets::python_base::set_safe_radius,
			"The planet safe radius (distance at which it is safe for spacecraft to fly-by)\n\n"
			"Example::\n\n"
			"  Rs = earth.safe_radius"
			"  earth.safe_radius = 1.05"
			)
		.add_property("mu_self", &planets::python_base::get_mu_self, &planets::python_base::set_mu_self,
			"The planet radius\n\n"
			"Example::\n\n"
			"  mu_pla = earth.mu_self"
			)
		.add_property("mu_central_body", &planets::python_base::get_mu_central_body, &planets::python_base::set_mu_central_body,
			"The planet radius\n\n"
			"Example::\n\n"
			"  mu = earth.mu_central_body"
			)
		.add_property("name", &planets::python_base::get_name, &planets::python_base::set_name,
			"The planet Name\n\n"
			"Example::\n\n"
			"  name = earth.name()"
		)
		.add_property("radius", &planets::python_base::get_radius, &planets::python_base::set_radius,
			"The planet radius\n\n"
			"Example::\n\n"
			"  R = earth.radius"
		)
		.add_property("period",&planets::python_base::compute_period,
			"The planet orbital period\n\n"
			"Example::\n\n"
			"  T = earth.period"
		)
		.def("osculating_elements", &planets::python_base::compute_elements,
			"Retruns a tuple containing the six osculating keplerian elements a,e,i,W,w,M at the reference epoch\n"
			"(SI units used). If no epoch is passed 0. is assumed\n\n"
			"Example::\n\n"
			"  elem = earth.osculating_elements()\n"
			"  elem = earth.osculating_elements(epoch(2345.3, 'mjd2000'))"
		)
		.def("eph",&eph_wrapper1,
			"PyKEP.planet.eph(when)\n\n"
			"- when: a :py:class:`PyKEP.epoch` indicating the epoch at which the ephemerides are needed\n"
			"        can also be a double in which case its interpreted as a mjd2000\n\n"
			"Retuns a tuple containing the planet position and velocity in SI units\n\n"
			"Example::\n\n"
			"  r,v = earth.eph(epoch(5433), 'mjd2000')\n"
			"  r,v = earth.eph(5433)"
		)
		.def("eph",&eph_wrapper2," ")
		// Virtual methods can be reimplemented
		.def("human_readable_extra", &planets::python_base::human_readable_extra, &planets::python_base::default_human_readable_extra)
		.def(repr(self))
		.def_pickle(python_class_pickle_suite<planets::python_base>());

		// We start exposing here the classes derived from base
		// Keplerian planet
		planet_wrapper<planets::keplerian>("keplerian","A planet with Keplerian ephemerides")
		.def(init<optional<const epoch&, const array6D&, double, double, double , double, const std::string &> >())
		.def(init<const epoch&, const array3D&, const array3D&, double, double, double, double , optional<const std::string &> >())
		.add_property("orbital_elements", &planets::keplerian::get_elements, &planets::keplerian::set_elements,
			"The keplerian elements of the planet\n\n"
			"Example::\n\n"
			"  el = earth.orbital_elements"
		)
		.add_property("ref_mjd2000", &planets::keplerian::get_ref_epoch, &planets::keplerian::set_ref_epoch,
			"The reference epoch at which the elements are given\n\n"
			"Example::\n\n"
			"  el = earth.ref_mjd2000"
		);

		planet_wrapper<planets::jpl_low_precision>("jpl_low_precision","A solar system planet that uses the JPL low-precision ephemerides")
		.def(init<optional<const std::string &> >());

		planet_kep_wrapper<planets::mpcorb>("mpcorb","A planet from the MPCORB database")
		.def(init<optional<const std::string &> >());
}
