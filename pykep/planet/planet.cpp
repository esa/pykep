/*****************************************************************************
 *   Copyright (C) 2004-2018 The pagmo development team,                     *
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

// Workaround for http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#ifdef _WIN32
#include <cmath>
#endif

#include <Python.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/self.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include "../../src/planet/base.h"
#include "../../src/planet/keplerian.h"
#include "../../src/planet/j2.h"
#include "../../src/planet/jpl_low_precision.h"
#include "../../src/planet/mpcorb.h"
#include "../../src/planet/tle.h"
#include "../../src/planet/spice.h"
#include "../../src/planet/gtoc2.h"
#include "../../src/planet/gtoc5.h"
#include "../../src/planet/gtoc6.h"
#include "../../src/planet/gtoc7.h"
#include "../utils.h"
#include "python_base.h"

using namespace boost::python;
using namespace kep_toolbox;

// Wrappers for the ephemerides computations
static inline tuple eph_wrapper1(const kep_toolbox::planet::base &p, const kep_toolbox::epoch &when)
{
	kep_toolbox::array3D r, v;
	p.eph(when.mjd2000(),r,v);
	return boost::python::make_tuple(r,v);
}

static inline tuple eph_wrapper2(const kep_toolbox::planet::base &p, double when)
{
	kep_toolbox::array3D r, v;
	p.eph(when,r,v);
	return boost::python::make_tuple(r,v);
}

// Wrapper to expose planet deriving from base
template <class Planet>
static inline class_<Planet, bases<planet::base> > planet_wrapper(const char *name, const char *descr)
{
	class_<Planet, bases<planet::base> > retval(name,descr,init<const Planet &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Planet>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Planet>);
	retval.def_pickle(python_class_pickle_suite<Planet>());
	return retval;
}

// Wrapper to expose planet deriving from keplerian
template <class Planet>
static inline class_<Planet, bases<planet::base>, bases<planet::keplerian> > planet_kep_wrapper(const char *name, const char *descr)
{
	class_<Planet, bases<planet::base>, bases<planet::keplerian> > retval(name,descr,init<const Planet &>());
	retval.def(init<>());
	retval.def("__copy__", &Py_copy_from_ctor<Planet>);
	retval.def("__deepcopy__", &Py_deepcopy_from_ctor<Planet>);
	retval.def_pickle(python_class_pickle_suite<Planet>());
	return retval;
}

BOOST_PYTHON_MODULE(_planet) {
	// Disable docstring c++ signature to allow sphinx autodoc to work properly
	docstring_options doc_options;
	doc_options.disable_signatures();

	register_ptr_to_python<kep_toolbox::planet::planet_ptr>();

	// Base planet class. This must be python_base as to allow the virtual methods handled in python
	class_<planet::python_base, boost::noncopyable>("_base", "The base class for all planets, it cannot be instantiated as it contains pure virtual methods", init<optional<double, double, double, double, const std::string&> >(
			"pykep.planet._base(mu_central body, mu_self, radius, self_radius, name)\n\n"
			"- mu_central_body: Gravity parameter of the central body (this is not used to compute the ephemerides)\n"
			"- mu_self: 		Gravity parameter of the target\n"
			"- radius: 			Radius of target body\n"
			"- self_radius: 	Safe radius of target body\n"
			"- name: 			Body name\n\n"
		))
		.add_property("safe_radius", &planet::python_base::get_safe_radius, &planet::python_base::set_safe_radius,
			"The body safe radius in [m](distance at which it is safe for spacecraft to fly-by)\n\n"
			"Example::\n\n"
			"  Rs = earth.safe_radius\n"
			"  earth.safe_radius = 1.05"
			)
		.add_property("mu_self", &planet::python_base::get_mu_self, &planet::python_base::set_mu_self,
			"The body gravity parameter in [m^3/s^2]\n\n"
			"Example::\n\n"
			"  mu_pla = earth.mu_self"
			)
		.add_property("mu_central_body", &planet::python_base::get_mu_central_body, &planet::python_base::set_mu_central_body,
			"The central body gravity parameter in [m^3/s^2]\n\n"
			"Example::\n\n"
			"  mu = earth.mu_central_body"
			)
		.add_property("name", &planet::python_base::get_name, &planet::python_base::set_name,
			"The body Name\n\n"
			"Example::\n\n"
			"  name = earth.name"
		)
		.add_property("radius", &planet::python_base::get_radius, &planet::python_base::set_radius,
			"The planet radius in [m]\n\n"
			"Example::\n\n"
			"  R = earth.radius"
		)
		.def("compute_period",&planet::python_base::compute_period,
			"The planet orbital period in [sec]\n\n"
			"Example::\n\n"
			"  T = earth.compute_period()\n"
			"  T = earth.compute_period(epoch(2345.3, 'mjd2000'))"
		)
		.def("osculating_elements", &planet::python_base::compute_elements,
			"Retruns a tuple containing the six osculating keplerian elements a,e,i,W,w,M at the reference epoch\n"
			"(SI units used). If no epoch is passed 0. is assumed\n\n"
			"Example::\n\n"
			"  elem = earth.osculating_elements()\n"
			"  elem = earth.osculating_elements(epoch(2345.3, 'mjd2000'))"
		)
		// Virtual methods that must be reimplemented
		.def("eph",&eph_wrapper1,
			"pykep.planet._base.eph(when)\n\n"
			"- when: a :py:class:`pykep.epoch` indicating the epoch at which the ephemerides are needed, it can also be a double in which case its interpreted as a mjd2000\n\n"
			"Retuns a tuple containing the planet position and velocity in SI units\n\n"
			".. note::\n\n"
			"   This is a pure virtual method and must be reimplemented in the derived class\n\n"
			"Example::\n\n"
			"  r,v = earth.eph(epoch(5433), 'mjd2000')\n"
			"  r,v = earth.eph(5433)"
		)
		.def("eph",&eph_wrapper2," ")
		// Virtual methods that can be reimplemented
		.def("human_readable_extra", &planet::python_base::human_readable_extra,
			"pykep.planet._base.human_readable_extra()\n\n"
			"Retuns a string with extra information on the problem.\n\n"
			".. note::\n\n"
			"   This is a virtual method and can be reimplemented in the derived class. In which case __repr__ will append its returned value\n"
		)
		.def(repr(self))
		.def_pickle(python_class_pickle_suite<planet::python_base>());

		// We start exposing here the classes derived from base
		// 1 - Planets deriving directly from base
		planet_wrapper<planet::keplerian>("keplerian","A planet with Keplerian ephemerides, derives from :py:class:`pykep.planet._base`")
		.def(init<optional<const epoch&, const array6D&, double, double, double , double, const std::string &> >())
		.def(init<const epoch&, const array3D&, const array3D&, double, double, double, double , optional<const std::string &> >())
		.add_property("orbital_elements", &planet::keplerian::get_elements, &planet::keplerian::set_elements,
			"The keplerian elements of the planet\n\n"
			"Example::\n\n"
			"  el = earth.orbital_elements"
		)
		.add_property("ref_epoch", &planet::keplerian::get_ref_epoch, &planet::keplerian::set_ref_epoch,
			"The reference epoch at which the elements are given\n\n"
			"Example::\n\n"
			"  el = earth.ref_epoch"
		)
		.add_property("ref_mjd2000", &planet::keplerian::get_ref_mjd2000, &planet::keplerian::set_ref_mjd2000,
			"The reference epoch at which the elements are given\n\n"
			"Example::\n\n"
			"  el = earth.ref_mjd2000"
		);

        planet_wrapper<planet::j2>("j2","An object with an orbit perturbed by J2, derives from :py:class:`pykep.planet._base`")
        .def(init<optional<const epoch&, const array6D&, double, double, double , double, double, const std::string &> >())
        .def(init<const epoch&, const array3D&, const array3D&, double, double, double, double, double, optional<const std::string &> >())
        .add_property("orbital_elements", &planet::j2::get_elements, &planet::j2::set_elements,
            "The keplerian, osculating, elements of the orbit at the reference epoch\n\n"
            "Example::\n\n"
            "  el = deb1.orbital_elements"
        )
        .add_property("ref_epoch", &planet::j2::get_ref_epoch, &planet::j2::set_ref_epoch,
            "The reference epoch at which the elements are given\n\n"
            "Example::\n\n"
            "  el = deb1.ref_epoch"
        )
        .add_property("ref_mjd2000", &planet::j2::get_ref_mjd2000, &planet::j2::set_ref_mjd2000,
            "The reference epoch at which the elements are given\n\n"
            "Example::\n\n"
            "  el = deb1.ref_mjd2000"
        );

		planet_wrapper<planet::jpl_lp>("jpl_lp","A solar system planet that uses the JPL low-precision ephemerides, derives from :py:class:`pykep.planet._base`")
		.def(init<optional<const std::string &> >(
			"pykep.planet.jpl_lp(name)\n\n"
			"- name: string containing the common planet name (e.g. 'earth', 'venus' etc.)\n\n"
			"Example::\n\n"
			"  earth = planet.jpl_lp('earth')"
		));

		planet_wrapper<planet::tle>("tle","An Earth satellite defined from the TLE format, derives from :py:class:`pykep.planet._base`")
        .def(init<optional<const std::string &, const std::string> >(
			"pykep.planet.tle(line1, line2)\n\n"
			"- line1: string containing the first line of a TLE (69 well formatted chars)\n"
			"- line2: string containing the second line of a TLE (69 well formatted chars)\n\n"
			"Example::\n\n"
			"  line1 = '1 23177U 94040C   06175.45752052  .00000386  00000-0  76590-3 0    95'\n"
			"  line2 = '2 23177   7.0496 179.8238 7258491 296.0482   8.3061  2.25906668 97438'\n"
			"  arianne = planet.tle(line1, line2)"
		))
        .add_property("line1", &planet::tle::get_line1, "Get 1st line of original TLE")
        .add_property("line2", &planet::tle::get_line2, "Get 2nd line of original TLE")
        .add_property("ref_mjd2000", &planet::tle::get_ref_mjd2000, "Get reference epoch (possibly modified from what is in line 1 of TLE)")
        .def("set_epoch",&planet::tle::set_epoch,
		    "PyKep.planet.tle.set_epoch(year, day)\n\n"
		    "- year: unsigned integer specifying the full 4 digit year\n"
		    "- day:  double representing the fractional day as in TLE format.\n\n"
		    "Set the epoch of the TLE to the given date without changing any other elements.\n"
		    "Used to work around the Y2056 bug in TLE definition\n\n");

#ifdef PYKEP_USING_SPICE
		planet_wrapper<planet::spice>("spice","A planet using the eph from the SPICE Toolbox, derives from :py:class:`pykep.planet._base`")
		.def(init<optional<const std::string &, const std::string &, const std::string &, const std::string &, double, double, double, double> >(
			"pykep.planet.spice(target, observer, ref_frame, aberrations, mu_central_body, mu_self, radius, self_radius)\n\n"
			"- target:			Target body (see `NAIF IDs <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html#NAIF%20Object%20ID%20numbers>`_)\n"
			"- observer:		Observer body (see `NAIF IDs <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html#NAIF%20Object%20ID%20numbers>`_)\n"
			"- ref_frame:		The reference frame (see `SPICE supported frames <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#Frames%20Supported%20in%20SPICE>`_)\n"
			"- aberrations: 	Aberration correction type (see spkezr_c docs <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html>`_)\n"
			"- mu_central_body: Gravity parameter of the central body (this is not used to compute the ephemerides). Note that this is needed to compute the orbital elements\n"
			"- mu_self: 		Gravity parameter of the target\n"
			"- radius: 			Radius of target body\n"
			"- self_radius: 	Safe radius of target bodyn\n\n"
			".. note::\n\n"
			"   The presence of the corresponding kernel files are only checked upon call to the ephemerides method. As a consequence this object can still be constructed with invalid names\n"
			"   or spelling mistakes. Only later the ephemerides call will fail throwing an excpetion\n\n"
			".. note::\n\n"
			"   mu_central_body must be set if the period or the orbital elements need to be computed"
			"Example::\n\n"
			"  planet = planet.spice('EARTH', 'SUN', 'ECLIPJ2000', 'NONE', MU_SUN, MU_EARTH, ERATH_R, EARTH_R * 1.05)"
		));
#endif
		// 2 - Planets deriving from keplerian
		planet_kep_wrapper<planet::mpcorb>("mpcorb","A planet from the MPCORB database, derives from :py:class:`pykep.planet.keplerian`")
		.def(init<optional<const std::string &> >(
			"pykep.planet.mpcorb(line)\n\n"
			"- line: a line from the MPCORB database file\n\n"
			"Example::\n\n"
			"  apophis = planet.mpcorb('99942   19.2   0.15 K107N 202.49545  126.41859  204.43202    3.33173  0.1911104  1.11267324   0.9223398  1 MPO164109  1397   2 2004-2008 0.40 M-v 3Eh MPCAPO     C802  (99942) Apophis            20080109')"
		))
		.add_property("H",&kep_toolbox::planet::mpcorb::get_H,
			"The asteroid absolute magnitude. This is assuming an albedo of 0.25 and using the formula at www.physics.sfasu.edu/astro/asteroids/sizemagnitude.html\n"
			"Example::\n\n"
			"  H = apophis.H"
		)
		.add_property("n_observations",&kep_toolbox::planet::mpcorb::get_n_observations,
			"Number of observations made on the asteroid\n"
			"Example::\n\n"
			"  R = apophis.n_observations"
		)
		.add_property("n_oppositions",&kep_toolbox::planet::mpcorb::get_n_oppositions,
			"The planet radius\n"
			"Example::\n\n"
			"  R = apophis.n_oppositions"
		)
		.add_property("year_of_discovery",&kep_toolbox::planet::mpcorb::get_year_of_discovery,
			"The year the asteroid was first discovered. In case the asteroid has been observed only once (n_observations), this number is, instead, the Arc Length in days\n"
			"Example::\n\n"
			"  R = apophis.year_of_discovery"
		);

		planet_kep_wrapper<planet::gtoc2>("gtoc2","An asteroid from gtoc2, derives from :py:class:`pykep.planet.keplerian`")
		.def(init<optional<int> >(
			"pykep.planet.gtoc2(ast_id)\n\n"
			"- ast_id: Construct from a consecutive id from 0 to 910 (Earth)."
			"The order is that of the original data file from JPL\n\n"
			"-    Group 1:   0 - 95\n"
			"-    Group 2:  96 - 271\n"
			"-    Group 3: 272 - 571\n"
			"-    Group 4: 572 - 909\n"
			"-    Earth:   910\n\n"
			"Example::\n\n"
			"  earth_gtoc2 = planet.gtoc2(910)"
		));

		planet_kep_wrapper<planet::gtoc5>("gtoc5","An asteroid from gtoc5, derives from :py:class:`pykep.planet.keplerian`")
		.def(init<optional<int> >(
			"pykep.planet.gtoc5(ast_id)\n\n"
			"- ast_id: a consecutive id from 1 to 7076 (Earth). The order is that of the original"
			"data file distributed by the Russian, see the (`gtoc portal <http://sophia.estec.esa.int/gtoc_portal/>`_). Earth is 7076\n\n"
			"Example::\n\n"
			"  russian_ast = planet.gtoc5(1)"
		));

		planet_kep_wrapper<planet::gtoc6>("gtoc6","A Jupiter moon from gtoc6, derives from :py:class:`pykep.planet.keplerian`")
		.def(init<optional<const std::string &> >(
			"pykep.planet.gtoc6(name)\n\n"
			"- name: string containing the common planet name (e.g. 'io', 'europa', 'callisto' or 'ganymede')\n\n"
			"Example::\n\n"
			"  io = planet.gtoc6('io')"
		));

		planet_kep_wrapper<planet::gtoc7>("gtoc7","An asteroid from gtoc7, derives from :py:class:`pykep.planet.keplerian`")
		.def(init<optional<int> >(
			"pykep.planet.gtoc7(ast_id)\n\n"
			"- ast_id: a consecutive id from 0 (Earth) to 16256. The order is that of the original file , see the (`gtoc portal <http://sophia.estec.esa.int/gtoc_portal/>`_).\n\n"
			"Example::\n\n"
			"  earth = planet.gtoc7(0)"
		));
}
