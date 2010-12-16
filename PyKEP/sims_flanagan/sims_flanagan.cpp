/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/self.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/overloads.hpp>
#include <boost/utility.hpp>


#include "../../src/keplerian_toolbox.h"
#include "../boost_python_container_conversions.h"

using namespace boost::python;


BOOST_PYTHON_MODULE(_sims_flanagan) {
	// Disable docstring c++ signature to allow sphinx autodoc to work properly
	docstring_options doc_options;
	doc_options.disable_signatures();

	// Exposing the arrays and vectors into python tuples
	to_tuple_mapping<kep_toolbox::array7D>();
	from_python_sequence<kep_toolbox::array7D,fixed_size_policy>();

	// Spacecraft class
	class_<kep_toolbox::sims_flanagan::spacecraft>("spacecraft","Contains design parameters of a NEP spacecraft",
		init<const double &,const double &,const double & >(
			"PyKEP.sims_flanagan.spacecraft(mass,thrust,isp)\n\n"
			"- mass: the spacecraft mass\n"
			"- thrust: the maximum thrust of the spacecraft propulsion system\n"
			"- isp: the specific impulse of the spacecarft propulsion system\n\n"
			"Examples::\n\n"
			" sc = sims_flanagan.spacecraft(4500,0.05,2500)"
		))
		.def("__repr__", &kep_toolbox::sims_flanagan::spacecraft::human_readable)
		.add_property("mass",&kep_toolbox::sims_flanagan::spacecraft::get_mass, &kep_toolbox::sims_flanagan::spacecraft::set_mass,
			"The spacecraft mass\n\n"
			"Example::\n\n"
			"  mass = sc.mass"
			"  sc.mass = 2500"
			)
		.add_property("thrust",&kep_toolbox::sims_flanagan::spacecraft::get_thrust, &kep_toolbox::sims_flanagan::spacecraft::set_thrust,
			"The spacecraft propulsion system maximum thrust\n\n"
			"Example::\n\n"
			"  T = sc.thrust"
			"  sc.thrust = 0.05"
			)
		.add_property("isp",&kep_toolbox::sims_flanagan::spacecraft::get_isp, &kep_toolbox::sims_flanagan::spacecraft::set_isp,
			"The spacecraft propulsion system specific impulse\n\n"
			"Example::\n\n"
			"  T = sc.isp"
			"  sc.isp = 2000"
			);

	// Spacecraft state class
	class_<kep_toolbox::sims_flanagan::sc_state>("sc_state","Represents the state of a spacecraft as position, velocity and mass",
		init<const kep_toolbox::array3D &, const kep_toolbox::array3D &, const double & >(
			"PyKEP.sims_flanagan.sc_state(r,v,m)\n\n"
			"- r: triplet containing the position vector in cartesian coordiantes\n"
			"- v: triplet containing the velocity vector in cartesian coordiantes\n"
			"- m: mass\n"
			"Examples::\n\n"
			" x0 = sims_flanagan.sc_state((1,0,0),(0,1,0),1000)"
		))
		.add_property("r",make_function(&kep_toolbox::sims_flanagan::sc_state::get_position,return_value_policy<copy_const_reference>()),&kep_toolbox::sims_flanagan::sc_state::set_position,
			"The spacecraft position in cartesian coordinates\n\n"
			"Example::\n\n"
			"  r = x0.r"
			"  x0.r = (1.0,2.0,0.0)"
		)
		.add_property("v",make_function(&kep_toolbox::sims_flanagan::sc_state::get_velocity,return_value_policy<copy_const_reference>()),&kep_toolbox::sims_flanagan::sc_state::set_velocity,
			"The spacecraft velocity in cartesian coordinates\n\n"
			"Example::\n\n"
			"  v = x0.v"
			"  x0.v = (0,1.2,0)"
		)
		.add_property("m",&kep_toolbox::sims_flanagan::sc_state::get_mass, &kep_toolbox::sims_flanagan::sc_state::set_mass,
			"The spacecraft mass\n\n"
			"Example::\n\n"
			"  m = x0.m"
			"  x0.m = 1200"
		)
		.def("set", &kep_toolbox::sims_flanagan::sc_state::set_state,
			"Sets the whole spacecraft state at once using a 7 dimensional tuple \n\n"
			"Example::\n\n"
			" x0.set((1,0,0,0,1,0,1000))"
			)
		.def("get", &kep_toolbox::sims_flanagan::sc_state::get_state,
			"Gets the whole spacecraft state at one p[utting it into a 7 dimensional tuple \n\n"
			"Example::\n\n"
			" state = get((1,0,0,0,1,0,1000))"
		);
}
