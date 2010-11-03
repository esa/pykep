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
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/self.hpp>
#include <boost/utility.hpp>

#include "../src/keplerian_toolbox.h"
#include "boost_python_container_conversions.h"
//#include "exceptions.h"
//#include "../utils.h"

using namespace boost::python;

static inline tuple planet_get_eph(const kep_toolbox::planet &p, const kep_toolbox::epoch &when)
{
	kep_toolbox::array3D r, v;
	p.get_eph(when,r,v);
	return boost::python::make_tuple(r,v);
}

BOOST_PYTHON_MODULE(_PyKEP) {
	// Translate exceptions for this module.
        //translate_exceptions();

	to_tuple_mapping<kep_toolbox::array6D>();
	from_python_sequence<kep_toolbox::array6D,fixed_size_policy>();
	to_tuple_mapping<kep_toolbox::array3D>();
	from_python_sequence<kep_toolbox::array3D,fixed_size_policy>();
	to_tuple_mapping<std::vector<kep_toolbox::array3D> >();
	from_python_sequence<std::vector<kep_toolbox::array3D>,variable_capacity_policy>();
	enum_<kep_toolbox::epoch::type>("epoch_type")
		.value("MJD", kep_toolbox::epoch::MJD)
		.value("MJD2000", kep_toolbox::epoch::MJD2000)
		.value("JD", kep_toolbox::epoch::JD);
	class_<kep_toolbox::epoch>("epoch","Epoch class.",init<const double &,kep_toolbox::epoch::type>())
		.def(repr(self));

	// Base planet class.
	class_<kep_toolbox::planet,boost::noncopyable>("_planet",no_init)
		.def("get_eph",&planet_get_eph)
		.def(repr(self))
		.def("get_elements",&kep_toolbox::planet::get_elements);
	
	// Solar system planet.
	class_<kep_toolbox::asteroid_gtoc5,bases<kep_toolbox::planet> >("asteroid_gtoc5",init<optional<const int &> >());

	// Lambert.
	class_<kep_toolbox::lambert_problem>("lambert_problem",init<const kep_toolbox::array3D &, const kep_toolbox::array3D &, const double &, optional<const double &, const int &> >())
		.def("get_v1",&kep_toolbox::lambert_problem::get_v1,return_value_policy<copy_const_reference>())
		.def("get_v2",&kep_toolbox::lambert_problem::get_v2,return_value_policy<copy_const_reference>())
		.def("get_a",&kep_toolbox::lambert_problem::get_a,return_value_policy<copy_const_reference>())
		.def("get_p",&kep_toolbox::lambert_problem::get_p,return_value_policy<copy_const_reference>())
		.def(repr(self))
		.def("is_reliable",&kep_toolbox::lambert_problem::is_reliable)
		.def("get_Nmax",&kep_toolbox::lambert_problem::get_Nmax);

	register_ptr_to_python<kep_toolbox::planet_ptr>();
}
