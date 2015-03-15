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
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>

#include "../../src/utils/spice_utils.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(_utils) {
	// Disable docstring c++ signature to allow sphinx autodoc to work properly
	docstring_options doc_options;
	doc_options.disable_signatures();
	
	// Spice utilities
	def("load_spice_kernel",&kep_toolbox::utils::load_spice_kernel,
			  "PyKEP.utils.load_spice_kernel(file_name)\n\n"
			  "- file_name: kernel file to load\n"
			  "Loads the SPICE kernel specified by the filename. \n"
			  "Example:: \n\n"
			  "  t = utils.load_spice_kernel('de432s.bsp')"
	);
}
