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

#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/module.hpp>

#ifdef PYKEP_USING_SPICE
#include "../../src/util/spice_utils.h"
#endif

using namespace boost::python;

BOOST_PYTHON_MODULE(_util)
{
    // Disable docstring c++ signature to allow sphinx autodoc to work properly
    docstring_options doc_options;
    doc_options.disable_signatures();

#ifdef PYKEP_USING_SPICE
    // Spice utilities
    def("load_spice_kernel", &kep_toolbox::util::load_spice_kernel,
        "pykep.util.load_spice_kernel(file_name)\n\n"
        "- file_name: string containing the kernel file to load\n\n"
        "Loads the SPICE kernel specified by the filename into memory. \n\n"
        ".. note::\n\n"
        "   The kernel must be in memory before its used, for example, when computing a planets.spice ephemerides\n\n"
        "Example:: \n\n"
        "  util.load_spice_kernel('de432s.bsp')");
#endif
}
