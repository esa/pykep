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

#ifndef M2E_H
#define M2E_H

#include<boost/bind.hpp>
#include<cmath>

#include"../astro_constants.h"
#include"../core_functions/kepler_equations.h"
#include"../numerics/newton_raphson.h"

namespace kep_toolbox {

	inline double m2e(const double& M, const double & eccentricity) {
		double E = M + eccentricity * cos(M);
		newton_raphson(E,boost::bind(kepE,_1,M, eccentricity),boost::bind(d_kepE,_1, eccentricity),100,ASTRO_TOLERANCE);
		return (E);
	}
	inline double e2m(const double& E, const double & eccentricity) {
		return (E - eccentricity * sin (E) );
	}
}
#endif // M2E_H
