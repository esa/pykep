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

#ifndef KEPLER_EQUATIONS_H
#define KEPLER_EQUATIONS_H

namespace kep_toolbox {
	inline double kepDE(const double& DE, const double& DM, const double& sigma0, const double& sqrta, const double& a, const double& R){
		return ( (DE - DM + sigma0 / sqrta * (1 - cos(DE)) - (1 - R / a) * sin(DE)));
	}

	inline double d_kepDE(const double& DE, const double& sigma0, const double& sqrta, const double& a, const double& R){
		return ( (1 + sigma0 / sqrta * sin(DE) - (1 - R / a) * cos(DE)) );
	}

	inline double kepDH(const double& DH, const double& DN, const double& sigma0, const double& sqrta, const double& a, const double& R){
		return ( -DH -DN + sigma0/sqrta * (cosh(DH) - 1) + (1 - R / a) * sinh(DH) );
	}

	inline double d_kepDH(const double& DH, const double& sigma0, const double& sqrta, const double& a, const double& R){
		return (-1 + sigma0 / sqrta * sinh(DH) + (1 - R / a) * cosh(DH));
	}
}
#endif // KEPLER_EQUATIONS_H
