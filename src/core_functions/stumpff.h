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

#ifndef STUMPFF_H
#define STUMPFF_H

#include<cmath>

namespace kep_toolbox {

inline double stumpff_s(const double x) {
	if (x > 0)
	{
		return (sqrt(x) - sin(sqrt(x)))/pow(sqrt(x),3);
	}
	else if (x < 0)
	{
		return (sinh(sqrt(-x)) - sqrt(-x))/pow(sqrt(-x),3);
	}
	else
	{
		return (1./6.);
	}
}


inline double stumpff_c(const double x) {
	if (x > 0)
	{
		return (1 - cos(sqrt(x)))/x;
	}
	else if (x < 0)
	{
		return (cosh(sqrt(-x)) - 1)/(-x);
	}
	else
	{
		return 0.5;
	}

}

}

#endif //STUMPFF_H
