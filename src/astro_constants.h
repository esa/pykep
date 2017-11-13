/*****************************************************************************
 *   Copyright (C) 2004-2015 The pykep development team,                     *
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

#ifndef KEP_TOOLBOX_ASTRO_CONSTANTS_H
#define KEP_TOOLBOX_ASTRO_CONSTANTS_H

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <string>
#include <vector>

// Redefining M_PI end M_E to avoid problems in windows
#ifndef M_PI
#define M_PI boost::math::constants::pi<double>()
#endif
#ifndef M_PI_4
#define M_PI_4 boost::math::constants::pi<double>() / 4
#endif
#ifndef M_E
#define M_E boost::math::constants::e<double>()
#endif

#define ASTRO_AU 149597870691.0
#define ASTRO_JR 71492000.0
#define ASTRO_MU_SUN 1.32712440018e20
#define ASTRO_MU_EARTH 398600441800000.0
#define ASTRO_EARTH_VELOCITY 29784.6905
#define ASTRO_EARTH_J2 1.08262668E-03
#define ASTRO_EARTH_RADIUS 6378137
#define ASTRO_DEG2RAD (M_PI / 180.0)
#define ASTRO_RAD2DEG (180.0 / M_PI)
#define ASTRO_DAY2SEC 86400.0 // needs to be a double
#define ASTRO_SEC2DAY (1. / 86400.0)
#define ASTRO_DAY2YEAR (1. / 365.25)
#define ASTRO_G0 9.80665
#define ASTRO_CAVENDISH 73.6687e-11

// This is used as a numerical proceure (e.g. newton-raphson or runge-kutta or
// regula-falsi) stopping criteria
#define ASTRO_TOLERANCE 1e-16

// This is used in the Lambert Problem and in the solution of kepler equation
#define ASTRO_MAX_ITER 50

// This needs to be set to the precision of the boost date library (microseconds
// is default,
// nanoseconds can be set when compiling boosts. Note that the code has not been
// checked in that case)
#define BOOST_DATE_PRECISION 1e-6

namespace kep_toolbox {
// Typedef for fixed size vectors
typedef std::array<double, 3> array3D;
typedef std::array<double, 6> array6D;
typedef std::array<double, 7> array7D;
}

namespace std {
/// Overload stream insertion operator for arrayND
template <size_t N>
inline ostream &operator<<(ostream &os, const std::array<double, N> &v) {
  os << '[';
  for (typename std::array<double, N>::size_type i = 0; i < v.size(); ++i) {
    os << boost::lexical_cast<std::string>(v[i]);
    if (i != v.size() - 1) {
      os << ", ";
    }
  }
  os << ']';
  return os;
}

inline ostream &operator<<(ostream &os, const std::vector<double> &v) {
  os << '[';
  for (kep_toolbox::array3D::size_type i = 0; i < v.size(); ++i) {
    os << boost::lexical_cast<std::string>(v[i]);
    if (i != v.size() - 1) {
      os << ", ";
    }
  }
  os << ']';
  return os;
}
}
#endif // KEP_TOOLBOX_ASTRO_CONSTANTS_H
