/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
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

#ifndef KEP_TOOLBOX_M2E_H
#define KEP_TOOLBOX_M2E_H

#include <cmath>

#include <boost/bind/bind.hpp>

#include <keplerian_toolbox/astro_constants.hpp>
#include <keplerian_toolbox/core_functions/kepler_equations.hpp>
#include <keplerian_toolbox/numerics/newton_raphson.hpp>

using namespace boost::placeholders;

namespace kep_toolbox
{

// mean to eccentric
inline double m2e(const double &M, const double &e)
{
    double E = M + e * sin(M);
    newton_raphson(E, boost::bind(kepE, _1, M, e), boost::bind(d_kepE, _1, e), 100, ASTRO_TOLERANCE);
    return (E);
}
// eccentric to mean
inline double e2m(const double &E, const double &e)
{
    return (E - e * sin(E));
}
// eccentric to true
inline double e2f(const double &E, const double &e)
{
    return 2 * std::atan(std::sqrt((1 + e) / (1 - e)) * std::tan(E / 2));
}
// true to eccentric
inline double f2e(const double &f, const double &e)
{
    return 2 * std::atan(std::sqrt((1 - e) / (1 + e)) * std::tan(f / 2));
}
// gudermannian to true
inline double zeta2f(const double &E, const double &e)
{
    return 2 * std::atan(std::sqrt((1 + e) / (e - 1)) * std::tan(E / 2));
}
// true to gudermannian
inline double f2zeta(const double &zeta, const double &e)
{
    return 2 * std::atan(std::sqrt((e - 1) / (1 + e)) * std::tan(zeta / 2));
}
} // namespace kep_toolbox
#endif // KEP_TOOLBOX_M2E_H
