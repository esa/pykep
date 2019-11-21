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
#ifndef KEP_TOOLBOX_THREE_IMPULSES_APPROX_H
#define KEP_TOOLBOX_THREE_IMPULSES_APPROX_H

#include <cmath>

#include <keplerian_toolbox/astro_constants.hpp>
#include <keplerian_toolbox/exceptions.hpp>
#include <keplerian_toolbox/planet/base.hpp>

/// Computes the orbital transfer cost using a 3 impulse model
/**
 * This function can be used to compute the orbital transfer cost between two generic keplerian orbits. It makes use
 * of a three impulse approximation, where the three impulses are used (in the best sequence) to match apogee, perigee
 * and inclination. The argument of perigee is not accounted for and thus the approximation is breaking down at high
 * eccentircities
 *
 * \param[in] a1  departure semi-major axis
 * \param[in] e1  departure eccentricity
 * \param[in] i1  departure inclination [rad]
 * \param[in] W1 departure Right Ascension of the Ascending Node (RAAN) [rad]
 * \param[in] a2  target semi-major axis
 * \param[in] e2  target eccentricity
 * \param[in] i3  target inclination [rad]
 * \param[in] W2 target Right Ascension of the Ascending Node (RAAN) [rad]
 * \param[in] mu  gravitational parameter of the central body
 *
 * \param[out] required DV  units compatible to a, mu]
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

namespace kep_toolbox
{

inline double three_impulses_approx(double a1, double e1, double i1, double W1, double a2, double e2, double i2,
                                    double W2, double mu)
{
    // Computes the apocenters and pericenters of target and departure
    double ra1 = a1 * (1 + e1);
    double ra2 = a2 * (1 + e2);
    double rp1 = a1 * (1 - e1);
    double rp2 = a2 * (1 - e2);

    // Computes the relative inclination between orbits
    double cosiREL = std::cos(i1) * std::cos(i2) + std::sin(i1) * std::sin(i2) * std::cos(W1 - W2);

    if (ra1 > ra2) // Strategy is Apocenter-Pericenter
    {
        // Change Inclination + Change Pericenter with a single burn at a2
        double Vi = std::sqrt(mu * (2.0 / ra1 - 1.0 / a1));
        double Vf = std::sqrt(mu * (2.0 / ra1 - 2.0 / (rp2 + ra1)));
        double DV1 = std::sqrt(Vi * Vi + Vf * Vf - 2.0 * Vi * Vf * cosiREL);
        // We change the apocenter with a single burn at p2
        double DV2
            = std::sqrt(mu) * std::abs(std::sqrt(2.0 / rp2 - 2.0 / (rp2 + ra1)) - std::sqrt(2.0 / rp2 - 1.0 / a2));
        return DV1 + DV2;
    } else // (ra1<ra2) Strategy is Pericenter-Apocenter
    {
        // We reverse the manouvres (could as well reverse the parameters)
        double DV1 = std::sqrt(mu)
                     * std::abs(std::sqrt(2.0 / rp1 - 2.0 / (rp1 + ra1)) - std::sqrt(2.0 / rp1 - 2.0 / (rp1 + ra2)));
        double Vi = std::sqrt(mu * (2.0 / ra2 - 2.0 / (rp1 + ra2)));
        double Vf = std::sqrt(mu * (2.0 / ra2 - 1.0 / a2));
        double DV2 = std::sqrt(std::abs(Vi * Vi + Vf * Vf - 2.0 * Vi * Vf * cosiREL));
        return DV1 + DV2;
    }
}

/// Computes the orbital transfer cost using a 3 impulse model
/**
 * This function can be used to compute the orbital transfer cost between two generic keplerian orbits. It makes use
 * of a three impulse approximation, where the three impulses are used (in the best sequence) to match apogee, perigee
 * and inclination. The argument of perigee is not accounted for and thus the approximation is breaking down at high
 * eccentircities
 *
 * \param[in] pl1 departure kep_toolbox::planet
 * \param[in] pl2 target kep_toolbox::planet
 *
 * \param[out] required DV [m/s]
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
inline double three_impulses_approx(const planet::base &pl1, const planet::base &pl2)
{
    if (pl1.get_mu_central_body() != pl2.get_mu_central_body()) {
        throw_value_error(
            "The departure and arrival planets do not have the same central body gravitational parameter");
    }
    array6D el1 = pl1.compute_elements();
    array6D el2 = pl2.compute_elements();
    return three_impulses_approx(el1[0], el1[1], el1[2], el1[3], el2[0], el2[1], el2[2], el2[3],
                                 pl1.get_mu_central_body());
}

/// Computes the orbital transfer cost using a 3 impulse model
/**
 * This function can be used to compute the orbital transfer cost between two generic keplerian orbits. It makes use
 * of a three impulse approximation, where the three impulses are used (in the best sequence) to match apogee, perigee
 * and inclination. The argument of perigee is not accounted for and thus the approximation is breaking down at high
 * eccentircities
 *
 * This overload allows to compute the orbital parameters at a given epoch, useful for non keplerian ephemerides as the
 * ones of
 * kep_toolbox::planet_tle
 *
 * \param[in] pl1 departure kep_toolbox::planet
 * \param[in] pl2 target kep_toolbox::planet
 * \param[in] pl1 departure kep_toolbox::epoch (note that this is only used to compute the orbital parameters of the
 * departure orbit)
 * \param[in] pl2 target kep_toolbox::epoch
 *
 * \param[out] required DV [m/s]
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
inline double three_impulses_approx(const planet::base &pl1, const planet::base &pl2, const epoch &ep1,
                                    const epoch &ep2)
{
    if (pl1.get_mu_central_body() != pl2.get_mu_central_body()) {
        throw_value_error(
            "The departure and arrival planets do not have the same central body gravitational parameter");
    }
    array6D el1 = pl1.compute_elements(ep1);
    array6D el2 = pl2.compute_elements(ep2);
    return three_impulses_approx(el1[0], el1[1], el1[2], el1[3], el2[0], el2[1], el2[2], el2[3],
                                 pl1.get_mu_central_body());
}

} // namespace kep_toolbox

#endif // KEP_TOOLBOX_THREE_IMPULSES_APPROX_H
