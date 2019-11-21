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

#include <numeric>
#include <vector>

#include <keplerian_toolbox/astro_constants.hpp>
#include <keplerian_toolbox/core_functions/array3D_operations.hpp>
#include <keplerian_toolbox/core_functions/propagate_lagrangian.hpp>
#include <keplerian_toolbox/exceptions.hpp>
#include <keplerian_toolbox/sims_flanagan/leg_s.hpp>
#include <keplerian_toolbox/sims_flanagan/sc_state.hpp>

namespace kep_toolbox
{
namespace sims_flanagan
{

std::string leg_s::human_readable() const
{
    std::ostringstream s;
    s << *this;
    return s.str();
}

/// Overload the stream operator for kep_toolbox::sims_flanagan::leg
/**
 * Streams out the leg object in a human readable format
 *
 * \param[in] s stream to which the planet will be sent
 * \param[in] in leg to be sent to stream
 *
 * \return reference to s
 *
 */
std::ostream &operator<<(std::ostream &s, const leg_s &in)
{
    s << "Leg in the Sundmann Variable dt = cr^(alpha) ds: " << std::endl << std::endl;
    s << std::setprecision(15);
    s << "Number of segments: " << in.m_throttles.size() << std::endl;
    s << "c: " << in.m_c << std::endl;
    s << "alpha: " << in.m_alpha << std::endl;
    s << "Taylor integration tol: " << in.m_tol << std::endl << std::endl;

    s << in.get_spacecraft() << std::endl;
    s << "Central body gravitational parameter: " << in.get_mu() << std::endl << std::endl;
    s << "Departure date: " << in.get_ti() << ", mjd2000: " << in.get_ti().mjd2000() << std::endl;
    s << "Arrival date: " << in.get_tf() << ", mjd2000: " << in.get_tf().mjd2000() << std::endl;
    s << "Initial mass: " << in.get_xi().get_mass() << " kg" << std::endl;
    s << "Final mass: " << in.get_xf().get_mass() << " kg" << std::endl;
    s << "State at departure: " << in.get_xi() << std::endl;
    s << "State at arrival: " << in.get_xf() << std::endl;

    s << std::endl << "Throttles values: " << std::endl;
    for (size_t i = 0; i < in.m_throttles.size(); i++) {
        s << "\t\t\t" << in.m_throttles[i].get_value()[0] << "\t" << in.m_throttles[i].get_value()[1] << "\t"
          << in.m_throttles[i].get_value()[2] << std::endl;
    }

    try {
        s << std::endl << "Mismatch at the midpoint: ";
        s << in.compute_mismatch_con() << std::endl;
    } catch (...) {
        s << std::endl
          << "Mismatch at the midpoint: ERROR!! COULD NOT CALCULATE THE STATE MISMATCH, CHECK YOUR DATA" << std::endl;
    }

    s << "Throttle magnitude constraints (if <=0 are satisfied): " << in.compute_throttles_con();
    return s;
}
}
} // namespaces
