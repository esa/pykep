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

#include <keplerian_toolbox/planet/spice.hpp>
#include <keplerian_toolbox/exceptions.hpp>

namespace kep_toolbox
{
namespace planet
{
/// Constructor
/** \param[in] mu_central_body gravitational parameter of the central attracting body [m^3/sec^2]
 * \param[in] mu_self gravitational parameter of the body [m^3/sec^2]
 * \param[in] radius body radius [m]
 * \param[in] safe_radius body safe radius [m]
 * \param[in] name body name
*/
spice::spice(const std::string &target, const std::string &observer, const std::string &reference_frame,
             const std::string &aberrations, double mu_central_body, double mu_self, double radius, double safe_radius)
    : base(mu_central_body, mu_self, radius, safe_radius, target + ", " + observer + ", " + reference_frame),
      m_target(target), m_observer(observer), m_reference_frame(reference_frame), m_aberrations(aberrations)
{
    /// Transferring error handling from SPICE to kep_toolbox
    erract_c("SET", 0, const_cast<char *>("RETURN"));
}

/// Polymorphic copy constructor.le::clone() const
planet_ptr spice::clone() const
{
    return planet_ptr(new spice(*this));
}

void spice::eph_impl(double mjd2000, array3D &r, array3D &v) const
{
    SpiceDouble spice_epoch = kep_toolbox::util::epoch_to_spice(mjd2000);
    spkezr_c(m_target.c_str(), spice_epoch, m_reference_frame.c_str(), m_aberrations.c_str(), m_observer.c_str(),
             m_state, &m_lt);
    r[0] = m_state[0] * 1000;
    r[1] = m_state[1] * 1000;
    r[2] = m_state[2] * 1000;
    v[0] = m_state[3] * 1000;
    v[1] = m_state[4] * 1000;
    v[2] = m_state[5] * 1000;
    /// Handling errors
    if (failed_c()) {
        std::ostringstream msg;
        msg << "SPICE cannot compute the ephemerides, have you loaded all needed Kernel files?" << std::endl;
        reset_c();
        throw_value_error(msg.str());
    }
}

/// Extra informations streamed in human readable format
std::string spice::human_readable_extra() const
{
    std::ostringstream s;
    s << "Target planet: " << m_target << std::endl;
    s << "Observer: " << m_observer << std::endl;
    s << "Reference frame: " << m_reference_frame << std::endl;
    s << "Aberrations: " << m_aberrations << std::endl;
    s << "Ephemerides type: SPICE Toolbox" << std::endl;
    return s.str();
}
}
} // namespace

BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::spice)
