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

#ifndef KEP_TOOLBOX_THROTTLE_H
#define KEP_TOOLBOX_THROTTLE_H

#include <iostream>
#include <numeric>

#include <keplerian_toolbox/astro_constants.hpp>
#include <keplerian_toolbox/detail/visibility.hpp>
#include <keplerian_toolbox/epoch.hpp>
#include <keplerian_toolbox/serialization.hpp>

namespace kep_toolbox
{
namespace sims_flanagan
{

/// A single throttle
/**
 * This class models a single throttle in the Sims-Flanagan model. It essentialy contains the cartesian
 * components of one throttle (non dimensional impulse)
 *impulse
 * @author David di Lorenzo
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class KEP_TOOLBOX_DLL_PUBLIC throttle
{
public:
    throttle() : m_start(), m_end()
    {
        m_value[0] = 0;
        m_value[1] = 0;
        m_value[2] = 0;
    }

    throttle(const epoch &_start, const epoch &_end, const array3D &_value)
        : m_start(_start), m_end(_end), m_value(_value)
    {
    }

    const epoch &get_start() const
    {
        return m_start;
    }

    const epoch &get_end() const
    {
        return m_end;
    }

    const array3D &get_value() const
    {
        return m_value;
    }

    void set_start(const epoch &e)
    {
        m_start = e;
    }
    void set_end(const epoch &e)
    {
        m_end = e;
    }
    void set_value(const array3D &e)
    {
        m_value = e;
    }

    double get_norm() const
    {
        return std::sqrt(std::inner_product(m_value.begin(), m_value.end(), m_value.begin(), 0.));
    }
    std::string human_readable() const
    {
        std::ostringstream s;
        s << "start = " << m_start << std::endl;
        s << "value = " << m_value << std::endl;
        s << "end = " << m_end;
        return s.str();
    }

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_start;
        ar &m_end;
        ar &m_value;
    }

    epoch m_start;
    epoch m_end;
    array3D m_value;
};
} // namespace sims_flanagan
} // namespace kep_toolbox

#endif // KEP_TOOLBOX_THROTTLE_H
