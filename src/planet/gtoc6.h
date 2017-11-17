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

#ifndef KEP_TOOLBOX_PLANET_GTOC_6_H
#define KEP_TOOLBOX_PLANET_GTOC_6_H

#include "../config.h"
#include "../serialization.h"
#include "keplerian.h"

namespace kep_toolbox
{
namespace planet
{

/// A Jupiter moon from GTOC6 (keplerian)
/**
 * This class allows to instantiate moons of
 * the Jupiter system by referring to their common names. Ephemerides are Keplerian
 * and elements are those defined for the GTOC6 competition
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE gtoc6 : public keplerian
{
public:
    gtoc6(const std::string & = "io");
    planet_ptr clone() const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &boost::serialization::base_object<keplerian>(*this);
    }
};
}
} /// End of namespaces

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet::gtoc6);

#endif // KEP_TOOLBOX_PLANET_GTOC_6_H
