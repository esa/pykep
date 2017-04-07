/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
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

#ifndef KEP_TOOLBOX_PLANET_TLE_H
#define KEP_TOOLBOX_PLANET_TLE_H

#include "base.h"
#include "../serialization.h"
#include "../config.h"
#include "../third_party/libsgp4/SGP4.h"
#include "../third_party/libsgp4/Tle.h"

namespace kep_toolbox{ namespace planet{
//using namespace boost::posix_time;

/// A planet from TLE format
/**
 * This class allows to instantiate Earth-orbiting
 * satellites from their Two Line Element format. The ephemerides will then be computed
 * using SGP4/SDP4 orbital model. The third party C++ library SGP4 Satellite Library is
 * used (source code in tp/libsgp4)
 *
 * NOTE: the constant used to initialize the data_members are not the pykep ones, rather the 
 * constants defined in the sgp4lib are used (tp/libsgp4/Globals.h)
 *
 * @see http://celestrak.com/columns/v04n03/#FAQ01
 * @see http://www.danrw.com/sgp4/
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE tle : public base
{
public:
	/**
	 * Construct a planet_tle from two strings containing the two line elements
	 * \param[in] line1 first line
	 * \param[in] line2 second line
	 */
	tle(const std::string & = "1 23177U 94040C   06175.45752052  .00000386  00000-0  76590-3 0    95", const std::string & = "2 23177   7.0496 179.8238 7258491 296.0482   8.3061  2.25906668 97438");
	planet_ptr clone() const;
	std::string human_readable_extra() const;

	double get_ref_mjd2000() const;
	void set_epoch(const unsigned int year, const double day);

private:
	void eph_impl(double mjd2000, array3D &r, array3D &v) const;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<std::string& >(m_line1);
		ar & const_cast<std::string& >(m_line2);
		ar & m_ref_mjd2000;
		boost::serialization::split_member(ar, *this, version);
	}

	template <class Archive>
		void save(Archive &, const unsigned int) const
		{}

	template <class Archive>
		void load(Archive &, const unsigned int)
		{
			// NOTE: the Tle and SGP4 data members are not saved during serialization. Hence, upon loading,
			// we are going to build them again from data. This set up was chosen to avoid implementing
			// serialization of the third-party library libsgp4 objects
			m_tle = Tle("TLE satellite", m_line1, m_line2);
			tm ep = to_tm(epoch(m_ref_mjd2000).get_posix_time());
			double day = ep.tm_yday+ep.tm_hour/24.0+ep.tm_min/(24.0*60.0)+ep.tm_sec/(24.0*60.0*60.0);
            m_tle.setEpoch(1900+ep.tm_year,day);
			m_sgp4_propagator = SGP4(m_tle);
		}

		const std::string m_line1;
		const std::string m_line2;
		Tle m_tle;
		SGP4 m_sgp4_propagator;
		double m_ref_mjd2000;
	};


}} /// End of namespace kep_toolbox

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet::tle)

#endif // KEP_TOOLBOX_PLANET_TLE_H
