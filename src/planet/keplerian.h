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

#ifndef KEP_TOOLBOX_PLANET_KEPLERIAN_H
#define KEP_TOOLBOX_PLANET_KEPLERIAN_H

#include <string>
#include <vector>

#include "base.h"
#include "../serialization.h"
#include "../config.h"
#include "../exceptions.h"
#include "../epoch.h"

namespace kep_toolbox{ namespace planet {

/// A Keplerian Planet
/**
 * This class allows to instantiate a planet having keplerian ephemerides
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE keplerian : public base
{

static const array6D default_elements;
public:

	keplerian(const epoch& ref_epoch  = kep_toolbox::epoch(0), const array6D& elem = default_elements, double mu_central_body = 0.1, double mu_self = 0.1, double radius = 0.1, double safe_radius = 0.1, const std::string &name = "Unknown");
	keplerian(const epoch& ref_epoch, const array3D& r0, const array3D& v0, double mu_central_body, double mu_self, double radius, double safe_radius, const std::string &name = "Unknown");

	virtual planet_ptr clone() const;
	std::string human_readable_extra() const;

	/** @name Getters */
	//@{
	array6D get_elements() const;
	kep_toolbox::epoch get_ref_epoch() const;
	double get_ref_mjd2000() const;
	double get_mean_motion() const;
	//@}

	/** @name Setters */
	//@{
	void set_elements(const array6D&);
	void set_ref_epoch(const kep_toolbox::epoch&);
	void set_ref_mjd2000(const double &);
	//@}

private:
	void eph_impl(double mjd2000, array3D &r, array3D &v) const;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_r;
		ar & m_v;
		ar & m_keplerian_elements;
		ar & m_mean_motion;
		ar & m_ref_mjd2000;
	}

protected:
	array6D m_keplerian_elements;
	array3D m_r, m_v;
	double m_mean_motion;
	double m_ref_mjd2000;
};

}} /// End of namespace kep_toolbox

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet::keplerian)

#endif // KEP_TOOLBOX_PLANET_KEPLERIAN_H
