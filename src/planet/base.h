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

#ifndef KEP_TOOLBOX_PLANET_BASE_H
#define KEP_TOOLBOX_PLANET_BASE_H

#include <boost/shared_ptr.hpp>
#include <string>

#include "../epoch.h"
#include "../exceptions.h"
#include "../serialization.h"
#include "../config.h"
#include "../astro_constants.h"

namespace kep_toolbox{ namespace planet {

// Forward declaration.
class __KEP_TOOL_VISIBLE base;
typedef boost::shared_ptr<base> planet_ptr;


/// Base class for planet
/**
 * A base planet in pykep is defined by its name, its radius, its safe radius (i.e. how close to it its considered to be safe)
 * its gravity parameter and the gravitational parameter of the attracting body. All classes deriving from planet::base
 * will have to implement the planet ephemerides planet::base::eph_impl which is the core method of this class.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE base
{
public:
	base(double mu_central_body = 0.1, double mu_self = 0.1, double radius = 0.1, double safe_radius = 0.1, const std::string &name = "Unknown");
	virtual planet_ptr clone() const = 0;
	virtual ~base() {};

	/// Ephemerides methods
	void eph(const epoch& when, array3D &r, array3D &v) const; 
	void eph(const double mjd2000, array3D &r, array3D &v) const;

	/// Simple basic keplerian mechanics computations
	array6D compute_elements(const epoch& when = kep_toolbox::epoch(0)) const;
	double compute_period(const epoch& when = kep_toolbox::epoch(0)) const;

	/// Methods to stream the base and derived class in a human readable format
	std::string human_readable() const;
	virtual std::string human_readable_extra() const {return std::string();}

	/** @name Getters */
	//@{
	double get_mu_central_body() const;
	double get_mu_self() const;
	double get_radius() const;
	double get_safe_radius() const;
	std::string get_name() const;
	//@}

	/** @name Setters */
	//@{
	void set_safe_radius(double sr);
	void set_mu_central_body(double mu);
	void set_mu_self(double mu);
	void set_radius(double radius);
	void set_name(const std::string &radius);
	//@}

protected:
	virtual void eph_impl(double mjd2000, array3D &r, array3D &v) const = 0;
	
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{	
		ar & m_mu_central_body;	
		ar & m_mu_self;
		ar & m_radius;
		ar & m_safe_radius;
		ar & m_name;
	}

	double m_mu_central_body;
	double m_mu_self;
	double m_radius;
	double m_safe_radius;
	std::string m_name;
};

__KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &s, const base &body);
}} /// End of namespace kep_toolbox planet

BOOST_SERIALIZATION_ASSUME_ABSTRACT(kep_toolbox::planet::base)

#endif // KEP_TOOLBOX_PLANET_BASE_H
