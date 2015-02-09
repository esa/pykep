/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://keptoolbox.sourceforge.net/index.html                            *
 *   http://keptoolbox.sourceforge.net/credits.html                          *
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

#ifndef KEP_TOOLBOX_PLANET_H
#define KEP_TOOLBOX_PLANET_H

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>

#include "../serialization.h"
#include "../config.h"
#include "../exceptions.h"
#include "../astro_constants.h"
#include "../epoch.h"

namespace kep_toolbox{

// Forward declaration.
class planet;

typedef boost::shared_ptr<planet> planet_ptr;


/// A Generic Planet
/**
 * This class is intended to represent bodies in a keplerian orbit around a primary body. A planet
 * is described by its gravity, the gravity of the body it is orbiting, its radius, its safe radius
 * and its orbit. The orbit is internally represented by the planet cartesian coordinates
 * at a given reference epoch.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE planet
{
	friend std::ostream &operator<<(std::ostream &, const planet &);
public:
	/// Default dummy constructor
	planet():mean_motion(0),ref_mjd2000(0), radius(0), safe_radius(0), mu_self(0), mu_central_body(0), cached_epoch(epoch(0)), cached_r(array3D()), cached_v(array3D()) {};
	/// Constructor
	/**
		* Constructs a planet from its elements and its phyisical parameters
		* \param[in] ref_epoch epoch to which the elements are referred to
		* \param[in] elem A STL vector containing the keplerian parameters (a,e,i,Om,om,M). (SI units)
		* \param[in] mu_central_body The gravitational parameter of the attracting body (SI units)
		* \param[in] mu_self The gravitational parameter of the planet (SI units)
		* \param[in] radius radius of the planet (SI units)
		* \param[in] safe_radius mimimual radius that is safe during a fly-by of the planet (SI units)
		* \param[in] name C++ string containing the planet name. Default value is "Unknown"
		*/
	planet(const epoch& ref_epoch, const array6D& elem, const double & mu_central_body, const double &mu_self, const double &radius, const double &safe_radius, const std::string &name = "Unknown");

	/// Constructor
	/**
		* Constructs a planet from its position at epoch and its phyisical parameters
		* \param[in] ref_epoch epoch to which the elements are referred to
		* \param[in] r0 A STL vector containing the planet position
		* \param[in] v0 A STL vector containing the planet velociy
		* \param[in] mu_central_body The gravitational parameter of the attracting body (SI units)
		* \param[in] mu_self The gravitational parameter of the planet (SI units)
		* \param[in] radius radius of the planet (SI units)
		* \param[in] safe_radius mimimual radius that is safe during a fly-by of the planet (SI units)
		* \param[in] name C++ string containing the planet name. Default value is "Unknown"
	*/
	planet(const epoch& ref_epoch, const array3D& r0, const array3D& v0, const double & mu_central_body, const double &mu_self, const double &radius, const double &safe_radius, const std::string &name = "Unknown");

	/// Polymorphic copy constructor.
	virtual planet_ptr clone() const;
	virtual ~planet();
	/** @name Getters */
	//@{
	/// Gets the planet position and velocity
	/**
		* \param[in] when Epoch in which ephemerides are required
		* \param[out] r Planet position at epoch (SI units)
		* \param[out] v Planet velocity at epoch (SI units)
		*/
	virtual void get_eph(const epoch& when, array3D &r, array3D &v) const;

	/// Getter for the central body gravitational parameter
	/**
	 * Gets the gravitational parameter of the central body
	 *
	 * @return mu_central_body (SI Units)
	 */
	double get_mu_central_body() const {return mu_central_body;}

	/// Getter for the planet gravitational parameter
	/**
	 * Gets the gravitational parameter of the planet
	 *
	 * @return mu_self (SI Units)
	 */
	double get_mu_self() const {return mu_self;}

	/// Getter for the planet radius
	/**
	 * Gets the radius of the planet
	 *
	 * @return const reference to radius (SI Units)
	 */
	double get_radius() const {return radius;}

	/// Getter for the planet safe-radius
	/**
	 * Gets the safe-radius of the planet. This is intended to be the minimum distance
	 * from the planet center that is safe ... It may be used, for example,  during fly-bys as a constarint
	 * on the spacecraft trajectory
	 *
	 * @return const reference to safe_radius (SI Units)
	 */
	double get_safe_radius() const {return safe_radius;}

	//@}

	/** @name Setters */
	/// Setter for the planet safe-radius
	/**
	 * Sets the safe-radius of the planet. This is intended to be the minimum distance
	 * from the planet center that is safe ... It may be used, for example,  during fly-bys as a constarint
	 * on the spacecraft trajectory
	 *
	 * \param[in] safe_radius Minimum allowed planetary distance (in planetary radius units)
	 * \throws value_error if safe_radius in < 1
	 */
	void set_safe_radius(double sr) {
		if (sr <=1) {
			throw_value_error("Trying to set a safe_radius that is smaller than the planetary radius");
		}
		safe_radius = sr * get_radius();
	}

	/** @name Ephemerides calculations */
	//@{
	/// Returns the planet position
	/**
	 * \param[in] when Epoch in which position is requested
	 *
	 * @return a boost array containing the planet position in epoch (SI Units)
	 */
	array3D get_position(const epoch& when) const;

	/// Returns the planet velocity
	/**
	  * \param[in] when Epoch in which velocity is requested
	  *
	  * @return a boost array containing the planet velocity in epoch (SI Units)
	  */

	array3D get_velocity(const epoch& when) const;

	/// Returns the planet orbital elements at a given epoch (a,e,i,Om,om,M)
	/**
	 * \param[in] when Epoch in which orbital elements are required
	 *
	 * @return a boost array containing the planet elements in epoch (SI Units) (a,e,i,Om,om,M). Mean anomaly is
	 * returned in range 0,2*pi
	 */
	array6D get_elements(const epoch& when) const;

	/// Returns the planet orbital elements at the reference epoch (a,e,i,Om,om,M)
	/**
	 * @return a boost array containing the planet elements in epoch (SI Units) (a,e,i,Om,om,M). Mean anomaly is
	 * returned in range 0,2*pi
	 */
	array6D get_elements() const;

	/// Returns the planet name
	std::string get_name() const;

	/// Returns the reference epoch
	epoch get_ref_epoch() const {return epoch(ref_mjd2000);}

	//@}

	/// Computes the orbital period
	double compute_period() const;

protected:
	/// Builds the planet assiging all values to members
	/**
	* Constructs a planet from its elements and its phyisical parameters
	* \param[in] ref_epoch epoch to which the elements are referred to
	* \param[in] elem A STL vector containing the keplerian parameters (a,e,i,Om,om,M). AU and degrees are assumed
	* \param[in] mu_central_body The gravitational parameter of the attracting body (SI units)
	* \param[in] mu_self The gravitational parameter of the planet (SI units)
	* \param[in] radius radius of the planet (SI units)
	* \param[in] safe_radius mimimual distance that is safe during a fly-by of the planet (SI units)
	* \param[in] name C++ string containing the planet name. Default value is "Unknown"
	*/
	void build_planet(const epoch& ref_epoch, const array6D& elem, const double & mu_central_body, const double &mu_self, const double & radius, const double & safe_radius, const std::string &name = "Unknown");
private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & keplerian_elements;
		ar & mean_motion;
		ar & ref_mjd2000;
		ar & radius;
		ar & safe_radius;
		ar & mu_self;
		ar & mu_central_body;
		ar & cached_epoch;
		ar & cached_r;
		ar & cached_v;
		ar & m_name;
	}
// Serialization code (END)

protected:
	array6D keplerian_elements;
	double mean_motion;
	double ref_mjd2000;
	double radius;
	double safe_radius;
	double mu_self;
	double mu_central_body;
	mutable epoch cached_epoch;
	mutable array3D cached_r;
	mutable array3D cached_v;

	std::string m_name;

};

__KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &s, const planet &body);
} /// End of namespace kep_toolbox

// Serialization code
BOOST_SERIALIZATION_ASSUME_ABSTRACT(kep_toolbox::planet)
// Serialization code (END)

#endif // KEP_TOOLBOX_PLANET_H
