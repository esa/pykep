/*****************************************************************************
*   Copyright (C) 2004-2009 The PaGMO development team,                     *
*   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
*   http://apps.sourceforge.net/mediawiki/pagmo                             *
*   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
*   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#ifndef SC_STATE_H
#define SC_STATE_H

#include <boost/lexical_cast.hpp>

#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
#include "../serialization.h"
#endif

#include "../astro_constants.h"

namespace kep_toolbox {
namespace sims_flanagan{

/// Spacecraft state
/** The state of a spacecraft as defined by its position, velocity and mass. The class is basically
 * a container of an array7D providing structured access to its components. In particular the
 * structure is as follows: \f$ \mathbf x = [\mathbf r, \mathbf v, m]\f$.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
class sc_state
{
public:
	/** @name Constructors*/
	//@{
	sc_state() {
		position[0] = 0;
		position[1] = 0;
		position[2] = 0;
		velocity[0] = 0;
		velocity[1] = 0;
		velocity[2] = 0;
		mass = 0;
	}

	/// Constructor.
	/**
	 * Constructs an sc_state object from position, velocity and mass.
	 *
	 * \param[in] r A three dimensional array containing the position
	 * \param[in] v A three dimensional array containing the velocity
	 * \param[in] m The mass
	 */
	sc_state(const array3D &r, const array3D &v, const double &m) {
	    position = r;
	    velocity = v;
	    mass = m;
	}
	//@}

	/** @name Getters*/
	//@{
	/// Gets the position
	const array3D& get_position() const { return position; }
	/// Gets the velocity
	const array3D& get_velocity() const { return velocity; }
	/// Gets the mass
	double get_mass() const {return mass;}
	/// Gets the entire state
	/**
	 * Gets the entire spacecraft state, that is the 7 dimensional array containing
	 * position, velocity and mass. \f$ \mathbf x = [\mathbf r, \mathbf v, m]\f$
	 * \return a const reference to the 7-dimensional array containing the state
	 */
	array7D get_state() const {
	    array7D x;
	    std::copy(position.begin(), position.end(), x.begin());
	    std::copy(velocity.begin(), velocity.end(), x.begin() + 3);
	    x[6] = mass;
	    return x;
	}
	
	//@}

	/** @name Setters*/
	//@{
	/// Sets the state all at once
	/**
	 * Sets the entire spacecraft state, that is the 7 dimensional array containing
	 * position, velocity and mass.
	 *
	 * \param[in] x_ An array7D containing \f$ \mathbf x\_ = [\mathbf r, \mathbf v, m]\f$
	 */
	void set_state(const array7D& x){
	    std::copy(x.begin(), x.begin() + 3, position.begin());
	    std::copy(x.begin() + 3, x.begin() + 6, velocity.begin());
	    mass = x[6];
	}
	/// Sets the position
	void set_position(const array3D& r_){ position = r_; }
	/// Sets the velocity
	void set_velocity(const array3D& v_){ velocity = v_; }
	/// Sets the mass
	void set_mass(const double& mass_){ mass = mass_; }
	//@}
private:
#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & position;
		ar & velocity;
		ar & mass;
	}
#endif
	array3D position;
	array3D velocity;
	double mass;
};
inline std::ostream &operator<<(std::ostream &s, const sc_state &in ){
	for (int i=0;i<3;i++) s << boost::lexical_cast<std::string>(in.get_position()[i]) << " ";
	for (int i=0;i<3;i++) s << boost::lexical_cast<std::string>(in.get_velocity()[i]) << " ";
	s << boost::lexical_cast<std::string>(in.get_mass());
	return s;
}
}} //Namespaces

#endif // SC_STATE_H
