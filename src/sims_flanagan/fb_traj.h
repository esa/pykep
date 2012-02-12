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

#ifndef FB_TRAJ_H
#define FB_TRAJ_H

#include <functional>
#include <vector>

#include "leg.h"
#include "../planet.h"
#include "sc_state.h"
#include "../exceptions.h"
#include "spacecraft.h"

// Serialization code
#include "../serialization.h"
// Serialization code (END)
#include "../config.h"

namespace kep_toolbox { namespace sims_flanagan{

/// A generic interplanetary trajectory
/**
 * A multiple fly-by low-thrust trajectory represented by a sequence of Sims-Flanagan legs
 * Fly-bys are considered instantaneous and the planet::safe_radius is used to calculate the fly-by feasibility.
 * The trajectory is built by specifying a sequence of planets
 * in an std:vector container, the number of segments in each leg and a spacecraft. After that, a call to
 * init_from_full_vector specifies the specific trajectory flight-plan. Care has to be taken to respect
 * the structure of the full_vector:
 *
 * \f$ [t_0, m^i_s, \mathbf v_{\infty_s}^i, \mathbf u^i, \mathbf v_{\infty_f}^i, m^i_f, T^i] \f$.
 * The index \f$i\f$ refers to the leg.
 * \f$ t_0\f$ is the starting date in mjd2000.
 * \f$m_s\f$ is the spacecraft mass at the leg beginning (kg).
 * \f$\mathbf v_{infty_s}\f$ are the three cartesian components of the relative velocity at the leg beginning (m/s).
 * \f$\mathbf u\f$ are the 3*n_seg cartesian components of the throttles (\f$\in [0,1]\f$).
 * \f$\mathbf v_{infty_f}\f$ are the three cartesian components of the relative velocity at the leg end (m/s).
 * \f$m_f\f$ is the spacecraft mass at the leg end (kg).
 * \f$T\f$ id the leg time duration (sec).
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
class __KEP_TOOL_VISIBLE fb_traj
{
	friend std::ostream &operator<<(std::ostream &s, const fb_traj &in );
public:
	/** @name Constructors*/
	//@{
	fb_traj(const std::vector<planet_ptr>& sequence, const std::vector<int>& n_seg, const spacecraft &sc_);
	fb_traj(const std::vector<planet_ptr>& sequence, const std::vector<int>& n_seg, const double &mass_, const double &thrust_, const double &isp_);
	fb_traj(const std::vector<planet_ptr>& sequence, const unsigned int& n_seg, const spacecraft &sc_);
	fb_traj(const std::vector<planet_ptr>& sequence, const unsigned int& n_seg, const double &mass_, const double &thrust_, const double &isp_);
	fb_traj():total_n_seg(0) {}
	fb_traj(const fb_traj &);
	fb_traj &operator=(const fb_traj &);
	//@}


	/**
	 * Calculates the state mismathces at the mid-point of each leg
	 */
	template<typename it_type>
			void evaluate_all_mismatch_con(it_type begin, it_type end) const {
		assert(end - begin == 7*legs.size());
		(void) end;
		for (size_t i=0; i<legs.size();i++){
			legs[i].get_mismatch_con(begin  + 7*i, begin + 7*(i+1));
		}
	}


	/** @name Flight-plan setters*/
	//@{


	/// Sets the trajectory flight plan from a full_vector
	/**
	     * Sets the trajectory flight plan from a full_vector. The function essentially loads into the object
	     * fb_traj the unstructured information contained in a full_vector and calculates all the different
	     * constraints related to such a vector storing them in the class private members
	     *
	     * \param[in] x Full-vector encoding the trajectory flight-plan. (check full_vector for the encoding rules)
	     */
	template<typename it_type, typename coding_type>
			void init_from_full_vector(it_type b, it_type e, const coding_type& coding) {

		int n = coding.n_legs();
		assert(legs.size() == n);
		if (coding.size() != e - b) {
			throw_value_error("The provided vector size to init the trajectory is inconsistent.");
		}

		array3D start_pos, start_vel, end_pos, end_vel;
		for(int i = 0; i < n; i++){
		        planets[i]->get_eph(coding.leg_start_epoch(i, b), start_pos, start_vel);		
		        planets[i + 1]->get_eph(coding.leg_end_epoch(i, b), end_pos, end_vel);

			legs[i].set_t_i(coding.leg_start_epoch(i, b));
			array3D dv = coding.leg_start_velocity(i, b);
			std::transform(dv.begin(), dv.end(),
				       start_vel.begin(), dv.begin(),
				       std::plus<double>());
			legs[i].set_x_i(sc_state(start_pos, dv, coding.leg_start_mass(i, b)));
			legs[i].set_throttles_size(coding.n_segments(i));
			for(size_t j = 0; j < legs[i].get_throttles_size(); ++j)
				legs[i].set_throttles(j, throttle(coding.segment_start_epoch(i, j, b),
								 coding.segment_end_epoch(i, j, b),
								 coding.segment_thrust(i, j, b)));
			legs[i].set_t_f(coding.leg_end_epoch(i, b));

			dv = coding.leg_end_velocity(i, b);
			std::transform(dv.begin(), dv.end(),
				       end_vel.begin(), dv.begin(),
				       std::plus<double>());
			legs[i].set_x_f(sc_state(end_pos, dv, coding.leg_end_mass(i, b)));
		}
	}

	/// Fly-by  constraints
	/**
	 * Calculates all the constarints related to a planetary fly-by storing them in the appropriate vectors,
	 * and in particular a) the constraint on the relative velocity at the entrance and exit of the sphere of influence (fb_relvel_con), b) the constraint
	 * on the fly-by altitude of the planetocentric hyperbola (fb_altitude_con), c) the consraint on the spacecraft
	 * mass not changing during the fly-by (fb_mass_con), d) the constraint on the fly-by duration (fb_epoch_con)
	 */
	template<typename it_type>
			void evaluate_fb_con(int fb_idx, it_type begin, it_type end) {
		assert(end - begin == 2);

		array3D vin,vout,vpla;
		double emax,alfa;
		//same relvel2
		vpla = planets[fb_idx+1]->get_velocity(legs[fb_idx].get_t_f());
		vin = legs[fb_idx].get_x_f().get_velocity();
		diff(vin,vin,vpla);
		double vin2 = std::inner_product(vin.begin(), vin.end(), vin.begin(), 0);
		vout = legs[fb_idx+1].get_x_i().get_velocity();
		diff(vout,vout,vpla);
		double vout2 = std::inner_product(vout.begin(), vout.end(), vout.begin(), 0);
		begin[0] = vin2 - vout2;

		//minimum altitude (when negative this constraint is satisfied)
		emax = 1 + planets[fb_idx+1]->get_safe_radius() / planets[fb_idx+1]->get_mu_central_body()*vin2;
		alfa = acos(dot(vin,vout) / vin2);
		begin[1] = alfa - 2 * asin(1/emax);
	}

	double evaluate_leg_vinf2_i(int leg_idx) {
		array3D vout,vpla;
		vpla = planets[leg_idx]->get_velocity(legs[leg_idx].get_t_i());
		vout = legs[leg_idx].get_x_i().get_velocity();
		diff(vout,vout,vpla);
		return std::inner_product(vout.begin(), vout.end(), vout.begin(), 0.);
	}


	const leg& get_leg(int index) const {
		return legs[index];
	}

	//@}
private:
	std::vector<leg> legs;
	std::vector<planet_ptr> planets;

	//This is here only for efficiency purposes
	unsigned int total_n_seg;
	/// Contains the trajectory flight-plan
	/**
	     * According to the number of legs and the number of segments in each leg, this STD vector contains
	     * all the necessary information to actually "fly" the interplanetary trajectory. The
	     * vector structure is as follows:
	     * \f$ [t_0, m^i_s, \mathbf v_{\infty_s}^i, \mathbf u^i, \mathbf v_{\infty_f}^i, m^i_f, T^i] \f$.
	     * The index \f$i\f$ refers to the leg.
	     * \f$ t_0\f$ is the starting date in mjd2000.
	     * \f$m_s\f$ is the spacecraft mass at the leg beginning (kg).
	     * \f$\mathbf v_{infty_s}\f$ are the three cartesian components of the relative velocity at the leg beginning (m/s).
	     * \f$\mathbf u\f$ are the 3*n_seg cartesian components of the throttles (\f$\in [0,1]\f$).
	     * \f$\mathbf v_{infty_f}\f$ are the three cartesian components of the relative velocity at the leg end (m/s).
	     * \f$m_f\f$ is the spacecraft mass at the leg end (kg).
	     * \f$T\f$ id the leg time duration (sec).
	     *
	     */

// Serialization code
        friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & legs;
		ar & planets;
		ar & total_n_seg;
        }
// Serialization code (END)
};

std::ostream &operator<<(std::ostream &s, const fb_traj &in );

}} // namespaces
#endif // FB_TRAJ_H

