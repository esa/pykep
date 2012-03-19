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

#ifndef LEG_S_H
#define LEG_S_H

#include <boost/utility.hpp>
#include <boost/type_traits/is_same.hpp>
#include <iterator>
#include <vector>

#include "spacecraft.h"
#include "../core_functions/array3D_operations.h"
#include "../core_functions/propagate_taylor_s.h"
#include "sc_state.h"
#include "../epoch.h"
#include "throttle.h"
#include "../exceptions.h"

// Serialization code
#include "../serialization.h"
// Serialization code (END)
#include "../config.h"

namespace kep_toolbox {
namespace sims_flanagan{

/// Single low-thrust leg (phase) using Sundmann variable
/**
* This class represents, generically, a low-thrust leg (phase) as a sequence of successive
* constant low-thrust segments. The segment duration is equally distributed in the pseudo
* time space \f$ds = c r^\alpha dt \f$
* The leg achieves to transfer a given spacecraft from an initial to a final state in the
* pseudo-time given (and can thus be considered as feasible) whenever the method evaluate_mismatch
* returns all zeros (8 values) and the method get_throttles_con returns all values less than zero.
* Th sequence of different thrusts is represented by the class throttles. These represent
* the cartesian components \f$ \mathbf u = (u_x,u_y,u_y) \f$ of a normalized thrust and are thus
* numbers that need to satisfy the constraint \f$|\mathbf u| \le 1\f$
*
* \image html s_leg.png "Visualization of a feasible leg (Earth-Mars)"
* \image latex s_leg.png "Visualization of a feasible leg (Earth-Mars)" width=5cm
*
* @author Dario Izzo (dario.izzo _AT_ googlemail.com)
*/
class __KEP_TOOL_VISIBLE leg_s
{
	friend std::ostream &operator<<(std::ostream &s, const leg &in );

public:
	std::string human_readable() const;
	/// Constructor.
	/**
	* Default constructor. Constructs a meaningless leg.
	*/
	leg_s(): m_ti(), m_xi(), m_throttles(), m_tf(), m_xf(), m_sf(0), m_sc(), m_mu(0), m_tol(-10), m_states(), m_ceq(8), m_cineq(), m_dv() {}

	/// Constructor
	/**
	 * Constructs an empty leg allocating memory for a given number of segments.
	*/
	leg_s(const unsigned int& n_seg, const double& tol=1e-10): m_ti(), m_xi(), m_throttles(n_seg), m_tf(), m_xf(), m_sf(0), m_sc(), m_mu(), m_tol(tol), m_states(n_seg+2), m_ceq(8), m_cineq(n_seg), m_dv(n_seg) {}

	/// Sets the leg's data
	/**
	* The throttles are provided via two iterators pointing to the beginning and to the end of
	* a sequence of doubles (\f$ x_1,y_1,z_1, ..., x_N,y_N,z_N \f$ containing the
	* cartesian components of each throttle \f$ x_i,y_i,z_i \in [0,1]\f$. The constructed leg
	* will have equally spaced segments in the pseudo-time
	*
	* \param[in] epoch_i Inital epoch
	* \param[in] state_i Initial sc_state (spacecraft state)
	* \param[in] throttles_start iterator pointing to the beginning of a cartesian throttle sequence.
	* \param[in] throttles_end iterator pointing to the end+1 of a cartesian throttle sequence.
	* \param[in] epoch_f Final epoch. Needs to be later than epoch_i
	* \param[in] state_f Final sc_state (spacecraft state)
	* \param[in] s_f Pseudo time of transfer
	* \param[in] mu_ Primary body gravitational constant
	* \param[in] sc_ Spacecraft
	*
	& \throws value_error if final epoch is before initial epoch, if mu_ not positive if the throttle size is not consistent
	*/
	template<typename it_type>
	void set_leg(const epoch& epoch_i, const sc_state& state_i,
		it_type throttles_start,
		it_type throttles_end,
		const epoch& epoch_f,
		const sc_state& state_f, const double& sf,
		const spacecraft &sc_, const double &mu_, typename boost::enable_if<boost::is_same<typename std::iterator_traits<it_type>::value_type,double> >::type * = 0)
	{
		//We check data consistency
		if (std::distance(throttles_start,throttles_end) % 3) {
			throw_value_error("The length of the throttles list must be a multiple of 3");
		}
		if (std::distance(throttles_start,throttles_end) / 3 != (int)m_throttles.size()) {
			throw_value_error("The number of segments in the leg do not match the length of the supplied throttle sequence");
		}

		if (mu_<=0)
		{
			throw_value_error("Gravity parameter must be larger than zero");
		}
		if (epoch_i.mjd() >= epoch_f.mjd()) {
			throw_value_error("Final epoch must be after the initial epoch");
		}

		//We fill up all leg's data member
		m_mu = mu_;
		m_sc = sc_;
		m_ti = epoch_i;
		m_xi = state_i;
		m_tf = epoch_f;
		m_xf = state_f;
		m_sf = sf;

		//note: the epochs of the throttles are meaningless at this point as pseudo-time is used
		for (size_t i = 0; i < m_throttles.size(); ++i) {
			kep_toolbox::array3D tmp = {{ *(throttles_start + 3 * i), *(throttles_start + 3 * i + 1), *(throttles_start + 3 * i + 2) }};
			m_throttles[i] = throttle(epoch(i),epoch(i+1),tmp);
		}
	}

	/// Sets the leg's data
	void set_leg(const epoch& epoch_i, const sc_state& state_i,const std::vector<double>& thrott,
		     const epoch& epoch_f, const sc_state& state_f, const double& sf) {
		set_leg(epoch_i, state_i,thrott.begin(),thrott.end(),epoch_f, state_f,sf, m_sc, m_mu);
	}

	/// Sets the leg's data
	void set_leg(const epoch& epoch_i, const sc_state& state_i,const std::vector<double>& thrott,
		     const epoch& epoch_f, const sc_state& state_f, const double& sf, const spacecraft &sc_, const double &mu_) {
		set_leg(epoch_i, state_i,thrott.begin(),thrott.end(),epoch_f, state_f,sf, sc_, mu_);
	}

	/** @name Setters*/
	//@{

	/// Sets the leg's spacecraft
	/**
	*
	* In order for the trajectory leg to be able to propagate the states, information on the
	* low-thrust propulsion system used needs to be available. This is provided by the object
	* spacecraft private member of the class and can be set using this setter.
	*
	* \param[in] sc The spacecraft object
	*/
	void set_sc(const spacecraft &sc) { m_sc = sc; }

	/// Sets the leg's primary body gravitational parameter
	/**
	*
	* Sets the leg's central body gravitational parameter
	*
	* \param[in] mu_ The gravitational parameter
	*/
	void set_mu(const double &mu_) { m_mu = mu_; }

	/** @name Getters*/
	//@{

	/// Gets the leg's spacecraft
	/**
	* Returns the spacecraft
	*
	* @return sc const reference to spacecraft object
	*/
	const spacecraft& get_spacecraft() const { return m_sc; }
	
	/// Gets the gravitational parameter
	/**
	* @return the gravitational parameter
	*/
	double get_mu() const { return m_mu; }

	/// Gets the leg' s number of segments
	/**
	* Returns the leg' s number of segments
	*
	* @return size_t the leg's number of segments
	*/
	size_t get_n_seg() const {return m_throttles.size(); }

	/// Gets the i-th throttle
	/**
	* Returns the i-th throttle
	*
	* @return const ref to the i-th throttle
	*/	
	const throttle& get_throttles(int index) { return m_throttles[index]; }
	
	/// Gets the throttles
	/**
	* Returns all throttles
	*
	* @return const ref to a vector of throttle
	*/	
	const std::vector<throttle>& get_throttles() { return m_throttles; }

	/// Gets the leg's initial epoch
	/**
	* Gets the epoch at the beginning of the leg
	*
	* @return const reference to the initial epoch
	*/
	const epoch& get_ti() const {return m_ti;}

	/// Gets the leg's final epoch
	/**
	* Gets the epoch at the end of the leg
	*
	* @return const reference to the final epoch
	*/
	const epoch& get_tf() const {return m_tf;}

	/// Gets the sc_state at the end of the leg
	/**
	* Gets the spacecraft state at the end of the leg
	*
	* @return const reference to the final sc_state
	*/
	const sc_state& get_xf() const {return m_xf;}

	/// Gets the initial sc_state
	/**
	* Gets the spacecraft state at the beginning of the leg
	*
	* @return const reference to the initial sc_state
	*/
	const sc_state& get_xi() const {return m_xi;}
	//@}

	/** @name Computations*/
	//@{

	/// Returns the computed state mismatch constraints (8 equality constraints)
	/**
	*/
	const std::vector<double>&  compute_mismatch_con() const {
		size_t n_seg = m_throttles.size();
		const int n_seg_fwd = (n_seg + 1) / 2, n_seg_back = n_seg / 2;

		//Aux variables
		double max_thrust = m_sc.get_thrust();
		double veff = m_sc.get_isp()*ASTRO_G0;
		array3D thrust;
		double ds = m_sf/n_seg; //pseudo-time interval for each segment

		//Initial state
		array3D rfwd = m_xi.get_position();
		array3D vfwd = m_xi.get_velocity();
		double mfwd = m_xi.get_mass();
		double tfwd = 0;

		//Forward Propagation
		for (int i = 0; i < n_seg_fwd; i++) {
			for (int j=0;j<3;j++){
				thrust[j] = max_thrust * m_throttles[i].get_value()[j];
			}
			propagate_taylor_s(rfwd,vfwd,mfwd,tfwd,thrust,ds,m_mu,veff,1.0,1.5,m_tol,m_tol);
		}

		//Final state
		array3D rback = m_xf.get_position();
		array3D vback = m_xf.get_velocity();
		double mback = m_xf.get_mass();
		double tback = 0;

		//Backward Propagation
		for (int i = 0; i < n_seg_back; i++) {
			for (int j=0;j<3;j++){
				thrust[j] = max_thrust * m_throttles[m_throttles.size() - i - 1].get_value()[j];
			}
			propagate_taylor_s(rback,vback,mback,tback,thrust,-ds,m_mu,veff,1.0,1.5,m_tol,m_tol);
		}

		//Return the mismatch
		diff(rfwd,rfwd,rback);
		diff(vfwd,vfwd,vback);

		std::copy(rfwd.begin(), rfwd.end(), m_ceq.begin());
		std::copy(vfwd.begin(), vfwd.end(), m_ceq.begin() + 3);
		m_ceq[6] = mfwd - mback;
		m_ceq[7] = tfwd + tback;
		return m_ceq;
	}
	/// Returns the computed throttles constraints (n_seg inequality constraints)
	/**
	*/
	const std::vector<double>&  compute_throttles_con() const {
		for (size_t i=0; i<m_throttles.size();++i){
			const array3D& t = m_throttles[i].get_value();
			m_cineq[i] = std::inner_product(t.begin(), t.end(), t.begin(), -1.);
		}
		return m_cineq;
	}
	//const std::vector<double>&  compute_dvs() const { return m_dv; }
	//const std::vector<boost::array<double,8> >& compute_states() const { return m_states; }

protected:
	void record_states(const double& t, const array3D& r, const array3D& v,const double& m,const unsigned int& idx) {
		assert(idx < m_states.size());
		m_states[idx][0] = t;
		std::copy(r.begin(),r.end(),m_states[idx].begin());
		std::copy(v.begin(),v.end(),m_states[idx].begin()+3);
		m_states[idx][7] = m;
	}
	//@}

private:
	// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & m_ti;
		ar & m_xi;
		ar & m_throttles;
		ar & m_tf;
		ar & m_xf;
		ar & m_sc;
		ar & m_mu;
		ar & m_tol;
		//ar & m_states;
		ar & m_ceq;
		ar & m_cineq;
		//ar & m_dv;
	}
	// Serialization code (END)
	epoch m_ti;
	sc_state m_xi;
	std::vector<throttle> m_throttles;
	epoch m_tf;
	sc_state m_xf;
	double m_sf;
	spacecraft m_sc;
	double m_mu;
	int m_tol;

	mutable std::vector<boost::array<double, 8> > m_states;
	mutable std::vector<double> m_ceq;
	mutable std::vector<double> m_cineq;
	mutable std::vector<double> m_dv;
};

std::ostream &operator<<(std::ostream &s, const leg &in );

}} //namespaces
#endif // LEG_H
