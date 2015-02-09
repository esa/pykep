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

#ifndef KEP_TOOLBOX_LEG_H
#define KEP_TOOLBOX_LEG_H

#include <boost/utility.hpp>
#include <boost/type_traits/is_same.hpp>
#include <iterator>
#include <vector>
#include "spacecraft.h"
#include "../core_functions/array3D_operations.h"
#include "../core_functions/propagate_lagrangian.h"
#include "../core_functions/propagate_taylor.h"
#include "sc_state.h"
#include "../epoch.h"
#include "throttle.h"
#include "../exceptions.h"

// Serialization code
#include "../serialization.h"
// Serialization code (END)
#include "../config.h"

namespace kep_toolbox {
/// Sims-Flanagan transcription of low-thrust trajectories
/**
* This namespace contains the routines that allow building and evaluating low-thrust trajectories using the
* Sims-Flanagan transcription method.
*/
namespace sims_flanagan{

/// Single low-thrust leg (phase)
/**
* This class represents, generically, a low-thrust leg (phase) as a sequence of successive
* impulses of magnitude compatible with the low-thrust propulsion system of a spacecraft.
* The leg achieves to transfer a given spacecraft from an initial to a final state in the
* time given (and can be considered as feasible) whenever the method evaluate_mismatch
* returns all zeros and the method get_throttles_con returns all values less than zero.
* Th sequence of different impulses is represented by the class throttles. These represent
* the cartesian components \f$ \mathbf x = (x_1,y_1,z_1) \f$ of a normalized \f$ \Delta V \f$ and are thus
* numbers that need to satisfy the constraint \f$|\mathbf x| \le 1\f$
*
* \image html sims_flanagan_leg.png "Visualization of a feasible leg (Earth-Mars)"
* \image latex sims_flanagan_leg.png "Visualization of a feasible leg (Earth-Mars)" width=5cm
*
* @author Dario Izzo (dario.izzo _AT_ googlemail.com)
*/
class __KEP_TOOL_VISIBLE leg
{
	friend std::ostream &operator<<(std::ostream &s, const leg &in );

public:
	std::string human_readable() const;
	/// Constructor.
	/**
	* Default constructor. Constructs a meaningless leg that will need to be properly initialized
	* using the various setters....
	*/
	leg():t_i(),x_i(),throttles(),t_f(),x_f(),m_sc(),m_mu(0),m_hf(false),m_tol(-10) {}

	/// Constructs the leg from epochs, sc_states and cartesian components of throttles
	/**
	* Constructs entirely a leg assuming high fidelity propagation switched off and equally spaced segments
	*
	* \param[in] epoch_i Inital epoch
	* \param[in] state_i Initial sc_state (spacecraft state)
	* \param[in] thrott sequence of doubles (\f$ x_1,y_1,z_1, ..., x_N,y_N,z_N \f$) representing the cartesian components of the throttles
	* \param[in] state_f Final sc_state (spacecraft state)
	* \param[in] sc Spacecraft object
	* \param[in] mu Primary body gravitational constant
	*
	*/
	leg(const epoch& epoch_i, const sc_state& state_i, const std::vector<double>& thrott,
	    const epoch& epoch_f, const sc_state& state_f, const spacecraft& sc, const double mu):m_sc(sc),m_hf(false),m_tol(-10) {
		set_leg(epoch_i, state_i,thrott.begin(),thrott.end(),epoch_f, state_f,mu);
	}

	/// Initialize a leg
	/**
	* Initialize a leg assuming that the user has or will initialize separately the spacecraft.
	* The throttles are provided via two iterators pointing
	* to the beginning and to the end of a throttle sequence.
	*
	* \param[in] epoch_i Inital epoch
	* \param[in] state_i Initial sc_state (spacecraft state)
	* \param[in] throttles_start iterator pointing to the beginning of a cartesian throttle sequence.
	* \param[in] throttles_end iterator pointing to the end+1 of a cartesian throttle sequence.
	* \param[in] epoch_f Final epoch. Needs to be later than epoch_i
	* \param[in] state_f Final sc_state (spacecraft state)
	* \param[in] mu_ Primary body gravitational constant
	*
	& \throws value_error if final epoch is before initial epoch, if mu_ not positive
	*/
	template<typename it_type>
	void set_leg(const epoch& epoch_i, const sc_state& state_i,
		     it_type throttles_start,
		     it_type throttles_end,
		     const epoch& epoch_f, const sc_state& state_f,
		     double mu_, typename boost::enable_if<boost::is_same<typename std::iterator_traits<it_type>::value_type,throttle> >::type * = 0)
	{
		if (epoch_f.mjd2000() <= epoch_i.mjd2000())
		{
			throw_value_error("Final epoch is before initial epoch");
		}

		t_i = epoch_i; x_i=state_i;
		t_f = epoch_f; x_f=state_f;

		throttles.assign(throttles_start, throttles_end);

		if (mu_<=0)
		{
			throw_value_error("Gravitational constant is less or equal to zero");
		}
		m_mu = mu_;
	}
	
	/// Initialize a leg
	/**
	* Initialize a leg assuming that the user has or will initialize separately the spacecraft.
	* The throttles are provided via two iterators pointing to the beginning and to the end of
	* a sequence of doubles (\f$ x_1,y_1,z_1, ..., x_N,y_N,z_N \f$ containing the
	* cartesian components of each throttle \f$ x_i,y_i,z_i \in [0,1]\f$. The constructed leg
	* will have by default equally spaced segments.
	*
	* \param[in] epoch_i Inital epoch
	* \param[in] state_i Initial sc_state (spacecraft state)
	* \param[in] throttles_start iterator pointing to the beginning of a cartesian throttle sequence.
	* \param[in] throttles_end iterator pointing to the end+1 of a cartesian throttle sequence.
	* \param[in] epoch_f Final epoch. Needs to be later than epoch_i
	* \param[in] state_f Final sc_state (spacecraft state)
	* \param[in] mu_ Primary body gravitational constant
	*
	& \throws value_error if final epoch is before initial epoch, if mu_ not positive
	*/
	
	template<typename it_type>
	void set_leg(const epoch& epoch_i, const sc_state& state_i,
		it_type throttles_start,
		it_type throttles_end,
		const epoch& epoch_f, const sc_state& state_f,
		double mu_, typename boost::enable_if<boost::is_same<typename std::iterator_traits<it_type>::value_type,double> >::type * = 0)
	{
		// 		if (epoch_f.mjd2000() <= epoch_i.mjd2000())
		// 		{
		// 			throw_value_error("Final epoch is before initial epoch");
		// 		}
		if (std::distance(throttles_start,throttles_end) % 3 || std::distance(throttles_start,throttles_end) <= 0) {
			throw_value_error("The length of the throttles list must be positive and a multiple of 3");
		}
		if (mu_<=0)
		{
			throw_value_error("Gravitational constant is less or equal to zero");
		}
		m_mu = mu_;

		t_i = epoch_i;
		x_i=state_i;
		t_f = epoch_f;
		x_f=state_f;

		const int throttles_vector_size = (throttles_end - throttles_start) / 3;
		throttles.resize(throttles_vector_size);
		const double seg_duration = (epoch_f.mjd() - epoch_i.mjd()) / throttles_vector_size;
		for (int i = 0; i < throttles_vector_size; ++i) {
			kep_toolbox::array3D tmp = {{ *(throttles_start + 3 * i), *(throttles_start + 3 * i + 1), *(throttles_start + 3 * i + 2) }};
			throttles[i] = throttle(epoch(epoch_i.mjd() + seg_duration * i,epoch::MJD),epoch(epoch_i.mjd() + seg_duration * (i + 1),epoch::MJD),tmp);
		}
	}

	/// Initialize a leg
	void set_leg(const epoch& epoch_i, const sc_state& state_i,const std::vector<double>& thrott,
		     const epoch& epoch_f, const sc_state& state_f) {
		set_leg(epoch_i, state_i,thrott.begin(),thrott.end(),epoch_f, state_f,m_mu);
	}

	/** @name Setters*/
	//@{

	/// Sets the leg's spacecraft
	/**
	*
	*In order for the trajectory leg to be able to propagate the states, information on the
	* low-thrust propulsion system used needs to be available. This is provided by the object
	* spacecraft private member of the class and can be set using this setter.
	*
	* \param[in] sc The spacecraft object
	*/
	void set_spacecraft(const spacecraft &sc) { m_sc = sc; }

	/// Sets the leg's primary body gravitational parameter
	/**
	*
	* Sets the leg's central body gravitational parameter
	*
	* \param[in] mu_ The gravitational parameter
	*/
	void set_mu(const double &mu_) { m_mu = mu_; }

	/// Sets the throttles
	/**
	*
	* \param[in] b iterator pointing to the begin of a throttles sequence
	* \param[in] e iterator pointing to the end of a throttles sequence
	*/
	template<typename it_type>
	void set_throttles(it_type b, it_type e) { throttles.assign(b, e); }

	/// Sets the throttles size
	/**
	* Resizes the throttles vector to a new size.
	*
	* \param[in] size The new size of the throttles vector (number of segments)
	*/
	void set_throttles_size(const int& size) { throttles.resize(size); }

	/**
	* Sets the throttles
	*
	* \param[in] t std::vector of throttles
	*/
	void set_throttles(const std::vector<throttle>& t) { throttles = t; }
	
	/**
	* Sets the ith throttle
	*
	* \param[in] index the index of the throttle
	* \param[in] t the throttle
	*/
	void set_throttles(int index, const throttle& t) { throttles[index] = t; }

	/// Sets the final sc_state
	/**
	* Sets the spacecraft state at the end of the leg
	*
	*/
	void set_x_f(const sc_state& s) { x_f = s;}

	/// Sets the initial sc_state
	/**
	* Sets the spacecraft state at the beginning of the leg
	*
	*/
	void set_x_i(const sc_state& s) { x_i = s;}
	//@}

	/// Sets the initial epoch
	/**
	* Sets the epoch at the beginning of the leg
	*
	*/
	void set_t_i(epoch e) { t_i = e; }

	/// Sets the final epoch
	/**
	* Sets the epoch at the end of the leg
	*
	*/
	void set_t_f(epoch e) { t_f = e; }

	/// Sets the leg high difelity state
	/**
	* Activates the evaluation of the state-mismatches using a high-fidelity model. Resulting leg
	* is a real low-thrust trajectory (i.e. no impulses approximation)
	*
	*/
	void set_high_fidelity(bool state) { m_hf = state; }


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

	/// Gets the throttle vector size
	/**
	* Returns the throttle vector size (number of segments)
	*
	* @return size_t containing the throttle vector size.
	*/
	size_t get_throttles_size() const {return throttles.size();}

	/// Gets the i-th throttle
	/**
	* Returns the i-th throttle
	*
	* @return const ref to the i-th throttle
	*/	
	const throttle& get_throttles(int index) { return throttles[index]; }
	
	/// Gets the throttles
	/**
	* Returns all throttles
	*
	* @return const ref to a vector of throttle
	*/	
	const std::vector<throttle>& get_throttles() { return throttles; }

	/// Gets the leg's initial epoch
	/**
	* Gets the epoch at the beginning of the leg
	*
	* @return const reference to the initial epoch
	*/
	const epoch& get_t_i() const {return t_i;}

	/// Gets the leg's final epoch
	/**
	* Gets the epoch at the end of the leg
	*
	* @return const reference to the final epoch
	*/
	const epoch& get_t_f() const {return t_f;}

	/// Gets the sc_state at the end of the leg
	/**
	* Gets the spacecraft state at the end of the leg
	*
	* @return const reference to the final sc_state
	*/
	const sc_state& get_x_f() const {return x_f;}

	/// Gets the initial sc_state
	/**
	* Gets the spacecraft state at the beginning of the leg
	*
	* @return const reference to the initial sc_state
	*/
	const sc_state& get_x_i() const {return x_i;}
	bool get_high_fidelity() const { return m_hf; }
	//@}

	/** @name Leg Feasibility*/
	//@{

	/// Evaluate the state mismatch
	/**
	* This is the main method of the class leg as it performs the orbital propagation from the initial sc_state, and
	* accounting for all the throttles, up to a mid-point. The same is done starting from the final sc_state up to
	* the same mid-point. The difference between the obtained values is then recorded at the memory location pointed by the iterators
	* If not all zero the leg is unfeasible. The values stored are \f$\mathbf r, \mathbf v, m\f$
	*
	* @param[in] begin iterator pointing to the beginning of the memory where the mismatches will be stored
	* @param[in] begin iterator pointing to the end of the memory where the mismatches will be stored
	*/

	template<typename it_type>
	void get_mismatch_con(it_type begin, it_type end) const
	{
		if (m_hf) {
			get_mismatch_con_low_thrust(begin, end);
		} else {
			get_mismatch_con_chemical(begin, end);
		}
	}

protected:
	template<typename it_type>
	void get_mismatch_con_chemical(it_type begin, it_type end) const
	{
		assert(end - begin == 7);
		(void)end;
		size_t n_seg = throttles.size();
		const int n_seg_fwd = (n_seg + 1) / 2, n_seg_back = n_seg / 2;

		//Aux variables
		double max_thrust = m_sc.get_thrust();
		double isp = m_sc.get_isp();
		double norm_dv;
		array3D dv;

		//Initial state
		array3D rfwd = x_i.get_position();
		array3D vfwd = x_i.get_velocity();
		double mfwd = x_i.get_mass();

		//Forward Propagation
		double current_time_fwd = t_i.mjd2000() * ASTRO_DAY2SEC;
		for (int i = 0; i < n_seg_fwd; i++) {
			double thrust_duration = (throttles[i].get_end().mjd2000() -
						  throttles[i].get_start().mjd2000()) * ASTRO_DAY2SEC;
			double manouver_time = (throttles[i].get_start().mjd2000() +
						throttles[i].get_end().mjd2000()) / 2. * ASTRO_DAY2SEC;
			propagate_lagrangian(rfwd, vfwd, manouver_time - current_time_fwd, m_mu);
			current_time_fwd = manouver_time;

			for (int j=0;j<3;j++){
				dv[j] = max_thrust / mfwd * thrust_duration * throttles[i].get_value()[j];
			}

			norm_dv = norm(dv);
			sum(vfwd,vfwd,dv);
			mfwd *= exp( -norm_dv/isp/ASTRO_G0 );
			//Temporary solution to the creation of NaNs when mass gets too small (i.e. 0)
			if (mfwd < 1) mfwd=1;
		}

		//Final state
		array3D rback = x_f.get_position();
		array3D vback = x_f.get_velocity();
		double mback = x_f.get_mass();

		//Backward Propagation
		double current_time_back = t_f.mjd2000() * ASTRO_DAY2SEC;
		for (int i = 0; i < n_seg_back; i++) {
			double thrust_duration = (throttles[throttles.size() - i - 1].get_end().mjd2000() -
						  throttles[throttles.size() - i - 1].get_start().mjd2000()) * ASTRO_DAY2SEC;
			double manouver_time = (throttles[throttles.size() - i - 1].get_start().mjd2000() +
						throttles[throttles.size() - i - 1].get_end().mjd2000()) / 2. * ASTRO_DAY2SEC;
			// manouver_time - current_time_back is negative, so this should propagate backwards
			propagate_lagrangian(rback, vback, manouver_time - current_time_back, m_mu);
			current_time_back = manouver_time;

			for (int j=0;j<3;j++){
				dv[j] = - max_thrust / mback * thrust_duration * throttles[throttles.size() - i - 1].get_value()[j];
			}
			norm_dv = norm(dv);
			sum(vback,vback,dv);
			mback *= exp( norm_dv/isp/ASTRO_G0 );
		}

		// finally, we propagate from current_time_fwd to current_time_back with a keplerian motion
		propagate_lagrangian(rfwd, vfwd, current_time_back - current_time_fwd, m_mu);

		//Return the mismatch
		diff(rfwd,rfwd,rback);
		diff(vfwd,vfwd,vback);

		std::copy(rfwd.begin(), rfwd.end(), begin);
		std::copy(vfwd.begin(), vfwd.end(), begin + 3);
		begin[6] = mfwd - mback;
	}


	template<typename it_type>
	void get_mismatch_con_low_thrust(it_type begin, it_type end) const
	{
		assert(end - begin == 7);
		(void)end;
		size_t n_seg = throttles.size();
		const int n_seg_fwd = (n_seg + 1) / 2, n_seg_back = n_seg / 2;

		//Aux variables
		double max_thrust = m_sc.get_thrust();
		double veff = m_sc.get_isp()*ASTRO_G0;
		array3D thrust;

		//Initial state
		array3D rfwd = x_i.get_position();
		array3D vfwd = x_i.get_velocity();
		double mfwd = x_i.get_mass();

		//Forward Propagation
		for (int i = 0; i < n_seg_fwd; i++) {
			double thrust_duration = (throttles[i].get_end().mjd2000() -
						  throttles[i].get_start().mjd2000()) * ASTRO_DAY2SEC;

			for (int j=0;j<3;j++){
				thrust[j] = max_thrust * throttles[i].get_value()[j];
			}
			propagate_taylor(rfwd,vfwd,mfwd,thrust,thrust_duration,m_mu,veff,m_tol,m_tol);
		}

		//Final state
		array3D rback = x_f.get_position();
		array3D vback = x_f.get_velocity();
		double mback = x_f.get_mass();

		//Backward Propagation
		for (int i = 0; i < n_seg_back; i++) {
			double thrust_duration = (throttles[throttles.size() - i - 1].get_end().mjd2000() -
						  throttles[throttles.size() - i - 1].get_start().mjd2000()) * ASTRO_DAY2SEC;
			for (int j=0;j<3;j++){
				thrust[j] = max_thrust * throttles[throttles.size() - i - 1].get_value()[j];
			}
			propagate_taylor(rback,vback,mback,thrust,-thrust_duration,m_mu,veff,m_tol,m_tol);
		}

		//Return the mismatch
		diff(rfwd,rfwd,rback);
		diff(vfwd,vfwd,vback);

		std::copy(rfwd.begin(), rfwd.end(), begin);
		std::copy(vfwd.begin(), vfwd.end(), begin + 3);
		begin[6] = mfwd - mback;
	}



public:
	/// Evaluate the state mismatch
	/**
	* This method overloads the same method using iterators but the mismatch values are stored
	* in a sc_state object
	*
	* @param[in] retval the state mismatch structured as a spacecraft state
	*/
	void get_mismatch_con(sc_state& retval) const
	{
		array7D tmp;
		get_mismatch_con(tmp.begin(), tmp.end());
		retval.set_state(tmp);
	}

	/// Evaluate the throttles magnitude
	/**
	* This methods loops on the vector containing the throttles \f$ (x_1,y_1,z_1,x_2,y_2,z_2,...,x_n,y_n,z_n) \f$
	* and stores the magnitudes \f$ x_i^2 + y_i^2 + z_i^2 - 1\f$ at the locations indicated by the iterators. The
	* iterators must have a distance of \f$ n\f$. If the stored values are not all \f$ \le 0 \f$ then the trajectory
	* is unfeasible.
	*
	* @param[out] start std::vector<double>iterator from the first element where to store the magnitudes
	* @param[out] start std::vector<double>iterator to the last+1 element where to store the magnitudes
	*/
	template<typename it_type>
	void get_throttles_con(it_type start, it_type end) const {
		if ( (end - start) != (int)throttles.size()) {
			throw_value_error("Iterators distance is incompatible with the throttles size");
		}
		int i=0;
		while(start!=end){
			const array3D& t = throttles[i].get_value();
			*start = std::inner_product(t.begin(), t.end(), t.begin(), -1.);
			++i; ++start;
		}
	}
	//@}

	/// Approximate the leg dv
	/**
	* This method returns the leg dv assuming a constant mass throughout the leg. Useful when
	* mass propagation is switched off by setting an infinite specific impulse
	*
	* @return the leg dv (mass is considered constant)
	*/
	double evaluate_dv() const
	{
		double tmp = 0;
		for (std::vector<double>::size_type i = 0; i < throttles.size(); ++i)
		{
			tmp += (throttles[i].get_end().mjd2000() -throttles[i].get_start().mjd2000())
					* ASTRO_DAY2SEC * throttles[i].get_norm() * m_sc.get_thrust() / m_sc.get_mass();
		}
		return tmp;
	}

private:
// Serialization code
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & t_i;
			ar & x_i;
			ar & throttles;
			ar & t_f;
			ar & x_f;
			ar & m_sc;
			ar & m_mu;
			ar & m_hf;
			ar & m_tol;
		}
// Serialization code (END)
		epoch t_i;
		sc_state x_i;
		std::vector<throttle> throttles;
		epoch t_f;
		sc_state x_f;
		spacecraft m_sc;
		double m_mu;
		bool m_hf;
		int m_tol;
	};

std::ostream &operator<<(std::ostream &s, const leg &in );

}} //namespaces
#endif // KEP_TOOLBOX_LEG_H
