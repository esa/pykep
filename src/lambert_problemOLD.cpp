/*****************************************************************************
 *   Copyright (C) 2004-2012 The PyKEP development team,                     *
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

#include "lambert_problemOLD.h"
#include "core_functions/array3D_operations.h"
#include "core_functions/lambert_find_N.h"
#include "core_functions/lambert_3d.h"

namespace kep_toolbox {

const array3D lambert_problemOLD::default_r1 = {{1,0,0}};
const array3D lambert_problemOLD::default_r2 = {{0,1,0}};

/// Constructor
/** It constructs and solves a Lambert problem.
 *
 * \param[in] r1 first cartesian position
 * \param[in] r2 second cartesian position
 * \param[in] tof time of flight
 * \param[in] mu gravity parameter
 * \param[in] cw when 1 a retrograde orbit is assumed
 * \param[in] multi_revs when true computes also all multiple revolutions soluitons
 */
lambert_problemOLD::lambert_problemOLD(const array3D &r1, const array3D &r2, const double &tof, const double& mu, const int &cw, const int &multi_revs) :
				m_r1(r1), m_r2(r2),m_tof(tof),m_mu(mu),m_has_converged(true), m_multi_revs(multi_revs)
{
	// 1 - Computing non dimensional units
	double R = norm(r1);
	double V = sqrt(mu / R);
	double T = R/V;

	// 2 - Computing geometry of transfer in non_dimensional units
	double R2 = norm(r2);
	double costheta = dot(r1,r2);
	costheta /= R*R2;
	double r2_mod = R2 / R;
	m_c = sqrt(1 + r2_mod*(r2_mod - 2.0 * costheta));
	m_s = (1 + r2_mod + m_c)/2.0;
	// 2a - long or short way?
	m_lw = ( (r1[0]*r2[1] - r1[1]*r2[0]) > 0 ) ? 0 : 1;	//prograde motion assumed
	if (cw) m_lw = (m_lw+1) % 2;				//changed to retrograde motion

	// 3 - computing maximum number of revolutions
	m_Nmax = 0;
	if (m_multi_revs>0) {
		m_Nmax = lambert_find_N(m_s,m_c,tof/T,m_lw);
	}
	m_Nmax = std::min(m_multi_revs,m_Nmax);
	// 4 - computing all solutions
	m_v1.resize(m_Nmax * 2 +1);
	m_v2.resize(m_Nmax * 2 +1);
	m_iters.resize(m_Nmax * 2 +1);
	m_a.resize(m_Nmax * 2 +1);
	m_p.resize(m_Nmax * 2 +1);
	//no rev solution
	m_iters[0] = lambert_3d(m_v1[0],m_v2[0],m_a[0],m_p[0],r1,r2,tof,mu,m_lw);
	//multirev solution
	for (int i=0;i<m_Nmax;++i)
	{
		m_iters[1+2*i] = lambert_3d(m_v1[1+2*i],m_v2[1+2*i],m_a[1+2*i],m_p[1+2*i],r1,r2,tof,mu,m_lw,i+1,'l');
		m_iters[2+2*i] = lambert_3d(m_v1[2+2*i],m_v2[2+2*i],m_a[2+2*i],m_p[2+2*i],r1,r2,tof,mu,m_lw,i+1,'r');
	}
	for (std::vector<int>::size_type i=0;i<m_iters.size();++i){
		if (m_iters[i] == ASTRO_MAX_ITER) m_has_converged = false;
	}
}

/// Reliability check
/** Checks that all lambert solver calls have terminated within the maximum allowed iteration
 * indicating convergence of the tof curve solution
 *
 * \return true if all solutions have converged
 */
bool lambert_problemOLD::is_reliable() const
{
	return m_has_converged;
}

/// Gets velocity at r1
/**
 *
 * \return an std::vector containing 3-d arrays with the cartesian components of the velocities at r1 for all 2N_max+1 solutions
 */
const std::vector<array3D>& lambert_problemOLD::get_v1() const
{
	return m_v1;
}

/// Gets velocity at r2
/**
 *
 * \return an std::vector containing 3-d arrays with the cartesian components of the velocities at r2 for all 2N_max+1 solutions
 */
const  std::vector<array3D>& lambert_problemOLD::get_v2() const
{
	return m_v2;
}

/// Gets r1
/**
 *
 * \return a 3-d array with the cartesian components of r1
 */
const array3D& lambert_problemOLD::get_r1() const
{
	return m_r1;
}

/// Gets r2
/**
 *
 * \return a 3-d array with the cartesian components of r2
 */
const array3D& lambert_problemOLD::get_r2() const
{
	return m_r2;
}

/// Gets the time of flight between r1 and r2
/**
 *
 * \return the time of flight
 */
const double& lambert_problemOLD::get_tof() const
{
	return m_tof;
}

/// Gets gravitational parameter
/**
 *
 * \return the gravitational parameter
 */
const double& lambert_problemOLD::get_mu() const
{
	return m_mu;
}


/// Gets semi-major axes
/**
 *
 * \return an std::vector containing the semi-major axes all 2N_max+1 solutions
 */
const  std::vector<double>& lambert_problemOLD::get_a() const
{
	return m_a;
}
/// Gets parameters
/**
 *
 * \return an std::vector containing the parameters of all 2N_max+1 solutions
 */
const  std::vector<double>& lambert_problemOLD::get_p() const
{
	return m_p;
}

/// Gets number of iterations
/**
 *
 * \return an std::vector containing the iterations taken to compute each one of the solutions
 */
const  std::vector<int>& lambert_problemOLD::get_iters() const
{
	return m_iters;
}

/// Gets N_max
/**
 *
 * \return the maximum number of revolutions. The number of solutions to the problem will be Nmax*2 +1
 */
int lambert_problemOLD::get_Nmax() const
{
	return m_Nmax;
}

/// Streaming operator
std::ostream &operator<<(std::ostream &s, const lambert_problemOLD &lp) {
	s << "Lambert's problem:" << std::endl;
	s << "r1 = " << lp.m_r1 << std::endl;
	s << "r2 = " << lp.m_r2 << std::endl;
	s << "Time of flight: " << lp.m_tof <<std::endl;
	s << "Gravity paramter: " << lp.m_mu << std::endl;
	(lp.m_lw ? s << "Long way selected" << std::endl : s << "Short way selected" << std::endl);
	s << "Non-dimensional 2D definitions (r1=1,mu=1):" << std::endl;
	s << "Semiperimeter: " << lp.m_s <<std::endl;
	s << "Chord: " << lp.m_c << std::endl;
	s << "Time of Flight: " << lp.m_tof / sqrt(pow(norm(lp.m_r1),3) / lp.m_mu) << std::endl;
	s << std::endl;
	s << "Multiple revolutions active: " << (lp.m_multi_revs? "True":"False") << std::endl;
	if (lp.m_multi_revs) {
		s << "Maximum number of revolutions: " << lp.m_Nmax << std::endl;
	}
	s << "Solutions: " << std::endl;
	s << "0 revs, Iters: " << lp.m_iters[0] << ", a: " << lp.m_a[0] << ", p: " << lp.m_p[0] <<std::endl;
	s <<"\tv1= " << lp.m_v1[0] << " v2= " << lp.m_v2[0] << std::endl;
	for (int i=0; i<lp.m_Nmax;++i)
	{
		s << i+1 << " revs,  left. Iters: " << lp.m_iters[1+2*i] << ", a: " << lp.m_a[1+2*i] << ", p: " << lp.m_p[1+2*i] <<std::endl;
		s << "\tv1= " << lp.m_v1[1+2*i] << " v2= " << lp.m_v2[1+2*i] << std::endl;
		s << i+1 << " revs, right. Iters: " << lp.m_iters[2+2*i] << ", a: " << lp.m_a[2+2*i] << ", p: " << lp.m_p[2+2*i] <<std::endl;
		s << "\tv1= " << lp.m_v1[2+2*i] << " v2= " << lp.m_v2[2+2*i] << std::endl;

	}
	return s;
}

} //namespaces

