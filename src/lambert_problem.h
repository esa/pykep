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

#ifndef KEP_TOOLBOX_LAMBERT_PROBLEM_H
#define KEP_TOOLBOX_LAMBERT_PROBLEM_H

#include <cmath>
#include <vector>

#include "astro_constants.h"
#include "serialization.h"
#include "config.h"




namespace kep_toolbox {

/// Lambert Problem
/**
 * This class represent a Lambert's problem. When instantiated it assumes a prograde orbit (unless otherwise stated)
 * and evaluates all the solutions up to a maximum number of multiple revolutions.
 * After the object is instantiated the solutions can be retreived using the appropriate getters. Note that the
 * number of solutions will be N_max*2 + 1, where N_max is the maximum number of revolutions.
 *
 * NOTE: The class has been tested extensively via monte carlo runs checked with numerical propagation. Compared
 * to the previous Lambert Solver in the keplerian_toolbox it is 1.7 times faster (on average as defined
 * by lambert_test.cpp). With respect to Gooding algorithm it is 1.3 - 1.5 times faster (zero revs - multi revs).
 * The algorithm is described in detail in the publication below and its original with the author.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE lambert_problem
{
	static const array3D default_r1;
	static const array3D default_r2;
public:
	friend std::ostream &operator<<(std::ostream &, const lambert_problem &);
	lambert_problem(const array3D &r1 = default_r1, const array3D &r2 = default_r2, const double &tof = boost::math::constants::pi<double>()/2, const double& mu = 1., const int &cw = 0, const int &multi_revs=5);
	const std::vector<array3D>& get_v1() const;
	const std::vector<array3D>& get_v2() const;
	const array3D& get_r1() const;
	const array3D& get_r2() const;
	const double& get_tof() const;
	const double& get_mu() const;
	const std::vector<double>& get_x() const;
	const std::vector<int>& get_iters() const;
	int get_Nmax() const;
private:
	int householder(const double T, double& x0, const int N, const double eps, const int itermax);
	void dTdx(double &DT,double &DDT,double &DDDT,const double x0, const double tof);
	void x2tof(double &tof,const double x0, const int N);
	void x2tof2(double &tof,const double x0, const int N);
	double hypergeometricF(double z, double tol);
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & const_cast<array3D&> (m_r1);
		ar & const_cast<array3D&> (m_r2);
		ar & const_cast<double&> (m_tof);
		ar & const_cast<double&> (m_mu);
		ar & m_v1;
		ar & m_v2;
		ar & m_iters;
		ar & m_x;
		ar & m_s;
		ar & m_c;
		ar & m_lambda;
		ar & m_iters;
		ar & m_Nmax;
		ar & m_has_converged;
		ar & m_multi_revs;
	}


	const array3D m_r1, m_r2;
	const double m_tof;
	const double m_mu;
	std::vector<array3D> m_v1;
	std::vector<array3D> m_v2;
	std::vector<int> m_iters;
	std::vector<double> m_x;
	double m_s,m_c,m_lambda;
	int m_Nmax;
	bool m_has_converged;
	int m_multi_revs;

};
__KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &, const lambert_problem &);
} //namespaces

#endif // KEP_TOOLBOX_LAMBERT_PROBLEM_H
