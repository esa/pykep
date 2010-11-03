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

#ifndef REGULA_FALSI_H
#define REGULA_FALSI_H
#include <math.h>

namespace kep_toolbox
{
	/// Regula-Falsi method
	/**
	 * Simple implementation of the Regula-Falsi method for zero finding.
	 *
	 * \param[in] a Starting point 1
	 * \param[in] a Starting point 2
	 * \param[in] F equation to be solved in the form F = 0. Needs to be callable as F(x)
	 *
	 * @author Dario Izzo (dario.izzo@esa.int)
	 */
	template <class my_float, class my_function>
			int regula_falsi(my_float &a, my_float &b, my_function F,
				   int max_loop, const double& accuracy)

{
	int n=0;
	double c=0;
	double Fa=(F)(a);
	double Fb=(F)(b);
	double Fc=Fa;
	while(n<max_loop){
		if(fabs(Fc)<accuracy) break; //Equation solved within accuracy
		c=(a*Fb-b*Fa)/(Fb-Fa);
		Fc = (F)(c);
		n++;
		a = b;
		Fa=Fb;
		b = c;
		Fb=Fc;
	}
	a=c;
	b=c;
	return n;
}
}

/* Main program to test the method
double func_1(double x) // root is 1.85792
{return (cosh(x) + cos(x) - 3.0); }

int main()
{
	double a = 1; double b=3; int iter;
	iter = ::kep_toolbox::regula_falsi(a,b,func_1,100,1e-4);
	std::cout << a << " " << b << " " << iter << std::endl;
	return 0;
}
*/

#endif // REGULA_FALSI_H
