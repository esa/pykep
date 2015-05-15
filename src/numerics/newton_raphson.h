/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
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

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H
#include <math.h>

namespace kep_toolbox
{
	/// Newton-Raphson method
	/**
	 * Standard implementation of the Newton-Raphson method to solve a non-linear equation
	 *
	 * \param[in] x Starting point
	 * \param[in] F equation to be solved in the form F = 0. Needs to be callable as F(x)
	 * \param[in] dF derivative of F with respect to x.  Needs to be callable as dF(x)
	 *
	 * @author Dario Izzo (dario.izzo@esa.int)
	 */
	template <class my_float, class my_functionA, class my_functionB>
			int newton_raphson(my_float &x, my_functionA F, my_functionB dF,
				   int max_loop, const double& accuracy)
	{
		my_float term;
		//double start = x; //DS: unused
		//main iteration
		do
		{
			term = (F)(x) / (dF)(x);
			x = x - term;
		}
		// check if term is within required accuracy or loop limit is exceeded
		while ((fabs(term / std::max(std::fabs(x), 1.)) > accuracy) && (--max_loop));
		return max_loop;
	}
}
/*
//----------------------------------------------------------------------------//
// test functions

double func_1(double x)                       // root is 1.85792
  {  return (cosh(x) + cos(x) - 3.0); }       // f(x) = cosh(x) + cos(x) - 3 = 0

double fdiv_1(double x)
  {  return (sinh(x) - sin(x));  }            // f'(x) = sinh(x) - sin(x)

double func_2(double x)                       // root is 5.0
  {  return (x*x - 25.0);  }                  // f(x) = x * x - 25 = 0

double fdiv_2(double x)
  {  return (2.0 * x);  }                     // f'(x) = 2x

complex<double> func_3(complex<double> x)     // roots 5 + or - 3
  { return  x*x - 10.0*x + 34.0;   }          // f(x) x^2 - 10x + 34

complex<double> fdiv_3(complex<double> x)
  { return  2.0*x -10.0;   }                  // f'(x) 2x - 10


struct f3 {
  complex<double> operator()(complex<double> x) const
  {
    return  x*x - 10.0*x + 34.0;
  }
};

struct fd3 {
  complex<double> operator()(complex<double> x) const
  {
    return  2.0*x -10.0;
  }
};


double func_4(double x)                             // three real roots 4, -3, 1
  {  return 2*x*x*x - 4*x*x - 22*x + 24 ;  }    // f(x) = 2x^3 - 4x^2 - 22x + 24

double fdiv_4(double x)
  {  return 6*x*x - 8*x - 22;  }                       // f'(x) = 6x^2 - 8x - 22

double func_5(double E, double M, double e)
    {        return (E - e * sin(E) - M);        }

double fdiv_5(double E, double e)
    {        return (1 - e * cos(E));        }

//----------------------------------------------------------------------------//
// Main program to test above function
int main()
{
    cout << "\nFind root of f(x) = cosh(x) + cos(x) - 3 = 0";
    double x = 1.0;                                   // initial 'guess' at root
    if (newton(x, func_1, fdiv_1, 100, 1.0e-8))
	  cout << "\n   root x = " << x << ", test of f(x) = " << func_1(x);
    else  cout << "\n   failed to find root ";

    cout << "\n\nFind root of f(x) = x * x - 25 = 0";
    x = 1.0;                                          // initial 'guess' at root
    if ( newton(x, func_2, fdiv_2, 100, 1.0e-8))
	  cout << "\n   root x = " << x << ", test of f(x) = " << func_2(x);
    else  cout << "\n   failed to find root ";

    cout << "\n\nFind root of f(x) = x^2 - 10x + 34  = 0";
    complex<double> xc = complex<double>(1.0, 1.0);   // initial 'guess' at root
    if ( newton(xc, func_3, fdiv_3, 100, 1.0e-8))
	  cout << "\n   root x = " << xc << ", test of f(x) = " << func_3(xc);
    else  cout << "\n   failed to find root ";

    cout << "\n\nFind root of f(x) = x^2 - 10x + 34  = 0";
    xc = complex<double>(1.0, -1.0);                  // initial 'guess' at root
    if ( newton(xc, func_3, fdiv_3, 100, 1.0e-8))
	  cout << "\n   root x = " << xc << ", test of f(x) = " << func_3(xc);
    else  cout << "\n   failed to find root ";

    cout << "\n\nFind root of f(x) = 2x^3 - 4x^2 - 22x + 24 = 0";
    x = 5.0;                                          // initial 'guess' at root
    if ( newton(x, func_4, fdiv_4, 100, 1.0e-8) )
	  cout << "\n   root x = " << x << ", test of f(x) = " << func_4(x);
    else  cout << "\n   failed to find root ";

    cout << "\n\nFind root of f(x) = E - e sin(E) - M";
    x = 0;                                          // initial 'guess' at root
    newton(x, boost::bind(func_5, _1, 0.2, 0.4), boost::bind(fdiv_5, _1, 0.4), 100, 1.0e-8);
    if ( newton(x, boost::bind(func_5, _1, 0.2, 0.4), boost::bind(fdiv_5, _1, 0.4), 100, 1.0e-8) )
    cout << "\n   root x = " << x << ", test of f(x) = " << boost::bind(func_5, _1, 0.2,0.4)(x);
    else  cout << "\n   failed to find root ";
    cout << endl;
    return 0;
}
*/

#endif // NEWTON_RAPHSON_H
