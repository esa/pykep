/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
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

#ifndef KEP_TOOLBOX_TAYLOR_H_
#define KEP_TOOLBOX_TAYLOR_H_
typedef double MY_FLOAT;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <keplerian_toolbox/detail/visibility.hpp>

/*
 *  MY_FLOAT is the data type to be used in computing derivatives.
 *  It may be 'float', 'double' or user defined private data types
 *  like 'doubledouble', 'complex' etc.
 */

/* for double or doubledouble, don't need to initialize */
#define InitMyFloat(r)
#define ClearMyFloat(r)

/* assign b to a */
#define AssignMyFloat(a, b)                                                                                            \
    {                                                                                                                  \
        (a) = (b);                                                                                                     \
    }

/* create a MY_FLOAT from a, assign to r. a is an integer or a float */
#define MakeMyFloatA(r, a) (r = (double)(a))

/* create a MY_FLOAT from string, a is an integer or a float, s is its string representation */
#define MakeMyFloatC(r, s, a) (r = (double)(a))

/* create a MY_FLOAT from a, assign to r and return r */
#define MakeMyFloatB(r, a) (r = (double)(a), r)

/* addition r=a+b   */
#define AddMyFloatA(r, a, b) (r = (a) + (b))

/* substraction r=a-b */
#define SubstractMyFloatA(r, a, b) (r = (a) - (b))

/* multiplication r=a*b */
#define MultiplyMyFloatA(r, a, b) (r = (a) * (b))

/* division r=a/b */
#define DivideMyFloatA(r, a, b) (r = (a) / (b))

/* division by an integer r=a/i */
#define DivideMyFloatByInt(r, a, i) (r = (a) / (double)(i))
/* negation r=-a*/
#define NegateMyFloatA(r, a) (r = -(a))

/* square root r=sqrt(a) */
#define sqrtMyFloatA(r, a) (r = sqrt(a))
/* exponentiation r=a^b */
#define ExponentiateMyFloatA(r, a, b) (r = pow((a), (b)))
/* exponentiation r=a^b, b is an integer */
#define ExponentiateMyFloatIA(r, a, b) (r = pow((a), (double)(b)))
/* sin(a)  r=sin(a) */
#define sinMyFloatA(r, a) (r = sin((a)))
/* cos(a)  r=cos(a) */
#define cosMyFloatA(r, a) (r = cos((a)))
/* tan(a)  r=tan(a) */
#define tanMyFloatA(r, a) (r = tan((a)))
/* atan(a) r=atan(a) */
#define atanMyFloatA(r, a) (r = atan((a)))
/* exp(a)  r=exp(a) */
#define expMyFloatA(r, a) (r = exp((a)))
/* log(a)  r=log(a) */
#define logMyFloatA(r, a) (r = log((a)))
/* sinh(a) r=sinh(a) */
#define sinhMyFloatA(r, a) (r = sinh(a))
/* cosh(a) r=cosh(a) */
#define coshMyFloatA(r, a) (r = cosh(a))
/* tanh(a) r=tanh(a) */
#define tanhMyFloatA(r, a) (r = tanh(a))

/* log10(a)  r=log10(a) */
#define log10MyFloatA(r, a) (r = log10((a)))
/* fabs(a) r=fabs(a) */
#define fabsMyFloatA(r, a) (r = fabs(a))

/* convert to int */
#define MyFloatToInt(ir, fa) (ir = (int)(fa))
/* convert to double */
#define MyFloatToDouble(ir, fa) (ir = (double)(fa))

/* boolean operation  */
#define MyFloatA_GE_B(a, b) ((a) >= (b))
#define MyFloatA_GT_B(a, b) ((a) > (b))
#define MyFloatA_LE_B(a, b) ((a) <= (b))
#define MyFloatA_LT_B(a, b) ((a) < (b))
#define MyFloatA_EQ_B(a, b) ((a) == (b))
#define MyFloatA_NEQ_B(a, b) ((a) != (b))

/************************************************************************/

#endif // KEP_TOOLBOX_TAYLOR_H_

MY_FLOAT **taylor_coefficients_fixed_thrust(MY_FLOAT t, MY_FLOAT *x, int order, double mu, double veff, double ux,
                                            double uy, double uz);

MY_FLOAT **taylor_coefficients_fixed_thrustA(MY_FLOAT t, MY_FLOAT *x, int order, int reuse_last_computation, double mu,
                                             double veff, double ux, double uy, double uz);

KEP_TOOLBOX_DLL_PUBLIC int taylor_step_fixed_thrust(MY_FLOAT *ti, MY_FLOAT *x, int dir, int step_ctl, double log10abserr,
                                                double log10relerr, MY_FLOAT *endtime, MY_FLOAT *ht, int *order,
                                                double mu, double veff, double ux, double uy, double uz);
