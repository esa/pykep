/*

-Procedure chbder_c ( Derivatives of a Chebyshev expansion )

-Abstract

   Given the coefficients for the Chebyshev expansion of a
   polynomial, this returns the value of the polynomial and its
   first nderiv derivatives evaluated at the input X.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS set_c FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   None.

-Keywords

    INTERPOLATION
    MATH
    POLYNOMIAL

*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #undef    chbder_c

   void chbder_c ( ConstSpiceDouble * cp,
                   SpiceInt           degp,
                   SpiceDouble        x2s[2],
                   SpiceDouble        x,
                   SpiceInt           nderiv,
                   SpiceDouble      * partdp,
                   SpiceDouble      * dpdxs )

/*

-Brief_I/O

    VARIABLE  I/O  DESCRIPTION
    --------  ---  --------------------------------------------------
    cp         I   degp+1 Chebyshev polynomial coefficients.
    degp       I   Degree of polynomial.
    x2s        I   Transformation parameters of polynomial.
    x          I   Value for which the polynomial is to be evaluated
    nderiv     I   The number of derivatives to compute
    partdp    I-O  Workspace provided for computing derivatives
    dpdxs      O   Array of the derivatives of the polynomial

-Detailed_Input

    cp         is an array of coefficients a polynomial with respect
               to the Chebyshev basis. The polynomial to be
               evaluated is assumed to be of the form:

                 cp(degp+1)*T(degp,S) + cp(degp)*T(degp-1,S) + ...

                                  ... + cp(2)*T(1,S) + cp(1)*T(0,S)

               where T(I,S) is the I'th Chebyshev polynomial
               evaluated  at a number S whose double precision
               value lies between -1 and 1. The value of S is
               computed from the input variables x2s[0],
               x2s[1] and X.

    degp       is the degree of the Chebyshev polynomial to be
               evaluated.

    x2s        is an array of two parameters. These parameters are
               used to transform the domain of the input variable X
               into the standard domain of the Chebyshev polynomial.
               x2s[0] should be a reference point in the domain of
               x; x2s[1] should be the radius by which points are
               allowed to deviate from the reference point and while
               remaining within the domain of x. The value of
               x is transformed into the value S given by

                         S = ( x - x2s[0] ) / x2s[1]

               Typically x2s[0] is the midpoint of the interval over
               which x is allowed to vary and x2s[1] is the radius
               of the interval.

               The main reason for doing this is that a Chebyshev
               expansion is usually fit to data over a span
               from A to B where A and B are not -1 and 1
               respectively. Thus to get the "best fit" the
               data was transformed to the interval [-1,1] and
               coefficients generated. These coefficients are
               not re-scaled to the interval of the data so that
               the numerical "robustness" of the Chebyshev fit will
               not be lost. Consequently, when the "best fitting"
               polynomial needs to be evaluated at an intermediate
               point, the point of evaluation must be transformed
               in the same way that the generating points were
               transformed.

    x          Value for which the polynomial is to be evaluated.

    nderiv     is the number of derivatives to be computed by the
               routine. nderiv should be non-negative.

    partdp     Is a work space used by the program to compute
               all of the desired derivatives. It should be declared
               in the calling program as

                      SpiceDouble partdp[3 * (nderiv+1)]

-Detailed_Output

    dpdxs(0)   The value of the polynomial to be evaluated. It
               is given by

                 cp(degp+1)*T(degp,S) + cp(degp)*T(degp-1,S) + ...

                                  ... + cp(2)*T(1,S) + cp(1)*T(0,S)

               where T(I,S) is the I'th Chebyshev polynomial
               evaluated  at a number S = ( x -  x2s[0] )/ x2s[1].

    dpdxs(i)   The value of the i'th derivative of the polynomial at
               x. (I ranges from 1 to nderiv) It is given by

                                           [i]
                 (1/x2s[1]**i) ( cp(degp+1)*T   (degp,S)

                                         [i]
                             + cp(degp)*T   (degp-1,S) + ...

                             .
                             .
                             .
                                      [i]
                         ... + cp(2)*T   (1,S)

                                      [i]
                             + cp(1)*T   (0,S) )

                                 [i]
               where T(k,S) and T   (i,S)  are the k'th Chebyshev
               polynomial and its i'th derivative respectively,
               evaluated at the number S = ( x - x2s[0] )/x2s[1].

-Parameters

    None.

-Exceptions

   Error free

   No tests are performed for exceptional values ( nderiv negative,
   degp negative, etc.) This routine is expected to be used at a low
   level in ephemeris evaluations. For that reason it has been
   elected as a routine that will not participate in error handling.

-Files

    None.

-Particulars

    This routine computes the value of a Chebyshev polynomial
    expansion and the derivatives of the expansion with respect to X.
    The polynomial is given by

         cp(degp+1)*T(degp,S) + cp(degp)*T(degp-1,S) + ...

                          ... + cp(2)*T(1,S) + cp(1)*T(0,S)

    where

         S  =  ( x - x2s[0] ) / x2s[1]

    and

         T(i,S) is the i'th Chebyshev polynomial of the first kind
         evaluated at S.

-Examples

   Depending upon the user's needs, there are 3 routines available
   for evaluating Chebyshev polynomials.

      chbval_c for evaluating a Chebyshev polynomial when no
             derivatives are desired.

      chbint_c for evaluating a Chebyshev polynomial and its
             first derivative.

      chbder_c for evaluating a Chebyshev polynomial and a user
             or application dependent number of derivatives.

   Of these 3 the one most commonly employed by NAIF software
   is chbint_c as it is used to interpolate ephemeris state
   vectors which requires the evaluation of a polynomial
   and its derivative. When no derivatives are desired one
   should use chbval_c, or when more than one or an unknown
   number of derivatives are desired one should use chbder_c.

   Example:

      #include <stdio.h>
      #include "SpiceUsr.h"
   
      int main()
         {
   
         /.
         Local variables.
         ./
         SpiceDouble      cp [] = {  1., 3., 0.5, 1., 0.5, -1., 1. };
         SpiceInt         degp    = 6;
         SpiceInt         nderiv = 3;
         SpiceDouble      x2s[]  = { .5, 3.};
         SpiceDouble      x      = 1.;
   
         /. Dimension partdp as 3 * (nderiv + 1) ./
         SpiceDouble      partdp[3 * 4];
   
        /. Dimension dpdxs as nderiv + 1. ./
        SpiceDouble       dpdxs [3+1];
   
         int              i;
   
         chbder_c ( cp, degp, x2s, x, nderiv, partdp, dpdxs );
   
         for ( i=0; i<=nderiv; i++ )
            {
            printf( "dpdxs = %lf\n", dpdxs[i] );
            }
   
         return(0);
         }

   The program outputs:

      dpdxs = -0.340878
      dpdxs = 0.382716
      dpdxs = 4.288066
      dpdxs = -1.514403

-Restrictions

    The user must be sure that the provided workspace is declared
    properly in the calling routine. The proper declaration is:

    SpiceInt         nderiv = the desired number of derivatives;
    SpiceDouble      partdp[3 * (nderiv + 1)];

    If for some reason a parameter is not passed to this routine in
    nderiv, the user should make sure that the value of nderiv is not
    so large that the work space provided is inadequate.

    One needs to be careful that the value (X-x2s[0]) / x2s[1] lies
    between -1 and 1. Otherwise, the routine may fail spectacularly
    (for example with a floating point overflow).

    While this routine will compute derivatives of the input
    polynomial, the user should consider how accurately the
    derivatives of the Chebyshev fit, match the derivatives of the
    function it approximates.

-Literature_References

    "Numerical Recipes -- The Art of Scientific Computing" by
     William H. Press, Brian P. Flannery, Saul A. Teukolsky,
     Willam T. Vetterling. (See Clenshaw's Recurrence Formula)

    "The Chebyshev Polynomials" by Theodore J. Rivlin.

    "CRC Handbook of Tables for Mathematics"

-Author_and_Institution

    N.J. Bachman    (JPL)
    W.L. Taber      (JPL)
    E.D. Wright     (JPL)

-Version

   -CSPICE Version 1.0.0, 24-AUG-2015 (EDW)

-Index_Entries

   derivatives of a chebyshev expansion

-&
*/

{ /* Begin chbder_c */

   /*
   Local constants
   */


   /*
   Local macros
   */


   /*
   Local variables
   */


   /*
   Static variables
   */


   /*
   Participate in error tracing.
   */

   chkin_c ( "chbder_c" );

   chbder_( (doublereal *) cp,
            (integer    *) &degp,
            (doublereal *) x2s,
            (doublereal *) &x,
            (integer    *) &nderiv,
            (doublereal *) partdp,
            (doublereal *) dpdxs);

   chkout_c ( "chbder_c" );

} /* End chbder_c */
