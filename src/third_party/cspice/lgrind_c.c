/*

-Procedure lgrind_c (Lagrange polynomial interpolation with derivative)

-Abstract

   Evaluate a Lagrange interpolating polynomial for a specified
   set of coordinate pairs, at a specified abscissa value.
   Return the value of both polynomial and derivative.

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
   #undef    lgrind_c

   void lgrind_c ( SpiceInt            n,
                   ConstSpiceDouble  * xvals,
                   ConstSpiceDouble  * yvals,
                   SpiceDouble       * work,
                   SpiceDouble         x,
                   SpiceDouble       * p,
                   SpiceDouble       * dp )

/*

-Brief_I/O

   Variable  I/O  Description
   --------  ---  --------------------------------------------------
   n          I   Number of points defining the polynomial.
   xvals      I   Abscissa values.
   yvals      I   Ordinate values.
   work      I-O  Work space array.
   x          I   Point at which to interpolate the polynomial.
   p          O   Polynomial value at x.
   dp         O   Polynomial derivative at x.

-Detailed_Input

   n              is the number of points defining the polynomial.
                  The arrays xvals and yvals contain n elements.

   xvals,
   yvals          are arrays of abscissa and ordinate values that
                  together define N ordered pairs. The set of points

                     ( xvals[i], yvals[i] )

                  define the Lagrange polynomial used for
                  interpolation. The elements of xvals must be
                  distinct and in increasing order.

   work           is an n x 2 work space array, where n is the same
                  dimension as that of xvals and yvals. It is used
                  by this routine as a scratch area to hold
                  intermediate results.

   x              is the abscissa value at which the interpolating
                  polynomial is to be evaluated.

-Detailed_Output

   p              is the value at x of the unique polynomial of
                  degree n-1 that fits the points in the plane
                  defined by xvals and yvals.

   dp             is the derivative at x of the interpolating
                  polynomial described above.

-Parameters

   None.

-Exceptions

   1)  The error SPICE(DIVIDEBYZERO) signals from a routine
       in the call tree if any two elements of the array
       xvals are equal.

   2)  The error SPICE(INVALIDSIZE) signals from a routine
       in the call tree if n is less than 1.

   3)  This routine does not attempt to ward off or diagnose
       arithmetic overflows.

-Files

   None.

-Particulars

   Given a set of n distinct abscissa values and corresponding
   ordinate values, there is a unique polynomial of degree n-1, often
   called the `Lagrange polynomial', that fits the graph defined by
   these values. The Lagrange polynomial can be used to interpolate
   the value of a function at a specified point, given a discrete
   set of values of the function.

   Users of this routine must choose the number of points to use
   in their interpolation method. The authors of Reference [1] have
   this to say on the topic:

      Unless there is solid evidence that the interpolating function
      is close in form to the true function f, it is a good idea to
      be cautious about high-order interpolation. We
      enthusiastically endorse interpolations with 3 or 4 points, we
      are perhaps tolerant of 5 or 6; but we rarely go higher than
      that unless there is quite rigorous monitoring of estimated
      errors.

   The same authors offer this warning on the use of the
   interpolating function for extrapolation:

      ...the dangers of extrapolation cannot be overemphasized:
      An interpolating function, which is perforce an extrapolating
      function, will typically go berserk when the argument x is
      outside the range of tabulated values by more than the typical
      spacing of tabulated points.

-Examples

   Example:

      Fit a cubic polynomial through the points

           ( -1, -2 )
           (  0, -7 )
           (  1, -8 )
           (  3, 26 )

      and evaluate this polynomial at x = 2.

      #include <stdio.h>
      #include "SpiceUsr.h"

      int main()
         {

         /.
         Local variables.
         ./
         SpiceDouble      p;
         SpiceDouble      dp;
         SpiceDouble      xvals [] = { -1., 0., 1., 3. };
         SpiceDouble      yvals [] = { -2., -7., -8., 26. };
         SpiceDouble      work  [4*2];
         SpiceInt         n = 4;

         lgrind_c ( n, xvals, yvals, work, 2., &p, &dp );

         printf( "p, dp = %lf %lf\n", p, dp);

         return(0);
        }

      The returned value of P should be 1., since the
      unique cubic polynomial that fits these points is

                     3       2
         f(x)   =   x   +  2x  - 4x  - 7

      The returned value of DP should be 16., since the
      derivative of f(x) is

          '          2
         f (x)  =  3x   +  4x  - 4

      We also could have invoked lgrind_c with the reference

         lgrind_c ( n, xvals, yvals, yvals, 2., &p, &dp );

      if we wished to; in this case yvals would have been
      modified on output.

-Restrictions

   None.

-Literature_References

   [1]  "Numerical Recipes---The Art of Scientific Computing" by
         William H. Press, Brian P. Flannery, Saul A. Teukolsky,
         William T. Vetterling (see sections 3.0 and 3.1).

-Author_and_Institution

    N.J. Bachman    (JPL)
    E.D. Wright     (JPL)

-Version

   -CSPICE Version 1.0.0, 24-AUG-2015 (EDW)

-Index_Entries

   interpolate function using Lagrange polynomial
   Lagrange interpolation

-&
*/

{ /* Begin lgrind_c */

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

   chkin_c ( "lgrind_c" );

   /*
   The f2c'd routine does the work.
   */
   lgrind_( (integer    *) &n,
            (doublereal *) xvals,
            (doublereal *) yvals,
            (doublereal *) work,
            (doublereal *) &x,
            (doublereal *) p,
            (doublereal *) dp);

   chkout_c ( "lgrind_c" );

} /* End lgrind_c */
