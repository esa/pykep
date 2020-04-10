/*

-Procedure oscltx_c ( Extended osculating elements from state )

-Abstract
 
   Determine the set of osculating conic orbital elements that 
   corresponds to the state (position, velocity) of a body at some 
   epoch. In additional to the classical elements, return the true 
   anomaly, semi-major axis, and period, if applicable. 
 
-Disclaimer
 
   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE 
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. 
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE 
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" 
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY 
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A 
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC 
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
 
   CONIC 
   ELEMENTS 
   EPHEMERIS 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #undef   oscltx_c

   void oscltx_c ( ConstSpiceDouble state [6],
                   SpiceDouble      et,
                   SpiceDouble      mu,
                   SpiceDouble      elts  [SPICE_OSCLTX_NELTS] ) 

/*

-Brief_I/O
 
   VARIABLE  I/O  DESCRIPTION 
   --------  ---  -------------------------------------------------- 
   state      I   State of body at epoch of elements. 
   et         I   Epoch of elements. 
   mu         I   Gravitational parameter (GM) of primary body. 
   elts       O   Extended set of classical conic elements. 
 
-Detailed_Input
 
   state      is the state (position and velocity) of the body 
              at some epoch. Components are x, y, z, dx/dt, dy/dt, 
              dz/dt. `state' must be expressed relative to an  
              inertial reference frame. Units are km and km/sec. 
 
 
   et         is the epoch of the input state, in ephemeris seconds 
              past J2000. 
 
                                                     3    2 
   mu         is the gravitational parameter (GM, km /sec ) of 
              the primary body. 
 
-Detailed_Output
 
   elts        are equivalent conic elements describing the orbit 
               of the body around its primary. The elements are, 
               in order: 
 
                  RP      Perifocal distance. 
                  ECC     Eccentricity. 
                  INC     Inclination. 
                  LNODE   Longitude of the ascending node. 
                  ARGP    Argument of periapsis. 
                  M0      Mean anomaly at epoch. 
                  T0      Epoch. 
                  MU      Gravitational parameter. 
                  NU      True anomaly at epoch. 
                  A       Semi-major axis. A is set to zero if 
                          it is not computable. 
                  TAU     Orbital period. Applicable only for 
                          elliptical orbits. Set to zero otherwise. 
 
               The epoch of the elements is the epoch of the input 
               state. Units are km, rad, rad/sec. The same elements 
               are used to describe all three types (elliptic, 
               hyperbolic, and parabolic) of conic orbits. 

               See the Parameters section for information on the
               declaration of `elts'.

 
-Parameters
 
   SPICE_OSCLTX_NELTS

               is the length of the output array `elts'.

               `elts' is intended to contain unused space to hold
               additional elements that may be added in a later version
               of this routine. In order to maintain forward
               compatibility, user applications should declare `elts'
               as follows:

                  SpiceDouble    elts[SPICE_OSCLTX_NELTS];

  
-Exceptions
 
   1) If MU is not positive, the error SPICE(NONPOSITIVEMASS) 
      is signaled. 
 
   2) If the specific angular momentum vector derived from `state' 
      is the zero vector, the error SPICE(DEGENERATECASE) 
      is signaled. 
 
   3) If the position or velocity vectors derived from `state' 
      is the zero vector, the error SPICE(DEGENERATECASE) 
      is signaled. 
 
   4) If the inclination is determined to be zero or 180 degrees, 
      the longitude of the ascending node is set to zero.   
 
   5) If the eccentricity is determined to be zero, the argument of 
      periapse is set to zero.      
    
   6) If the eccentricity of the orbit is very close to but not 
      equal to zero, the argument of periapse may not be accurately 
      determined. 
 
   7) For inclinations near but not equal to 0 or 180 degrees, 
      the longitude of the ascending node may not be determined 
      accurately.  The argument of periapse and mean anomaly may 
      also be inaccurate. 
 
   8) For eccentricities very close to but not equal to 1, the 
      results of this routine are unreliable.  
 
   9) If the specific angular momentum vector is non-zero but 
      "close" to zero, the results of this routine are unreliable. 
 
  10) If `state' is expressed relative to a non-inertial reference 
      frame, the resulting elements are invalid.  No error checking 
      is done to detect this problem. 
 
  11) The semi-major axis and period may not be computable for 
      orbits having eccentricity too close to 1. If the semi-major 
      axis is not computable, both it and the period are set to zero. 
      If the period is not computable, it is set to zero. 
 
-Files
 
   None. 
 
-Particulars
 
   This routine returns in the first 8 elements of the array `elts' 
   the outputs computed by oscelt_c, and in addition returns in 
   elements 9-11 the quantities: 
 
      elts[8]    true anomaly at `et', in radians. 
 
      elts[9]   orbital semi-major axis at `et', in km. Valid 
                if and only if this value is non-zero.
 
                The semi-major axis won't be computable if the 
                eccentricity of the orbit is too close to 1. 
                In this case A is set to zero. 
  
      elts[10]  orbital period. If the period is not computable, 
                TAU is set to zero. 
 
   The CSPICE routine conics_c is an approximate inverse of this 
   routine: conics_c maps a set of osculating elements and a time to a 
   state vector. 
 
-Examples
 
   Let vinit contain the initial state of a spacecraft relative to
   the center of a planet at epoch `et', and let GM be the gravitation
   parameter of the planet. The call

      oscltx_c ( vinit, et, gm, elts );

   produces a set of osculating elements describing the nominal
   orbit that the spacecraft would follow in the absence of all
   other bodies in the solar system.

   Now let state contain the state of the same spacecraft at some
   other epoch, later. The difference between this state and the
   state predicted by the nominal orbit at the same epoch can be
   computed as follows.

      conics_c ( elts,    later, nominal );
      vsubg_c  ( nominal, state, 6, diff );

      printf( "Perturbation in x, dx/dt = %e %e\n", diff[0], diff[3] );
      printf( "                y, dy/dt = %e %e\n", diff[1], diff[4] );
      printf( "                z, dz/dt = %e %e\n", diff[2], diff[5] );
 
-Restrictions
 
   1) The input state vector must be expressed relative to an 
      inertial reference frame. 
 
   2) Osculating elements are generally not useful for 
      high-accuracy work. 
 
   3) Accurate osculating elements may be difficult to derive for 
      near-circular or near-equatorial orbits. Osculating elements 
      for such orbits should be used with caution. 
 
   4) Extracting osculating elements from a state vector is a  
      mathematically simple but numerically challenging task.  The 
      mapping from a state vector to equivalent elements is 
      undefined for certain state vectors, and the mapping is 
      difficult to implement with finite precision arithmetic for 
      states near the subsets of R6 where singularities occur. 
 
      In general, the elements found by this routine can have 
      two kinds of problems: 
 
         - The elements are not accurate but still represent 
           the input state accurately.  The can happen in 
           cases where the inclination is near zero or 180 
           degrees, or for near-circular orbits. 
 
         - The elements are garbage.  This can occur when 
           the eccentricity of the orbit is close to but 
           not equal to 1. In general, any inputs that cause 
           great loss of precision in the computation of the 
           specific angular momentum vector or the eccentricity 
           vector will result in invalid outputs. 
 
      For further details, see the Exceptions section. 
 
      Users of this routine should carefully consider whether 
      it is suitable for their applications.  One recommended 
      "sanity check" on the outputs is to supply them to the 
      CSPICE routine conics_c and compare the resulting state  
      vector with the one supplied to this routine. 
 
-Literature_References
 
   [1] Roger Bate, Fundamentals of Astrodynamics, Dover, 1971. 
 
-Author_and_Institution
 
   N.J. Bachman    (JPL) 
   K.R. Gehringer  (JPL) 
   I.M. Underwood  (JPL) 
   E.D. Wright     (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 25-JAN-2017 (NJB) (KRG) (IMU) (EDW)
   
      Original version 11-NOV-2014 (NJB) (KRG) (IMU) (EDW)

-Index_Entries
 
   extended conic elements from state 
   extended osculating elements from state 
   convert state to extended osculating elements 
 
-&
*/

{ /* Begin oscltx_c */


   /*
   Participate in error tracing.
   */
   chkin_c ( "oscltx_c" );


   /*
   Let the f2c'd routine do the work. 
   */
   oscltx_ (  (doublereal *) state,
              (doublereal *) &et,
              (doublereal *) &mu,
              (doublereal *) elts  );


   chkout_c ( "oscltx_c" );

} /* End oscltx_c */
