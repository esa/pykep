/*

-Procedure pltnrm_c ( DSK, compute outward normal of plate )

-Abstract
 
   Compute an outward normal vector of a triangular plate. 
   The vector does not necessarily have unit length. 
 
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
 
   DSK 
 
-Keywords
 
   DSK 
   FILES 
   TOPOGRAPHY 
  
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZmc.h"
   #undef   pltnrm_c


   void pltnrm_c ( ConstSpiceDouble    v1[3],
                   ConstSpiceDouble    v2[3],
                   ConstSpiceDouble    v3[3],
                   SpiceDouble         normal[3] ) 
/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   v1, 
   v2, 
   v3         I   Vertices of a plate. 
   normal     O   Plate's outward normal vector. 
 
-Detailed_Input
 
   v1, 
   v2, 
   v3             are vertices of a triangular plate. 
                   
-Detailed_Output
 
   normal         is an outward normal vector of the plate defined by 
                  the input vertices. The order of the vertices is 
                  used to determine the choice of normal direction: 
                  the normal vector is 
 
                     ( V2 - V1 ) x ( V3 - V2 ) 
                    
-Parameters
 
   None. 
 
-Exceptions
 
   1) The input plate may be degenerate: it may be a line segment 
      or a point. These are not considered to be erroneous inputs. 
 
-Files
 
   None. 
 
-Particulars
 
   This routine saves computation time by not scaling the output 
   vector to unit length. The caller can scale the vector using 
   the routine vhat_c. 
 
-Examples
 
 
   1) Compute an upward normal of an equilateral triangle lying 
      in the X-Y plane and centered at the origin. 
 
      Example code begins here. 
 
         #include <stdio.h>
         #include <math.h>
         #include "SpiceUsr.h"

         int main()
         {
            ./
            Local variables 
            ./ 
            SpiceDouble             normal[3];
            SpiceDouble             s;
            SpiceDouble             v1[3];
            SpiceDouble             v2[3];
            SpiceDouble             v3[3];

            s = sqrt(3.0)/2;

            vpack_c (   s,  -0.5,  0.0, v1 );
            vpack_c ( 0.0,   1.0,  0.0, v2 );
            vpack_c (  -s,  -0.5,  0.0, v3 );

            pltnrm_c ( v1, v2, v3, normal );

            printf ( "NORMAL = %18.11e %18.11e %18.11e\n", 
                     normal[0], normal[1], normal[2]      );

            return ( 0 );
         }

  
      When run on a PC/Linux/gcc/64-bit platform, the output 
      from this program was: 
 
 
   NORMAL =  0.000000000000E+00  0.000000000000E+00  0.259807621135E+01 
 
 
-Restrictions
 
   None. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman    (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 26-JAN-2016 (NJB)

-Index_Entries

   compute normal vector of triangular plate from vertices
 
-&
*/

{ /* Begin pltnrm_c */

   
   /*
   Participate in error tracing.
   */
   chkin_c ( "pltnrm_c" );

   pltnrm_ ( (SpiceDouble *) v1,
             (SpiceDouble *) v2,
             (SpiceDouble *) v3,
             (SpiceDouble *) normal );

   chkout_c ( "pltnrm_c" );

} /* End pltnrm_c */
