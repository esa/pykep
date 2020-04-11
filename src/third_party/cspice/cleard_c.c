/*

-Procedure cleard_c ( Clear a double precision array )

-Abstract
 
   Fill a double precision array with zeros. 
 
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
 
   ARRAY
   ASSIGNMENT 
 
*/
   #include "SpiceUsr.h"


   void cleard_c ( SpiceInt       ndim,
                   SpiceDouble  * array ) 
/*

-Brief_I/O
 
   VARIABLE  I/O   DESCRIPTION 
   --------  ---  --------------------------------------- 
   ndim       I   The number of elements of `array' which are to be 
                  set to zero. 
   array      O   Double precision array to be filled. 
 
-Detailed_Input
 
   ndim       is the number of elements in `array' which are to be 
              set to zero. 
 
-Detailed_Output
 
   array      is the double precision array which is to be filled 
              with zeros. 
 
-Parameters
 
   None. 
 
-Exceptions
 
   Error free. 
 
   1) If ndim < 1, the array is not modified. 
 
-Files
 
   None. 
 
-Particulars
 
   None. 
 
-Examples
 
   If ndim = 4, then the contents of `array' are: 
 
   array[0] = 0.0 
   array[1] = 0.0 
   array[2] = 0.0 
   array[3] = 0.0 

-Restrictions
 
   None. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman    (JPL)
   W.M. Owen       (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 21-MAR-2016 (NJB) (WMO)

-Index_Entries
 
   clear a d.p. array 
 
-&
*/

{ /* Begin cleard_c */

   /*
   Local variables
   */
   SpiceInt                i;


   for ( i = 0;  i < ndim;  i++ )
   {
      array[i] = 0.0;
   }


} /* End cleard_c */
