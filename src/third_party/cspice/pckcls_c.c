/*

-Procedure pckcls_c ( PCK, close file )

-Abstract
 
   Close an open PCK file. 
 
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
 
   PCK 
 
-Keywords
 
   PCK 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZst.h"

   void pckcls_c ( SpiceInt handle ) 

/*

-Brief_I/O
 
   VARIABLE  I/O  DESCRIPTION 
   --------  ---  -------------------------------------------------- 
   handle     I   Handle of the PCK file to be closed. 
 
-Detailed_Input
 
   handle   The handle of the PCK file that is to be closed. 
 
-Detailed_Output
 
   None. 
 
-Parameters
 
   None. 
 
-Exceptions
 
   1) If there are no segments in the file the error 
      SPICE(NOSEGMENTSFOUND) will be signaled. 
 
-Files
 
   See argument `handle'. 
 
-Particulars
 
   None. 
 
-Examples
 
   Suppose that you want to create a new PCK file called "new.PCK"
   that contains a single type 2 PCK segment and has room for at 
   least 5000 comment characters. The following code fragment should 
   take care of this for you, assuming that all of the variables 
   passed to the PCK type 2 segment writer have appropriate values. 
 
      name   = "new.pck";
      ifname = "Test PCK file";
 
      pckopn_c ( name, ifname, 5000, &handle );

      pckw02_c ( handle, body,   frame_c, first,   last, 
                 segid,  intlen, n,       polydg,  cdata, 
                 btime                                   );

      pckcls_c ( handle ); 

 
-Restrictions
 
    None. 
 
-Literature_References
 
    None. 
 
-Author_and_Institution
 
    N.J. Bachman      (JPL)
    K.R. Gehringer    (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 16-DEC-2016 (NJB) (KRG)

-Index_Entries
 
   close a pck file 
 
-&
*/

{ /* Begin pckcls_c */

 
   /*
   Participate in error tracing.
   */
   chkin_c ( "pckcls_c" );


   pckcls_ ( (integer *) &handle );


   chkout_c ( "pckcls_c" );

} /* End pckcls_c */
