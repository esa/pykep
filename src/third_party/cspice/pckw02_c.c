/*

-Procedure pckw02_c ( PCK, write type 2 segment )

-Abstract
 
  Write a type 2 segment to a PCK binary file given the file handle,
  frame class ID, base frame, time range covered by the segment, and
  the Chebyshev polynomial coefficients.
 
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

   DAF
   FRAMES 
   NAIF_IDS 
   PCK 
   SPC 
 
-Keywords
 
   PCK 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZmc.h"
   #include "SpiceZst.h"
   #undef   pckw02_c


   void pckw02_c ( SpiceInt           handle,
                   SpiceInt           clssid,
                   ConstSpiceChar   * frame,
                   SpiceDouble        first,
                   SpiceDouble        last,
                   ConstSpiceChar   * segid,
                   SpiceDouble        intlen,
                   SpiceInt           n,
                   SpiceInt           polydg,
                   SpiceDouble        cdata  [],
                   SpiceDouble        btime      )
/*

-Brief_I/O
 
 Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   handle     I   Handle of binary PCK file open for writing. 
   clssid     I   Frame class ID of body-fixed frame.
   frame      I   Name of base reference frame.
   first      I   Start time of interval covered by segment. 
   last       I   End time of interval covered by segment. 
   segid      I   Segment identifier. 
   intlen     I   Length of time covered by logical record. 
   n          I   Number of logical records in segment. 
   polydg     I   Chebyshev polynomial degree. 
   cdata      I   Array of Chebyshev coefficients. 
   btime      I   Begin time of first logical record. 
 
-Detailed_Input
 
   handle         is the DAF handle of an PCK file to which a type 2 
                  segment is to be added. The PCK file must be open 
                  for writing. 
 
   clssid         is the frame class ID of a body-fixed reference frame
                  whose orientation is described by the segment to be
                  created.
 
   frame          is the NAIF name for a reference frame relative to 
                  which the orientation information for the frame
                  designated by `clssid' is specified. The frame
                  designated by `frame' is called the "base frame."
 
   first, 
   last           are, respectively, the start and stop times of 
                  the time interval over which the segment defines 
                  the orientation of clssid. 
 
   segid          is the segment identifier. A PCK segment 
                  identifier may contain up to 40 characters. 
 
   intlen         Length of time, in seconds, covered by each set of 
                  Chebyshev polynomial coefficients (each logical 
                  record). Each set of Chebyshev coefficients must 
                  cover this fixed time interval, `intlen'. 
 
   n              is the number of sets of Chebyshev polynomial 
                  coefficients (number of logical records) 
                  to be stored in the segment. There is one set 
                  of Chebyshev coefficients for each time interval. 
 
   polydg         Degree of each set of Chebyshev polynomials. 
 
   cdata          Array containing all the sets of Chebyshev 
                  polynomial coefficients to be contained in the 
                  segment of the PCK file.  The coefficients are 
                  stored in cdata in order as follows: 
 
                     the (degree + 1) coefficients for the first 
                     Euler angle of the first logical record 
 
                     the coefficients for the second Euler angle 
 
                     the coefficients for the third Euler angle 
 
                     the coefficients for the first Euler angle for 
                     the second logical record, ... 
 
                     and so on. 
 
   btime          Begin time (seconds past J2000 TDB) of first set 
                  of Chebyshev polynomial coefficients (first 
                  logical record). 
 
-Detailed_Output
 
    None. 
 
-Parameters
 
    None. 
 
-Exceptions
 
   1) If the number of sets of coefficients is not positive 
      SPICE(NUMCOEFFSNOTPOS) is signaled. 
 
   2) If the interval length is not positive, SPICE(INTLENNOTPOS) 
      is signaled. 
 
   3) If the integer code for the reference frame is not recognized, 
      SPICE(INVALIDREFFRAME) is signaled. 
 
   4) If segment stop time is not greater then the begin time, 
       SPICE(BADDESCRTIMES) is signaled. 
 
   5) If the time of the first record is not greater than 
      or equal to the descriptor begin time, SPICE(BADDESCRTIMES) 
      is signaled. 
 
   6) If the end time of the last record is not greater than 
      or equal to the descriptor end time, SPICE(BADDESCRTIMES) is 
      signaled. 
 
   7) If any input string pointers are null, the error 
      SPICE(NULLPOINTER) will be signaled.
      
   8) If any input strings have length zero, the error 
      SPICE(EMPTYSTRING) will be signaled.
 

-Files
 
   A new type 2 PCK segment is written to the PCK file attached 
   to HANDLE. 
 
-Particulars
 
   This routine writes an PCK type 2 data segment to the designated 
   PCK file, according to the format described in the PCK Required 
   Reading. 
 
   Each segment can contain data for only one body-fixed frame and base
   reference frame. The Chebyshev polynomial degree and length of time
   covered by each logical record are also fixed.  However, an
   arbitrary number of logical records of Chebyshev polynomial
   coefficients can be written in each segment. Minimizing the number
   of segments in a PCK file will help optimize how the SPICE system
   accesses the file.

 
-Examples
 
 
   Suppose that you have sets of Chebyshev polynomial coefficients 
   in an array `cdata' pertaining to the position of the moon (NAIF ID 
   = 301) in the J2000 reference frame, and want to put these into a 
   type 2 segment in an existing PCK file. The following code could 
   be used to add one new type 2 segment. To add multiple segments, 
   put the call to pckw02_c in a loop. 
 
      /.
      First open the PCK file and get a handle for it. 
      ./
      dafopw_c ( pcknam, &handle );
 
      /.
      Create a segment identifier. 
      ./ 
      segid = "MY_SAMPLE_PCK_TYPE_2_SEGMENT";
 
      /.
      Write the segment. 
      ./
      pckw02_c ( handle, 301,    "j2000", 
                 first,  last,   segid,   intlen, 
                 n,      polydg, cdata,   btime  );
 
      /.
      Close the file. 
      ./ 
      pckcls_c ( handle );

 
-Restrictions
 
   None. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman      (JPL)
   K.S. Zukor        (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 16-DEC-2016 (NJB) (KSZ)

-Index_Entries
 
   write pck type_2 data segment 
 
-&
*/

{ /* Begin pckw02_c */
 

   /*
   Participate in error tracing.
   */
   chkin_c ( "pckw02_c" );


   /*
   Check the input string frame to make sure the pointer is non-null
   and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "pckw02_c", frame );

   /*
   Check the input string segid to make sure the pointer is
   non-null and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "pckw02_c", segid );


   /*
   Call the f2c'd routine. 
   */
   pckw02_ ( (integer     *) &handle,
             (integer     *) &clssid,
             (char        *) frame,
             (doublereal  *) &first,
             (doublereal  *) &last,
             (char        *) segid,
             (doublereal  *) &intlen,
             (integer     *) &n,
             (integer     *) &polydg,
             (doublereal  *) cdata,
             (doublereal  *) &btime,
             (ftnlen       ) strlen(frame),
             (ftnlen       ) strlen(segid)  );

   chkout_c ( "pckw02_c" );

} /* End pckw02_c */
