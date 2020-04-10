/*

-Procedure pckopn_c ( PCK, open new file )

-Abstract
 
   Create a new PCK file, returning the handle of the opened file. 
 
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
   #include "SpiceZmc.h"
   #include "SpiceZst.h"

   void pckopn_c ( ConstSpiceChar   * name,
                   ConstSpiceChar   * ifname,
                   SpiceInt           ncomch,
                   SpiceInt         * handle  ) 

/*

-Brief_I/O
 
   VARIABLE  I/O  DESCRIPTION 
   --------  ---  -------------------------------------------------- 
   name       I   The name of the PCK file to be opened. 
   ifname     I   The internal filename for the PCK. 
   ncomch     I   The number of characters to reserve for comments. 
   handle     O   The handle of the opened PCK file. 
 
-Detailed_Input
 
   name     The name of the PCK file to be created. 
 
   ifname   The internal filename for the PCK file that is being 
            created. The internal filename may be up to 60 characters 
            long. `ifname' may not contain non-printing characters;
            otherwise there are no restrictions on its contents. 
 
   ncomch   This is the space, measured in characters, to be 
            initially set aside for the comment area when a new PCK 
            file is opened. The amount of space actually set aside 
            may be greater than the amount requested, due to the 
            manner in which comment records are allocated in an PCK 
            file. However, the amount of space set aside for comments 
            will always be at least the amount that was requested. 
 
            The value of ncomch should be greater than or equal to 
            zero, i.e., 0 <= ncomch. A negative value, should one 
            occur, will be assumed to be zero. 
 
-Detailed_Output
 
   handle   The handle of the opened PCK file. If an error occurs 
            when opening the file, the value of this variable should 
            not be used, as it will not represent a valid handle. 
 
-Parameters
 
   None. 
 
-Exceptions
 
   1) If the value of `ncomch' is negative, a value of zero (0) will 
      be used for the number of comment characters to be set aside 
      for comments. 
 
   2) If an error occurs while attempting to open a CK file the 
      value of `handle' will not represent a valid file handle. 
 
   3) If any input string pointers are null, the error 
      SPICE(NULLPOINTER) will be signaled.
      
   4) If any input strings have length zero, the error 
      SPICE(EMPTYSTRING) will be signaled.
 
-Files
 
   See descriptions of `name' and `handle'.
 
-Particulars
 
   Open a new PCK file, reserving room for comments if requested. 
 
-Examples
 
   Suppose that you want to create a new PCK file called 'new.PCK' 
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
 
   open a new pck file 
 
-&
*/

{ /* Begin pckopn_c */

 
   /*
   Participate in error tracing.
   */
   chkin_c ( "pckopn_c" );

   /*
   Check the input string name to make sure the pointer is non-null 
   and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "pckopn_c", name );
   
   /*
   Check the input string ifname to make sure the pointer is 
   non-null and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "pckopn_c", ifname );


   /*
   Call the f2c'd Fortran routine.
   */
   pckopn_ ( (char       *) name,
             (char       *) ifname,
             (integer    *) &ncomch,
             (integer    *) handle,
             (ftnlen      ) strlen(name),
             (ftnlen      ) strlen(ifname) );


   chkout_c ( "pckopn_c" );

} /* End pckopn_c */
