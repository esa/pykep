/*

-Procedure dskopn_c ( DSK, open new file )

-Abstract
 
   Open a new DSK file for subsequent write operations. 
 
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
 
   DAS 
   DSK 
 
-Keywords
 
   DAS 
   DSK 
   FILES 
 
*/

   #include <string.h>
   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZst.h"
   #include "SpiceZmc.h"
   #include "SpiceDSK.h"

   void dskopn_c ( ConstSpiceChar  * fname,
                   ConstSpiceChar  * ifname,
                   SpiceInt          ncomch,
                   SpiceInt       *  handle ) 

/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   fname      I   Name of a DSK file to be opened. 
   ifname     I   Internal file name. 
   ncomch     I   Number of comment characters to allocate. 
   handle     O   Handle assigned to the opened DSK file. 
 
-Detailed_Input
 
   fname       is the name of a new DSK file to be created.  The 
               file will be left opened for write access. 
 
   ifname      is the internal file name for the new file.  The name 
               may contain as many as 60 characters.  All characters 
               of `ifname' should be printing characters (ASCII codes 
               32-126 decimal). This name should uniquely identify 
               the file. 
 
   ncomch      is the number of comment characters to allocate. 
               Allocating comment characters at file creation time 
               may reduce the likelihood of having to expand the 
               comment area later. 
 
-Detailed_Output
 
   handle      is the file handle associated with the file. This 
               handle is used to identify the file in subsequent 
               calls to other DSK routines. 
 
-Parameters
 
   None. 
 
-Exceptions
 
   1) If the input filename is blank, the error will be diagnosed by 
      routines in the call tree of this routine.  No file will be 
      created. 
 
   2) If the specified file cannot be opened without exceeding 
      the maximum allowed number of open DAS files, the error 
      will be diagnosed by routines in the call tree of this  
      routine.  No file will be created. 
 
   3) If the file cannot be opened properly, the error will be 
      diagnosed by routines in the call tree of this routine.  No 
      file will be created. 
 
   4) If the initial records in the file cannot be written, the 
      error is diagnosed by routines in the call tree of this 
      routine.  No file will be created. 
 
   5) If no logical units are available, the error will be diagnosed 
      by routines in the call tree of this routine. No file will be 
      created. 
 
   6) If the internal file name contains nonprinting characters (ASCII 
      codes decimal 0-31 and 127-255), the error will be diagnosed 
      by routines in the call tree of this routine.  No file will be 
      created. 
 
   7) If the number of comment characters allocated NCOMCH is  
      negative, the error will be diagnosed by routines in the call 
      tree of this routine.  No file will be created. 
 
-Files
 
   See argument FNAME. 
 
-Particulars
 
   DSK files are built using the DLA low-level format and 
   the DAS architecture; DLA files are a specialized type of DAS 
   file in which data are organized as a doubly linked list of 
   segments.  Each segment's data belong to contiguous components of 
   character, double precision, and integer type. 
 
   This routine creates a new DSK file and sets the type of the 
   file to the mnemonic code passed to it. 
 
   DSK files created by this routine have initialized file records. 
   The ID word in a DSK file record has the form 
 
      DAS/DSK 
 
   where the characters following the slash are supplied by the 
   caller of this routine. 
 
-Examples
 
   1)  Create a new DSK file, using an internal file name that 
       attempts to serve as an unique identifier.  No room for 
       comments will be reserved. 
 
          fname  =  "TEST.DSK"; 
          ifname =  "TEST.DSK/NAIF/NJB/20-OCT-2006/14:37:00";
          ncomch =   0;
  
          dskopn_c ( fname, ifname, ncomch, &handle );
 
-Restrictions
 
   None. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman    (JPL) 
 
-Version
 
   -CSPICE Version 1.0.1, 23-JAN-2016 (NJB)

      Corrected spelling errors in comments.

   -DSKLIB_C Version 1.0.0, 12-FEB-2010 (NJB)

-Index_Entries
 
   open a new dsk file 
   open a new dsk file with write access 
 
-&
*/

{ /* Begin dskopn_c */

 


   /*
   Participate in error tracing.
   */

   chkin_c ( "dskopn_c" );

   /*
   Check the input string fname to make sure the pointer is non-null 
   and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "dskopn_c", fname );
   
   /*
   Check the input string ifname to make sure the pointer is 
   non-null and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "dskopn_c", ifname );


   dskopn_  (  ( char        * ) fname,
               ( char        * ) ifname,
               ( integer     * ) &ncomch,
               ( integer     * ) handle,
               ( ftnlen        ) strlen(fname),
               ( ftnlen        ) strlen(ifname)  );



   chkout_c ( "dskopn_c" );

} /* End dskopn_c */
