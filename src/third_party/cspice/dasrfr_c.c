/*

-Procedure dasrfr_c ( DAS, read file record )

-Abstract
 
   Return the contents of the file record of a specified DAS 
   file. 
 
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
 
-Keywords
 
   DAS 
   FILES 
   UTILITY 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZst.h"
   #include "SpiceZmc.h"
   #undef dasrfr_c


   void dasrfr_c ( SpiceInt            handle,
                   SpiceInt            idwlen,
                   SpiceInt            ifnlen,
                   SpiceChar         * idword,
                   SpiceChar         * ifname,
                   SpiceInt          * nresvr,
                   SpiceInt          * nresvc,
                   SpiceInt          * ncomr,
                   SpiceInt          * ncomc  )

/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   handle     I   DAS file handle. 
   idwlen     I   Length of ID word string.
   ifnlen     I   Length of internal file name string.
   idword     O   ID word. 
   ifname     O   DAS internal file name. 
   nresvr     O   Number of reserved records in file. 
   nresvc     O   Number of characters in use in reserved rec. area. 
   ncomr      O   Number of comment records in file. 
   ncomc      O   Number of characters in use in comment area. 
 
-Detailed_Input
 
   handle         is a file handle for a previously opened DAS file. 
 
   idwlen         is the number of characters available in the output
                  ID word string. Space for the null terminator is
                  included in this count.

   ifnlen         is the number of characters available in the output
                  internal file name string. Space for the null
                  terminator is included in this count.

-Detailed_Output
 
   idword      is the `ID word' contained in the first eight 
               characters of the file record. 
 
   ifname      is the internal file name of the DAS file.  The 
               maximum length of the internal file name is 60 
               characters. 
 
   nresvr      is the number of reserved records in the DAS file 
               specified by `handle'. 
 
   nresvc      is the number of characters in use in the reserved 
               record area of the DAS file specified by `handle'. 
 
   ncomr       is the number of comment records in the DAS file 
               specified by `handle'. 
 
   ncomc       is the number of characters in use in the comment area 
               of the DAS file specified by `handle'. 
 
 
-Parameters
 
   None. 
 
-Exceptions
 
   1) If the file read attempted by this routine fails, the error 
      SPICE(DASFILEREADFAILED) will be signaled. 
 
-Files
 
   See the description of `handle' under $Detailed_Input. 
 
-Particulars
 
   This routine provides a convenient way of retrieving the 
   information contained in the file record of a DAS file. 
 
-Examples
 
   1)  Obtain the internal file name, comment record count, and comment 
       character count of an existing DAS file. 
 
          #include <stdio.h>
          #include "SpiceUsr.h"

          int main()
          {
             /.
             Local constants 
             ./
             #define IDWLEN         9
             #define IFNLEN         61
             #define FILSIZ         256

             /.
             Local variables
             ./
             SpiceChar               fname [FILSIZ];             
             SpiceChar               idword[IDWLEN];
             SpiceChar               ifname[IFNLEN];

             SpiceInt                handle;
             SpiceInt                nresvr;
             SpiceInt                nresvc;
             SpiceInt                ncomr;
             SpiceInt                ncomc;
      

             /.
             Obtain the file name.
             ./
             prompt_c ( "Enter DAS file name > ", FILSIZ, fname );

             /.          
             Open the file for reading. 
             ./
             dasopr_c ( fname, &handle  );

             /.
             Retrieve the internal file name and print it. 
             ./ 
 
             dasrfr_c ( handle, IDWLEN, IFNLEN, idword,
                        ifname, &nresvr, &nresvc, &ncomr, &ncomc ); 

             printf ( "Internal file name is:           %s\n", ifname     );
             printf ( "Number of comment records is:    %d\n", (int)ncomr );
             printf ( "Number of comment characters is: %d\n", (int)ncomc );

             return(0);
          }

 
-Restrictions
 
   None. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman   (JPL) 
   K.R. Gehringer (JPL) 
   W.L. Taber     (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 11-NOV-2016 (NJB) (KRG) (WLT)

-Index_Entries
 
   read DAS file record 
   read DAS internal file name 
 
-&
*/

{ /* Begin dasrfr_c */


   /*
   Participate in error tracing.
   */
   chkin_c ( "dasrfr_c" );

   /*
   Check output string pointers and lengths. 
   */
   CHKOSTR ( CHK_STANDARD, "dasrfr_c", idword, idwlen );
   CHKOSTR ( CHK_STANDARD, "dasrfr_c", ifname, ifnlen );

   
   dasrfr_ ( (integer    *) &handle,
             (char       *) idword, 
             (char       *) ifname, 
             (integer    *) nresvr,
             (integer    *) nresvc,
             (integer    *) ncomr,
             (integer    *) ncomc,
             (ftnlen      ) idwlen,
             (ftnlen      ) ifnlen  );

   /*
   Convert the output strings to C style strings. 
   */
   F2C_ConvertStr ( idwlen, idword );
   F2C_ConvertStr ( ifnlen, ifname );


   chkout_c ( "dasrfr_c" );

} /* End dasrfr_c */
