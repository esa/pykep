/*

-Procedure dsksrf_c ( DSK, get surface IDs for body )

-Abstract
 
   Find the set of surface ID codes for all surfaces associated with 
   a given body in a specified DSK file. 
 
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
 
   CELLS 
   DAS 
   DSK 
   NAIF_IDS 
   SETS 
    
-Keywords
 
   COVERAGE 
   SURFACE 
   TOPOGRAPHY 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZst.h"
   #include "SpiceZmc.h"
   #undef dsksrf_c


   void dsksrf_c ( ConstSpiceChar  * dsk,
                   SpiceInt          bodyid,
                   SpiceCell       * srfids ) 
/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   dsk        I   Name of DSK file. 
   bodyid     I   Integer body ID code. 
   srfids    I-O  Set of ID codes of surfaces in DSK file. 
 
-Detailed_Input
 
   dsk            is the name of a DSK file. This file will be  
                  opened for read access by this routine. 
    
 
   bodyid         is the integer ID code of a body for which 
                  topographic data are present in the specified DSK 
                  file. 
 
   srfids         is an initialized CSPICE set data structure. 
 
                  `srfids' optionally may contain a set of surface ID 
                  codes on input; on output, the ID codes already 
                  present in `srfids' will be combined with surface ID 
                  code set found for the body designated by 
                  `bodyid' in the file DSK. 
 
                  If `srfids' contains no data on input, its size and 
                  cardinality still must be initialized. 
 
-Detailed_Output
 
   srfids         is a CSPICE set data structure that contains 
                  the union of its contents upon input with the set 
                  of ID codes of the surfaces associated with the 
                  body designated by `bodyid', for which segments were 
                  found in the indicated DSK file. 
 
                  The elements of CSPICE sets are unique; each ID 
                  code in `srfids' appears only once, even if the DSK 
                  file contains multiple segments for that ID code. 
 
                  See the Examples section below for a complete 
                  example program showing how to retrieve body and 
                  surface ID codes from a DSK file. 
 
-Parameters
 
   None. 
 
-Exceptions
 
   1)  If the input file has transfer format, the error  
       SPICE(INVALIDFORMAT) is signaled. 
 
   2)  If the input file is not a transfer file but has architecture 
       other than DAS, the error SPICE(BADARCHTYPE) is signaled. 
 
   3)  If the input file is a binary DAS file of type other than 
       DSK, the error SPICE(BADFILETYPE) is signaled. 
 
   4)  If the DSK file cannot be opened or read, the error will 
       be diagnosed by routines called by this routine. 
 
   5)  If the size of the output set argument `srfids' is insufficient 
       to contain the actual number of ID codes of objects covered 
       by the indicated DSK file, the error will be diagnosed by 
       routines called by this routine. 
 
-Files
 
   See the description of `dsk' above.
 
-Particulars
 
   This routine provides an API via which applications can determine
   the set of surfaces associated with a given body in a specified
   DSK file. This routine is normally used together with dskobj_c.

-Examples

   The formatting of the results shown for this example may differ 
   across platforms. 
 
 
   1)  Display the coverage for each object and surface in a specified
       DSK file. Find the set of objects in the file. Loop over the
       contents of the ID code set: find the surface ID for each item
       in the set and display the surface ID.
 
 
   Example code begins here. 
 

      #include <stdio.h>
      #include "SpiceUsr.h"

      /.
      Program EX1 

      Examine a DSK file and identify the set of
      central bodies associated with the segments
      in the file. For each body, find the
      set of surfaces associated with that body.
      ./
      int main()
      {
         /.
         Local constants 
         ./
         #define FILSIZ          256
         #define MAXID           10000

         /.
         Local variables 
         ./
         SPICEINT_CELL         ( bodids, MAXID );
         SPICEINT_CELL         ( srfids, MAXID );

         SpiceChar               dsk [ FILSIZ ];

         SpiceInt                i;
         SpiceInt                j;

         /.
         Prompt for the name of a DSK file.
         ./
         prompt_c ( "Enter name of DSK file > ", FILSIZ, dsk );

         dskobj_c ( dsk, &bodids );

         for ( i = 0;  i < card_c( &bodids );  i++ )
         {
            printf ( "\n"
                     "Body ID:     %d\n",  
                     (int)SPICE_CELL_ELEM_I( &bodids, i) );

            /.
            Get the surface IDs for the Ith body.
            ./
            dsksrf_c ( dsk, SPICE_CELL_ELEM_I( &bodids, i), &srfids );

            for ( j = 0;  j < card_c( &srfids );  j++ )
            {
               printf ( "   Surface ID:  %d\n",  
                        (int)SPICE_CELL_ELEM_I( &srfids, j) );   
            }
         }    
         return ( 0 );
      } 

 
   When this program was executed on a PC/Linux/gcc/64-bit 
   platform, using as input the name of a DSK created by a code 
   example in the header of dskw02_c, the output was: 

      Enter name of DSK file > ../phobos_3_3_3seg.bds

      Body ID:     401
         Surface ID:  1
         Surface ID:  2
         Surface ID:  3


-Restrictions
 
   1) If an error occurs while this routine is updating the set 
      `srfids', the set may be corrupted. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman   (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 28-JAN-2016 (NJB)

-Index_Entries
 
   find id codes of surfaces in dsk file 
 
-&
*/

{ /* Begin dsksrf_c */

   /*
   Participate in error tracing.
   */
   chkin_c ( "dsksrf_c" );

   /*
   Check the input string to make sure the pointer is non-null and
   the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "dsksrf_c", dsk );

   /*
   Make sure cell data type is integer.
   */
   CELLTYPECHK ( CHK_STANDARD, "dsksrf_c", SPICE_INT, srfids );

   /*
   Initialize the cell if necessary. 
   */
   CELLINIT ( srfids );

   /*
   The f2c'd routine will operate on the integer array associated
   with the input cell. 
   */
   dsksrf_ ( (char     *) dsk,
             (integer  *) &bodyid,
             (integer  *) srfids->base,
             (ftnlen    ) strlen(dsk)   );

   /*
   Sync the cell. 
   */      
   if ( !failed_c() )
   {
      zzsynccl_c ( F2C, srfids );
   }

   chkout_c ( "dsksrf_c" );

} /* End dsksrf_c */
