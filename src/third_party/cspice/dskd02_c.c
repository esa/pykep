/*

-Procedure dskd02_c ( DSK, fetch d.p. type 2 data )

-Abstract
 
   Fetch double precision data from a type 2 DSK segment. 
 
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

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #undef dskd02_c


   void dskd02_c ( SpiceInt               handle,
                   ConstSpiceDLADescr   * dladsc,
                   SpiceInt               item,
                   SpiceInt               start,
                   SpiceInt               room,
                   SpiceInt             * n,
                   SpiceDouble          * values ) 
/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   handle     I   DSK file handle. 
   dladsc     I   DLA descriptor. 
   item       I   Keyword identifying item to fetch. 
   start      I   Start index. 
   room       I   Amount of room in output array. 
   n          O   Number of values returned. 
   values     O   Array containing requested item. 
 
-Detailed_Input
 
   handle         is the handle of a DSK file containing a type 2 
                  segment from which data are to be fetched. 
 
   dladsc         is the DLA descriptor associated with the segment 
                  from which data are to be fetched. 
 
   item           is an integer "keyword" parameter designating the
                  double precision data item to fetch.
 
                  Names, meanings, and value of keyword parameters
                  supported by this routine are given in the header
                  file

                     SpiceDSK.h
 
                  The keyword parameters for double precision data
                  listed there are supported by this routine.
 
   start          is the start index within specified data item from
                  which data are to be fetched.  The index of the first
                  element of each data item is 0. This convention
                  applies uniformly to all data, even if the data are
                  associated with a set of 1-based indices.  For
                  example, the vertex ID range starts at 1 (this fact
                  is language-independent), but a caller would use a
                  `start' value of 0 to fetch the first vertex.
 
   room           is the amount of room in the output array.  It is
                  permissible to provide an output array that has too
                  little room to fetch an item in one call.
 
 
-Detailed_Output
 
   n              is the number of elements fetched to the output 
                  array `values'. `n' is normally in the range  
                  1:room; if an error occurs on the call, `n' is 
                  undefined. 
 
   values         is a contiguous set of elements of the item 
                  designated by `item'. The correspondence of  
                  `values' with the elements of the data item is: 
 
                     values[0]      item[start] 
                       ...             ... 
                     values[n-1]    item[start+n-1]
                   
                  If an error occurs on the call, `values' is  
                  undefined. 
 
-Parameters
 
   See the header files  
 
      SpiceDLA.h 
      SpiceDSK.h 
 
-Exceptions
 
   1) If the input handle is invalid, the error will be diagnosed by 
      routines in the call tree of this routine.  
 
   2) If a file read error occurs, the error will be diagnosed by 
      routines in the call tree of this routine. 
 
   3) If the input DLA descriptor is invalid, the effect of this 
      routine is undefined. The error *may* be diagnosed by routines 
      in the call tree of this routine, but there are no 
      guarantees. 
 
   4) If `room' is non-positive, the error SPICE(VALUEOUTOFRANGE) 
      is signaled. 
 
   5) If the coarse voxel scale read from the designated segment 
      is less than 1, the error PICE(VALUEOUTOFRANGE) is signaled. 
 
   6) If the input keyword parameter is not recognized, the error 
      SPICE(NOTSUPPORTED) is signaled. 
 
   7) If `start' is less than 0 or greater than or equal to the size of
      the item to be fetched, the error SPICE(INDEXOUTOFRANGE) is
      signaled.
 
-Files
 
   See input argument `handle'. 
 
-Particulars
 
   Most SPICE applications will not need to call this routine. The
   routines dskv02_c, dskp02_c, and dskz02_c provide a higher-level
   interface for fetching DSK type 2 vertex and plate data.

   DSK files are built using the DLA low-level format and 
   the DAS architecture; DLA files are a specialized type of DAS 
   file in which data are organized as a doubly linked list of 
   segments.  Each segment's data belong to contiguous components of 
   character, double precision, and integer type. 

   Note that the DSK descriptor for the segment is not needed by this
   routine; the DLA descriptor contains the base address and size
   information for the integer, double precision, and character
   components of the segment, and these suffice for the purpose of
   fetching data.

-Examples
 

   The numerical results shown for this example may differ across
   platforms. The results depend on the SPICE kernels used as
   input, the compiler and supporting libraries, and the machine
   specific arithmetic implementation.


   1) Look up all the vertices associated with each plate 
      of the model contained in a specified type 2 segment. 
      For this example, we'll show the context of this look-up: 
      opening the DSK file for read access, traversing a trivial, 
      one-segment list to obtain the segment of interest. 


   Example code begins here.


         #include <stdio.h>
         #include "SpiceUsr.h"
         #include "SpiceDLA.h"
         #include "SpiceDSK.h"

         int main()
         {

            /.
            Local parameters
            ./
            #define FILSIZ          256 

            /.
            Local variables
            ./
            SpiceBoolean            found;

            SpiceChar               dsk     [ FILSIZ ];

            SpiceDLADescr           dladsc;

            SpiceDouble             vrtces  [3][3];

            SpiceInt                handle;
            SpiceInt                i;
            SpiceInt                j;
            SpiceInt                n;
            SpiceInt                np;
            SpiceInt                start;
            SpiceInt                vrtids  [3];


            /. 
            Prompt for the name of the DSK to read.
            ./
            prompt_c ( "Enter DSK name > ", FILSIZ, dsk );

            /.
            Open the DSK file for read access. We use the DAS-level
            interface for this function.
            ./
            dasopr_c ( dsk, &handle );

            /.
            Begin a forward search through the kernel, treating the
            file as a DLA. In this example, it's a very short search.
            ./
            dlabfs_c ( handle, &dladsc, &found );

            if ( !found  )
            {
               /.
               We arrive here only if the kernel
               contains no segments.  This is
               unexpected, but we're prepared for it.
               ./
               setmsg_c ( "No segments found in DSK file #." );
               errch_c  ( "#", dsk                           );
               sigerr_c ( "SPICE(NODATA)"                    );
            }

            /.
            If we made it this far, `dladsc' is the
            DLA descriptor of the first segment.

            Find the number of plates in the model.
            ./    
            dski02_c ( handle, &dladsc, SPICE_DSK02_KWNP, 
                       0,      1,       &n,             &np );

            /.
            For each plate, look up the desired data.
            ./
            for ( i = 1;  i <= np;  i++ )
            {
               /.
               For the Ith plate, find the associated
               vertex IDs.  We must take into account
               the fact that each plate has three
               vertices when we compute the start
               index.     
               ./
               start = 3*(i-1);

               dski02_c ( handle, &dladsc, SPICE_DSK02_KWPLAT, start,  
                          3,      &n,      vrtids                   );

               for ( j = 0;  j < 3;  j++  )
               {
                  /.
                  Fetch the vertex associated with
                  the jth vertex ID.  Again, each
                  vertex is a 3-vector.  Note that
                  the vertices are double-precision
                  data, so we fetch them using
                  dskd02_c.          
                  ./
                  start = (vrtids[j]-1)*3;

                  dskd02_c ( handle, &dladsc, SPICE_DSK02_KWVERT, start,  
                             3,      &n,      vrtces[j]               );
               }

               /.
               Display the vertices of the ith plate:
               ./
               printf ( "\n"
                        "Plate number:  %d\n"
                        "   Vertex 1: ( %+e   %+e   %+e )\n"
                        "   Vertex 2: ( %+e   %+e   %+e )\n"
                        "   Vertex 3: ( %+e   %+e   %+e )\n",
                        (int)i,
                        vrtces[0][0], vrtces[0][1], vrtces[0][2], 
                        vrtces[1][0], vrtces[1][1], vrtces[1][2], 
                        vrtces[2][0], vrtces[2][1], vrtces[2][2]  );
            }

            /.
            Close the kernel.  This isn't necessary in a stand-
            alone program, but it's good practice in subroutines
            because it frees program and system resources.
            ./
            dascls_c ( handle );

            return ( 0 );
         }


   When this program was run on a PC/Linux/gcc/64-bit platform,
   using a DSK file containing data representing a regular
   icosahedron, the output was:


      Plate number:  1
         Vertex 1: ( +0.000000e+00   +0.000000e+00   +1.175570e+00 )
         Vertex 2: ( +1.051460e+00   +0.000000e+00   +5.257310e-01 )
         Vertex 3: ( +3.249200e-01   +1.000000e+00   +5.257310e-01 )

      Plate number:  2
         Vertex 1: ( +0.000000e+00   +0.000000e+00   +1.175570e+00 )
         Vertex 2: ( +3.249200e-01   +1.000000e+00   +5.257310e-01 )
         Vertex 3: ( -8.506510e-01   +6.180340e-01   +5.257310e-01 )

         ...

      Plate number:  20
         Vertex 1: ( +8.506510e-01   -6.180340e-01   -5.257310e-01 )
         Vertex 2: ( +0.000000e+00   +0.000000e+00   -1.175570e+00 )
         Vertex 3: ( +8.506510e-01   +6.180340e-01   -5.257310e-01 )
 
 
-Restrictions

   1) The underlying f2c'd routine

         dskd02_ 

      called by this routine uses discovery check-in to boost execution
      speed.  However, that routine is in violation of NAIF standards
      for use of discovery check-in:  routines called from that routine
      may signal errors.  If errors are signaled in routines called
      from dskd02_, that routine's name will be missing from the
      traceback message.
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman    (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 04-APR-2017 (NJB)

      Updated parameter references in example program.
      Removed unnecessary include statements.
      Updated header.

      DSKLIB_C Version 1.0.0, 12-FEB-2010 (NJB)

-Index_Entries
 
   fetch d.p. data from a type 2 dsk segment 
 
-&
*/

{ /* Begin dskd02_c */

   /*
   Local variables
   */
   SpiceInt                fDLADescr  [ SPICE_DLA_DSCSIZ ];

   integer                 fItem;
   integer                 fRoom;
   integer                 fStart;


   /*
   Participate in error tracing.
   */
   chkin_c ( "dskd02_c" );

   /*
   Populate the Fortran DLA descriptor array fCurrent with the contents
   of the input descriptor.
   */
   fDLADescr[SPICE_DLA_BWDIDX] = dladsc->bwdptr;
   fDLADescr[SPICE_DLA_FWDIDX] = dladsc->fwdptr;
   fDLADescr[SPICE_DLA_IBSIDX] = dladsc->ibase;
   fDLADescr[SPICE_DLA_ISZIDX] = dladsc->isize;
   fDLADescr[SPICE_DLA_DBSIDX] = dladsc->dbase;
   fDLADescr[SPICE_DLA_DSZIDX] = dladsc->dsize;
   fDLADescr[SPICE_DLA_CBSIDX] = dladsc->cbase;
   fDLADescr[SPICE_DLA_CSZIDX] = dladsc->csize;

   /*
   Store the requested item in an integer variable. 
   */
   fItem  = item;

   /*
   Adjust the start value:  convert to Fortran-style indexing. 
   */
   fStart = start + 1;

   /*
   Store the room value in an integer variable. 
   */
   fRoom  = room;
   
   /*
   Call the f2c'd routine. 
   */
   dskd02_ ( ( integer      * ) &handle,      
             ( integer      * ) fDLADescr,
             ( integer      * ) &fItem,
             ( integer      * ) &fStart,
             ( integer      * ) &fRoom,
             ( integer      * ) n,
             ( doublereal   * ) values     );


   chkout_c ( "dskd02_c" );

} /* End dskd02_c */
