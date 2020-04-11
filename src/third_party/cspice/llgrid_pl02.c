/*
 
-Procedure llgrid_pl02 ( Lon/lat grid using DSK type 2 plate model )
 
-Abstract
 
   Deprecated: This routine has been superseded by the CSPICE routine
   latsrf_c. This routine is supported for purposes of backward
   compatibility only.

   Given the planetocentric longitude and latitude values of a set of
   surface points on a specified target body, compute the corresponding
   rectangular coordinates of those points.  The target body's
   surface is represented by a triangular plate model contained in a
   type 2 DSK segment.

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
 
   FRAMES
   PCK
   SPK
   TIME
 
-Keywords
 
   GEOMETRY
 
*/
 
   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZmc.h"


   void llgrid_pl02 ( SpiceInt               handle,
                      ConstSpiceDLADescr   * dladsc,
                      SpiceInt               npoints,
                      ConstSpiceDouble       grid     [][2],
                      SpiceDouble            spoints  [][3],
                      SpiceInt               plateIDs []     )
/*
 
-Brief_I/O
 
   Variable  I/O  Description
   --------  ---  --------------------------------------------------
   handle     I   DSK handle.
   dlasdc     I   DLA descriptor of target body segment.
   npoints    I   Number of grid coordinate pairs.
   grid       I   Lon/lat values of surface points (radians).
   spoints    O   Rectangular coordinates of surface points.
   plateIDs   O   DSK plate IDs of surface points.
 
-Detailed_Input
 
   handle      is the DAS file handle of a DSK file open for read
               access.  This kernel must contain a type 2 segment
               that provides a plate model representing the entire
               surface of the target body.

   dladsc      is the DLA descriptor of a DSK segment representing
               the surface of a target body.
 
   npoints     is the number of longitude/latitude pairs in the array
               of grid points `grid'.

   grid        is an array of planetocentric longitude/latitude pairs
               to be mapped to surface points on the target body.
               `grid' should be declared by the caller

                   SpiceDouble grid [npoints][2];

               Elements

                  grid[i][0]
                  grid[i][1]

               are, respectively, the planetocentric longitude and
               latitude of the ith grid point.

               Units are radians.
 
-Detailed_Output

   spoints     is an array containing the rectangular (Cartesian)
               coordinates of the surface points on the target body,
               expressed relative to the body-fixed reference frame of
               the target body, corresponding to the input grid points.

               `spoints' should be declared by the caller

                   SpiceDouble spoints [npoints][3];

 
   plateIDs    is an array of integer ID codes of the plates on which
               the surface points are located.  The ith plate ID
               corresponds to the ith surface point. These ID codes can
               be use to look up data associated with the plate, such
               as the plate's vertices or outward normal vector.
 
               `plateIDs' should be declared by the caller

                  SpiceInt plateIDs [npoints];

-Parameters
 
   None.
 
-Exceptions
 
   If any of the listed errors occur, the output arguments are
   left unchanged.
 
 
   1) If a DSK providing a DSK type 2 plate model has not been
      loaded prior to calling llgrid_pl02, the error will be
      diagnosed and signaled by a routine in the call tree of this
      routine.

   2) If the segment associated with the input DLA descriptor is not
      of data type 2, the error SPICE(WRONGDATATYPE) is signaled.

   3) If a surface point cannot be computed because the ray corresponding
      to a longitude/latitude pair fails to intersect the target
      surface as defined by the plate model, the error
      SPICE(NOINTERCEPT) is signaled.


-Files
 
   The following data are required:
 
      - DSK data:  a DSK file containing a plate model representing the
        target body's surface must be loaded. This kernel must contain
        a type 2 segment that contains data for the entire surface of
        the target body.
 
   In all cases, kernel data are normally loaded once per program
   run, NOT every time this routine is called.
 
-Particulars
 
   See the headers of the CSPICE routines

      reclat_c
      latrec_c

   for detailed definitions of Planetocentric coordinates. 

-Examples
 
   The numerical results shown for this example may differ across
   platforms. The results depend on the SPICE kernels used as input,
   the compiler and supporting libraries, and the machine specific
   arithmetic implementation.
 
   In the following example program, the file

      phobos.3.3.bds

   is a DSK file containing a type 2 segment that provides a plate model 
   representation of the surface of Phobos. 

   Find the surface points on a target body corresponding to a given
   planetocentric longitude/latitude grid.  In order to duplicate
   the example output, the kernel name

      phobos.3.3.bds

   should be supplied at the prompt.


      #include <stdio.h>
      #include <math.h>
      #include "SpiceUsr.h"

      int main()
      {
         /.
         Local parameters
         ./
         #define  FILSIZ         256 
         #define  NAMLEN         33
         #define  NLAT           9
         #define  NLON           9 
         #define  MAXGRID        ( NLAT * NLON ) 
         #define  TOL            ( 1.e-12 )

         /.
         Local variables
         ./
         SpiceBoolean            found;

         SpiceChar               dsk      [ FILSIZ ];

         SpiceDLADescr           dladsc;

         SpiceDouble             grid     [ MAXGRID ][2];
         SpiceDouble             lat;
         SpiceDouble             lon;
         SpiceDouble             spoints  [ MAXGRID ][3];
         SpiceDouble             xlat;
         SpiceDouble             xlon;
         SpiceDouble             xr;

         SpiceInt                handle;
         SpiceInt                i;
         SpiceInt                j;
         SpiceInt                n;
         SpiceInt                npoints;
         SpiceInt                plateIDs [ MAXGRID ];



         /.
         Prompt for the name of the DSK to read and the name of
         the target body.
         ./
         prompt_c ( "Enter DSK name    > ", FILSIZ, dsk    );

         /.
         Open the DSK file for read access.
         We use the DAS-level interface for
         this function.
         ./
         dasopr_c ( dsk, &handle );

         /.
         Begin a forward search through the
         kernel, treating the file as a DLA.
         In this example, it's a very short
         search.
         ./
         dlabfs_c ( handle, &dladsc, &found );

         if ( !found  )
         {
            /.
            We arrive here only if the kernel
            contains no segments.  This is
            unexpected, but we're prepared for it.
            ./
            setmsg_c ( "No segments found in DSK file #.");
            errch_c  ( "#",  dsk                         );
            sigerr_c ( "SPICE(NODATA)"                   );
         }

         /.
         If we made it this far, DLADSC is the
         DLA descriptor of the first segment.

         Now generate the grid points.  We generate
         points along latitude bands, working from
         north to south.  The latitude range is selected
         to range from +80 to -80 degrees.  Longitude
         ranges from 0 to 350 degrees.  The increment
         is 20 degrees for latitude and 40 degrees for
         longitude.
         ./

         n = 0;

         for ( i = 0;  i < NLAT;  i++  )
         {
            lat = rpd_c() * ( 80.0 - 20.0*i );

            for ( j = 0;  j < NLON;  j++  )
            {
               lon = rpd_c() * 40.0*j;

               grid[n][0] = lon;
               grid[n][1] = lat;

               ++n;
            }
         }
         npoints = n - 1;

         /.
         Find the surface points corresponding to the grid points. 
         ./
         llgrid_pl02 ( handle,   
                       &dladsc, 
                       npoints,  
                       (ConstSpiceDouble (*)[2])grid,    
                       spoints, 
                       plateIDs                        );

         /.
         Print out the surface points in latitudinal
         coordinates and compare the derived lon/lat values
         to those of the input grid.
         ./
         for ( i = 0;  i < npoints;  i++  )
         {
            /.
            Use recrad_c rather than reclat_c to produce
            non-negative longitudes.
            ./
            recrad_c ( spoints[i], &xr, &xlon, &xlat );

            printf ( "\n"
                     "Intercept for grid point %d:\n"
                     "   Plate ID:                 %d\n"
                     "   Cartesian Coordinates:    (%e %e %e)\n"
                     "   Latitudinal Coordinates:\n"
                     "   Longitude (deg): %f\n"
                     "   Latitude  (deg): %f\n"
                     "   Radius     (km): %f\n"
                     "\n"
                     "Original grid coordinates:\n"
                     "   Longitude (deg): %f\n"
                     "   Latitude  (deg): %f\n"
                     "\n",
                     (int)i,
                     (int)plateIDs[i],
                     spoints[i][0], spoints[i][1], spoints[i][2],
                     xlon * dpr_c(),
                     xlat * dpr_c(),
                     xr,
                     grid[i][0] * dpr_c(),
                     grid[i][1] * dpr_c()                         );

            /.
            Perform sanity checks on the intercept
            coordinates.  Stop the program if any error
            is larger than our tolerance value.
            ./
            lon = grid[i][0];
            lat = grid[i][1];

            if ( fabs(xlat-lat) > TOL ) 
            {
               sigerr_c ( "Latitude error!" );
            }

            if (  (xlon - lon)  > pi_c()  )
            {
               xlon -= twopi_c();
            }

            if (  (xlon - lon)  > TOL  )
            {
               sigerr_c ( "Longitude error!" );
            }

         }

         /.
         Close the kernel.  This isn't necessary in a stand-
         alone program, but it's good practice in subroutines
         because it frees program and system resources.
         ./
         dascls_c ( handle );

         return ( 0 );
      }


   When this program was executed on a PC/Linux/gcc platform, 
   the output for the first 3 points (the rest of the output 
   is not shown due to its large volume) was:


      Enter DSK name    > phobos.3.3.bds

      Intercept for grid point 0:
         Plate ID:                 306238
         Cartesian Coordinates:    (1.520878e+00 -7.129713e-16 8.625327e+00)
         Latitudinal Coordinates:
         Longitude (deg): 360.000000
         Latitude  (deg): 80.000000
         Radius     (km): 8.758387

      Original grid coordinates:
         Longitude (deg): 0.000000
         Latitude  (deg): 80.000000


      Intercept for grid point 1:
         Plate ID:                 317112
         Cartesian Coordinates:    (1.189704e+00 9.982799e-01 8.807772e+00)
         Latitudinal Coordinates:
         Longitude (deg): 40.000000
         Latitude  (deg): 80.000000
         Radius     (km): 8.943646

      Original grid coordinates:
         Longitude (deg): 40.000000
         Latitude  (deg): 80.000000


      Intercept for grid point 2:
         Plate ID:                 324141
         Cartesian Coordinates:    (2.777752e-01 1.575341e+00 9.072029e+00)
         Latitudinal Coordinates:
         Longitude (deg): 80.000000
         Latitude  (deg): 80.000000
         Radius     (km): 9.211980

      Original grid coordinates:
         Longitude (deg): 80.000000
         Latitude  (deg): 80.000000

 


-Restrictions
 
   1) This routine assumes that the origin of the body-fixed reference 
      frame associated with the target body is located in the interior
      of that body.

   2) The results returned by this routine may not be meaningful
      if the target surface has multiple surface points associated
      with some (longitude, latitude) coordinates.
 
-Literature_References
 
   None.
 
-Author_and_Institution

   N.J. Bachman   (JPL)
 
-Version
 
   -CSPICE Version 2.0.0, 23-JUL-2016 (NJB) 
 
       Added failed_c calls.

       Include file references have been updated. Now calls
       zzdsksgr_ instead of dsksgr_. Integer output format
       in example program has been updated.

   -Beta Version 1.3.0, 30-APR-2014 (NJB) (BVS)

       Now includes dsk_proto.h.

       Last update was 1.2.1, 07-APR-2014 (BVS)

       Changed FRAME to FRAMES in the Required_Reading section.

   -Beta Version 1.2.0, 14-MAY-2010 (NJB)

       Updated for compatibility with new DSK design.
       DSK name was updated in example program.
 
   -Beta Version 1.1.0, 09-FEB-2007 (NJB)

       Bug fix:  type of local variable fDLADescr was changed to SpiceInt.

   -Beta Version 1.0.0, 06-NOV-2006 (NJB)
 
-Index_Entries
 
   map latitudinal grid to DSK type 2 plate model surface
 
-&
*/
 
{ /* Begin llgrid_pl02 */
 

   /*
   Prototypes 
   */


   /*
   Local parameters 
   */
   #define CTRIDX          1
   #define TYPIDX          3

   /*
   Local variables 
   */
   SpiceBoolean            found;
   
   SpiceDouble             fDSKDescr  [ SPICE_DSK_DSCSIZ ];
   SpiceDouble             maxrad;
   SpiceDouble             raydir     [3];
   SpiceDouble             scale;
   SpiceDouble             vertex     [3];

   SpiceInt                centerID;   
   SpiceInt                dataType;
   SpiceInt                fDLADescr  [ SPICE_DLA_DSCSIZ ];
   SpiceInt                i;

 
   /*
   Participate in error tracing.
   */
   if ( return_c()  )
   {
      return;
   }
   chkin_c ( "llgrid_pl02" );

   /*
   Fetch and examine the DSK descriptor of the segment from which we're
   getting the shape data.  Make sure this is a type 2 segment for the
   specified target.  Since we don't know a priori that the segment
   has type 2, we can't use the type 2 fetch routine for this operation.
   Instead, we use the f2c'd version of the Fortran routine dskgd.
   We must prepare a Fortran-style DLA descriptor in order to call
   this routine.

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
   Call the f2c'd "get DSK descriptor" routine.
   */
   dskgd_ ( ( integer    * ) &handle,
            ( integer    * ) fDLADescr,
            ( doublereal * ) fDSKDescr  );

   if ( failed_c() )
   {
      chkout_c ( "llgrid_pl02" );
      return;
   }

   /*
   Extract object ID and DSK data type from the descriptor. 
   */
   centerID = (SpiceInt)fDSKDescr[CTRIDX];
   dataType = (SpiceInt)fDSKDescr[TYPIDX];


   if ( dataType != 2 ) 
   {
      setmsg_c ( "Input segment has DSK data type #.  A segment of "
                 "type 2 is required."                               );
      errint_c ( "#", dataType                                       );
      sigerr_c ( "SPICE(WRONGDATATYPE)"                              );
      chkout_c ( "llgrid_pl02"                                       );
      return;
   }

   /*
   We're done with error checks on our inputs.

   Get the maximum radius value associated with the target body.  
   We'll use this later to compute a numerically safe ray vertex.
   */
   maxrad = zzdsksgr_ ( fDSKDescr );

   if ( failed_c() )
   {
      chkout_c ( "llgrid_pl02" );
      return;
   }
   
   for ( i = 0;  i < npoints;  i++  )
   {
      /*
      Find the outward unit vector corresponding to the ith lon/lat pair.
      Note:  longitude comes first!
      */
      latrec_c ( 1.0, grid[i][0], grid[i][1], vertex );

      vminus_c ( vertex, raydir );

      /*
      To avoid numerical problems, we pick a vertex that is guaranteed
      to be a reasonable distance away from the target's surface.
      */ 
      scale = maxd_c ( 2,  1.0,  2.0*maxrad );

      vscl_c ( scale, vertex, vertex );


      /*
      Find the surface intercept defined by the vertex, ray direction,
      and surface plate model.
      */
      dskx02_c ( handle,       dladsc,      vertex,  raydir, 
                 plateIDs+i,   spoints[i],  &found          );

      if ( failed_c() )
      {
         chkout_c ( "llgrid_pl02" );
         return;
      }

      if ( !found  )
      { 
         setmsg_c ( "Ray from vertex number # having longitude # "
                    "and latitude # (radians) to center "
                    "of target # did not intersect the surface  "
                    "defined by the input handle and descriptor."  );    
         errint_c ( "#", i                                         );
         errdp_c  ( "#", grid[i][0]                                );
         errdp_c  ( "#", grid[i][1]                                );
         errint_c ( "#", centerID                                  );
         sigerr_c ( "SPICE(NOINTERCEPT)"                           );
         chkout_c ( "llgrid_pl02"                                  );
         return;
      }
   }


   /*
   At this point, `spoints' and `plateIDs' are set.
   */

   chkout_c ( "llgrid_pl02" );
 
} /* End llgrid_pl02 */
