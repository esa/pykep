/*

-Procedure dskcls_c ( DSK, close file )

-Abstract
 
   Close a DSK file. 
 
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
   TOPOGRAPHY

*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZst.h"

   void dskcls_c ( SpiceInt      handle,
                   SpiceBoolean  optmiz ) 

/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   handle     I   Handle assigned to the opened DSK file. 
   optmiz     I   Flag indicating whether to segregate the DSK. 
 
-Detailed_Input
 
   handle      is the DAS file handle associated with the file. 
               The file may be open for read or write access. 
 
   optmiz      is a logical flag indicating whether the DSK  
               should be segregated before it is closed. This 
               option applies only to files open for write  
               access. The value of `optmiz' has no effect for 
               files opened for read access. 
 
               See the DAS Required Reading das.req for a  
               discussion of segregation of DAS files. 
 
-Detailed_Output
 
   None. This routine operates by side effects. 
  
-Parameters
 
   None. 
 
-Exceptions
 
   1) If an error occurs when the file is closed, the error will be  
      diagnosed by routines in the call tree of this routine. 
 
-Files
 
   See argument `handle'. 
 
-Particulars
 
   This routine provides a DSK-level interface for closing DSK files. 
 
   In cases where DSKs opened for write access are to be closed 
   without segregation, this interface is slightly simpler than that 
   available at the DAS level. 
 
-Examples
 
   1) Close a new DSK file using DAS segregation. `handle' 
      is the DAS file handle of the DSK.  
 
      This is the normal choice for DSK creation. 
  
         dskcls_c ( HANDLE, SPICETRUE ) 
 
   2) Close a new DSK file without using DAS segregation. The  
      close operation will be fast, but reading the file will be  
      less efficient than if the file had been segregated. 
 
         dskcls_c ( HANDLE, SPICETRUE ) 
 
   3) Close an existing DSK file that had been opened 
      for read access. In this case `optmiz' is ignored: 
  
         dskcls_c ( HANDLE, SPICEFALSE ) 
 
      or 
 
         dskcls_c ( HANDLE, SPICETRUE ) 
    
-Restrictions
 
   1) This routine should not be called by user applications 
      that have loaded a DSK file via furnsh_c. Such applications 
      should call the functions unload_c or kclear_c instead. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman    (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 23-JAN-2016 (NJB)

-Index_Entries
 
   close a dsk file 
 
-&
*/

{ /* Begin dskcls_c */

   /*
   Local constants
   */
   logical                 optflg;



   /*
   Participate in error tracing.
   */
   chkin_c ( "dskcls_c" );

   optflg = (logical) optmiz;

   dskcls_ ( (integer  *) &handle,
             (logical  *) &optflg  );

   chkout_c ( "dskcls_c" );

} /* End dskcls_c */
