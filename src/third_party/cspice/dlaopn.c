/* dlaopn.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c_b8 = 1000000;
static integer c_n1 = -1;

/* $Procedure DLAOPN ( DLA, open new file ) */
/* Subroutine */ int dlaopn_(char *fname, char *ftype, char *ifname, integer *
	ncomch, integer *handle, ftnlen fname_len, ftnlen ftype_len, ftnlen 
	ifname_len)
{
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer ncomr;
    extern /* Subroutine */ int dasadi_(integer *, integer *, integer *), 
	    sigerr_(char *, ftnlen), dasonw_(char *, char *, char *, integer *
	    , integer *, ftnlen, ftnlen, ftnlen), chkout_(char *, ftnlen), 
	    setmsg_(char *, ftnlen), errint_(char *, integer *, ftnlen);
    extern logical return_(void);

/* $ Abstract */

/*     Open a new DLA file and set the file type. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     DAS */
/*     DLA */

/* $ Keywords */

/*     DAS */
/*     DLA */
/*     FILES */

/* $ Declarations */

/*     Include file dla.inc */

/*     This include file declares parameters for DLA format */
/*     version zero. */

/*        Version 3.0.1 17-OCT-2016 (NJB) */

/*           Corrected comment: VERIDX is now described as a DAS */
/*           integer address rather than a d.p. address. */

/*        Version 3.0.0 20-JUN-2006 (NJB) */

/*           Changed name of parameter DSCSIZ to DLADSZ. */

/*        Version 2.0.0 09-FEB-2005 (NJB) */

/*           Changed descriptor layout to make backward pointer */
/*           first element.  Updated DLA format version code to 1. */

/*           Added parameters for format version and number of bytes per */
/*           DAS comment record. */

/*        Version 1.0.0 28-JAN-2004 (NJB) */


/*     DAS integer address of DLA version code. */


/*     Linked list parameters */

/*     Logical arrays (aka "segments") in a DAS linked array (DLA) file */
/*     are organized as a doubly linked list.  Each logical array may */
/*     actually consist of character, double precision, and integer */
/*     components.  A component of a given data type occupies a */
/*     contiguous range of DAS addresses of that type.  Any or all */
/*     array components may be empty. */

/*     The segment descriptors in a SPICE DLA (DAS linked array) file */
/*     are connected by a doubly linked list.  Each node of the list is */
/*     represented by a pair of integers acting as forward and backward */
/*     pointers.  Each pointer pair occupies the first two integers of a */
/*     segment descriptor in DAS integer address space.  The DLA file */
/*     contains pointers to the first integers of both the first and */
/*     last segment descriptors. */

/*     At the DLA level of a file format implementation, there is */
/*     no knowledge of the data contents.  Hence segment descriptors */
/*     provide information only about file layout (in contrast with */
/*     the DAF system).  Metadata giving specifics of segment contents */
/*     are stored within the segments themselves in DLA-based file */
/*     formats. */


/*     Parameter declarations follow. */

/*     DAS integer addresses of first and last segment linked list */
/*     pointer pairs.  The contents of these pointers */
/*     are the DAS addresses of the first integers belonging */
/*     to the first and last link pairs, respectively. */

/*     The acronyms "LLB" and "LLE" denote "linked list begin" */
/*     and "linked list end" respectively. */


/*     Null pointer parameter. */


/*     Segment descriptor parameters */

/*     Each segment descriptor occupies a contiguous */
/*     range of DAS integer addresses. */

/*        The segment descriptor layout is: */

/*           +---------------+ */
/*           | BACKWARD PTR  | Linked list backward pointer */
/*           +---------------+ */
/*           | FORWARD PTR   | Linked list forward pointer */
/*           +---------------+ */
/*           | BASE INT ADDR | Base DAS integer address */
/*           +---------------+ */
/*           | INT COMP SIZE | Size of integer segment component */
/*           +---------------+ */
/*           | BASE DP ADDR  | Base DAS d.p. address */
/*           +---------------+ */
/*           | DP COMP SIZE  | Size of d.p. segment component */
/*           +---------------+ */
/*           | BASE CHR ADDR | Base DAS character address */
/*           +---------------+ */
/*           | CHR COMP SIZE | Size of character segment component */
/*           +---------------+ */

/*     Parameters defining offsets for segment descriptor elements */
/*     follow. */


/*     Descriptor size: */


/*     Other DLA parameters: */


/*     DLA format version.  (This number is expected to occur very */
/*     rarely at integer address VERIDX in uninitialized DLA files.) */


/*     Characters per DAS comment record. */


/*     End of include file dla.inc */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     FNAME      I   Name of a DLA file to be opened. */
/*     FTYPE      I   Mnemonic code for type of data in the DLA file. */
/*     IFNAME     I   Internal file name. */
/*     NCOMR      I   Number of comment records to allocate. */
/*     HANDLE     O   Handle assigned to the opened DLA file. */

/* $ Detailed_Input */

/*     FNAME       is the name of a new DLA file to be created.  The */
/*                 file will be left opened for write access. */

/*     FTYPE       is a code for type of data placed into a DLA file. The */
/*                 non-blank part of FTYPE is used as the "file type" */
/*                 portion of the ID word in the DLA file. */

/*                 The first nonblank character and the three, or fewer, */
/*                 characters immediately following it, giving four */
/*                 characters, are used to represent the type of the */
/*                 data placed in the DLA file. This is provided as a */
/*                 convenience for higher level software. It is an error */
/*                 if this string is blank. Also, the file type may not */
/*                 contain any nonprinting characters. When written to */
/*                 the DLA file, the value for the type IS case */
/*                 sensitive. */

/*                 NAIF has reserved for its own use file types */
/*                 consisting of the upper case letters (A-Z) and the */
/*                 digits 0-9. NAIF recommends lower case or mixed case */
/*                 file types be used by all others in order to avoid */
/*                 any conflicts with NAIF file types. */

/*     IFNAME      is the internal file name for the new file.  The name */
/*                 may contain as many as 60 characters.  This name */
/*                 should uniquely identify the file. */

/*     NCOMR       is the number of comment records to allocate. */
/*                 Allocating comment records at file creation time may */
/*                 reduce the likelihood of having to expand the */
/*                 comment area later. */

/* $ Detailed_Output */

/*     HANDLE      is the file handle associated with the file. This */
/*                 handle is used to identify the file in subsequent */
/*                 calls to other DLA routines. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the input filename is blank, the error will be diagnosed by */
/*        routines in the call tree of this routine.  No file will be */
/*        created. */

/*     2) If the specified file cannot be opened without exceeding */
/*        the maximum allowed number of open DAS files, the error */
/*        will be diagnosed by routines in the call tree of this */
/*        routine.  No file will be created. */

/*     3) If the file cannot be opened properly, the error will be */
/*        diagnosed by routines in the call tree of this routine.  No */
/*        file will be created. */

/*     4) If the initial records in the file cannot be written, the */
/*        error is diagnosed by routines in the call tree of this */
/*        routine.  No file will be created. */

/*     5) If no logical units are available, the error will be diagnosed */
/*        by routines in the call tree of this routine. No file will be */
/*        created. */

/*     6) If the file type is blank, the error will be diagnosed by */
/*        routines in the call tree of this routine.  No file will be */
/*        created. */

/*     7) If the file type contains nonprinting characters, decimal */
/*        0-31 and 127-255, the error will be diagnosed by routines in */
/*        the call tree of this routine.  No file will be created. */

/*     8) If the number of comment records allocated NCOMR is negative, */
/*        the error SPICE(BADRECORDCOUNT) will be signaled. No file will */
/*        be created. */

/* $ Files */

/*     See argument FNAME. */

/* $ Particulars */

/*     DLA files are built using the DAS low-level format; DLA files are */
/*     a specialized type of DAS file in which data are organized as a */
/*     doubly linked list of segments. Each segment's data belong to */
/*     contiguous components of character, double precision, and integer */
/*     type. */

/*     This routine creates a new DLA file and sets the type of the */
/*     file to the mnemonic code passed to it. */

/*     DLA files created by this routine have initialized file records. */
/*     The ID word in a DLA file record has the form */

/*        DAS/+*** */

/*     where the characters following the slash are supplied by the */
/*     caller of this routine. */

/* $ Examples */

/*     1)  Create a new DLA file, using an internal file name that */
/*         attempts to serve as an unique identifier, and give the file a */
/*         type of 'TEST'.  No room for comments will be reserved. */

/*            FNAME  =  'TEST.DLA' */
/*            FTYPE  =  'TEST' */
/*            IFNAME =  'TEST.DLA/NAIF/NJB/07-FEB-2005-02:57:00' */
/*            NCOMCH =   0 */

/*            CALL DLAOPN ( FNAME, FTYPE, IFNAME, NCOMCH, HANDLE ) */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     K.R. Gehringer  (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) */

/*        Updated version info. */

/*        01-APR-2016 (NJB) */

/*           Changed short error message for invalid comment */
/*           count. Corrected reference to "DASCLU" in comments. */

/*        08-OCT-2009 (NJB) */

/*           Updated header. */

/*        09-FEB-2005 (NJB) (KRG) */

/* -& */
/* $ Index_Entries */

/*     open a new dla file */
/*     open a new dla file with write access */

/* -& */
/* $ Revisions */

/*     None. */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("DLAOPN", (ftnlen)6);

/*     Compute the number of comment records required. */

    if (*ncomch > 0) {
	ncomr = (*ncomch - 1) / 1024 + 1;
    } else if (*ncomch == 0) {
	ncomr = 0;
    } else {
	setmsg_("Requested number of comment characters must be non-negative"
		" but was #.", (ftnlen)70);
	errint_("#", ncomch, (ftnlen)1);
	sigerr_("SPICE(BADRECORDCOUNT)", (ftnlen)21);
	chkout_("DLAOPN", (ftnlen)6);
	return 0;
    }

/*     Let the DAS "open new" routine do the work. */

    dasonw_(fname, ftype, ifname, &ncomr, handle, fname_len, ftype_len, 
	    ifname_len);

/*     Write the format version. */

    dasadi_(handle, &c__1, &c_b8);

/*     Initialize the forward and backward segment list pointers. */

    dasadi_(handle, &c__1, &c_n1);
    dasadi_(handle, &c__1, &c_n1);

/*     We leave the file open, since further writes to the file */
/*     should occur next.  The file will eventually be closed */
/*     by a call to DASCLS or DASLLC, if all goes well. */

    chkout_("DLAOPN", (ftnlen)6);
    return 0;
} /* dlaopn_ */

