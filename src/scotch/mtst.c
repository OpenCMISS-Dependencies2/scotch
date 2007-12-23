/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : mtst.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program gives statistics on source **/
/**                meshes.                                 **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 25 feb 2003     **/
/**                                 to     27 jan 2004     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 23 dec 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MTST

#include "common.h"
#include "scotch.h"
#include "mtst.h"

/*
**  The static definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "-", NULL, "r" },
                              { "-", NULL, "w" } };

static const char *         C_usageList[] = {
  "mtst [<input mesh file> [<output data file>]] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Mesh         meshdat;
  SCOTCH_Num          velmnbr;
  SCOTCH_Num          vnodnbr;
  SCOTCH_Num          vnlomin;
  SCOTCH_Num          vnlomax;
  SCOTCH_Num          vnlosum;
  double              vnloavg;
  double              vnlodlt;
  SCOTCH_Num          edegmin;
  SCOTCH_Num          edegmax;
  double              edegavg;
  double              edegdlt;
  SCOTCH_Num          ndegmin;
  SCOTCH_Num          ndegmax;
  double              ndegavg;
  double              ndegdlt;
  SCOTCH_Num          edgenbr;
  int                 i;

  errorProg ("mtst");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    exit       (0);
  }

  for (i = 0; i < C_FILENBR; i ++)                /* Set default stream pointers */
    C_fileTab[i].pntr = (C_fileTab[i].mode[0] == 'r') ? stdin : stdout;
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        C_fileTab[C_fileNum ++].name = argv[i];
      else {
        errorPrint ("main: too many file names given");
        exit       (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          exit       (0);
        case 'V' :
          fprintf (stderr, "mtst, version %s - F. Pellegrini\n", SCOTCH_VERSION);
          fprintf (stderr, "Copyright 2004,2007 ENSEIRB, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option (\"%s\")", argv[i]);
          exit       (1);
      }
    }
  }

  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      if ((C_fileTab[i].pntr = fopen (C_fileTab[i].name, C_fileTab[i].mode)) == NULL) { /* Open the file */
        errorPrint ("main: cannot open file (%d)", i);
        exit       (1);
      }
    }
  }

  SCOTCH_meshInit  (&meshdat);
  SCOTCH_meshLoad  (&meshdat, C_filepntrsrcinp, -1);
  SCOTCH_meshCheck (&meshdat);

  SCOTCH_meshSize  (&meshdat, &velmnbr, &vnodnbr, &edgenbr);
  SCOTCH_meshStat  (&meshdat, &vnlomin, &vnlomax, &vnlosum, &vnloavg, &vnlodlt,
                    &edegmin, &edegmax, &edegavg, &edegdlt,
                    &ndegmin, &ndegmax, &ndegavg, &ndegdlt);

  fprintf (C_filepntrdatout, "S\tElements\tnbr=%ld\n",
           (long) velmnbr);
  fprintf (C_filepntrdatout, "S\tNodes\tnbr=%ld\n",
           (long) vnodnbr);
  fprintf (C_filepntrdatout, "S\tNode load\tmin=%ld\tmax=%ld\tsum=%ld\tavg=%g\tdlt=%g\n",
           (long) vnlomin, (long) vnlomax, (long) vnlosum, vnloavg, vnlodlt);
  fprintf (C_filepntrdatout, "S\tElement degree\tmin=%ld\tmax=%ld\tsum=%ld\tavg=%g\tdlt=%g\n",
           (long) edegmin, (long) edegmax, (long) (edgenbr / 2), edegavg, edegdlt);
  fprintf (C_filepntrdatout, "S\tNode degree\tmin=%ld\tmax=%ld\tsum=%ld\tavg=%g\tdlt=%g\n",
           (long) ndegmin, (long) ndegmax, (long) (edgenbr / 2), ndegavg, ndegdlt);
  fprintf (C_filepntrdatout, "S\tEdge\tnbr=%ld\n",
           (long) (edgenbr / 2));

#ifdef SCOTCH_DEBUG_ALL
  SCOTCH_meshExit (&meshdat);

  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      fclose (C_fileTab[i].pntr);                 /* Close the stream */
    }
  }
#endif /* SCOTCH_DEBUG_ALL */

  exit (0);
}
