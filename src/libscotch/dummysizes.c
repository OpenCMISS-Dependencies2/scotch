/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id$
*/
/************************************************************/
/**                                                        **/
/**   NAME       : dummysizes.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of the libScotch compilation job.  **/
/**                This small program processes files that **/
/**                are in fact pattern header files for    **/
/**                the libScotch library, and replaces     **/
/**                symbolic sizes of the opaque libScotch  **/
/**                by the proper integer values according  **/
/**                to the machine on which it is run.      **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 22 oct 2001     **/
/**                                 to   : 22 nov 2001     **/
/**                # Version 4.0  : from : 25 nov 2001     **/
/**                                 to   : 06 jan 2006     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DUMMYSIZES

#define CHARMAX                     2048          /* Maximum line size */

#define SUBSMAX                     32            /* Maximum number of substitutions */

#define C_FILENBR                   2             /* Number of files in list                */
#define C_FILEARGNBR                2             /* Number of files which can be arguments */

#define C_filenamehedinp            C_fileTab[0].name /* Source graph input file name */
#define C_filenamehedout            C_fileTab[1].name /* Statistics output file name  */

#define C_filepntrhedinp            C_fileTab[0].pntr /* Source graph input file */
#define C_filepntrhedout            C_fileTab[1].pntr /* Statistics output file  */

#define EXPAND(s) EXPANDTWO(s)
#define EXPANDTWO(s) #s

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "geom.h"
#include "mesh.h"
#include "arch.h"
#include "mapping.h"
#include "order.h"
#include "library_mapping.h"
#include "library_order.h"

/*
**  The static definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "-", NULL, "r" },
                              { "-", NULL, "w" } };

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

void
subsFill (
char *                      substab[2],
char *                      origptr,
int                         subsval)
{
  char *              subsptr;

  subsptr = malloc (32 * sizeof (char));
  sprintf (subsptr, "%d", (int) (subsval + sizeof (double) - 1) / sizeof (double));

  substab[0] = origptr;
  substab[1] = subsptr;
}

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
  char                chartab[CHARMAX];
  char                chartmp[CHARMAX];
  char *              substab[SUBSMAX][2];        /* Substitution array */
  int                 subsnbr;
  int                 i;

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    printf ("Usage is:\ndummysizes [<input pattern header file> [<output header file>]]\n");
    return (((argv[1][0] == '?') && argv[1][1] == '\0') ? 0 : 1);
  }

  for (i = 0; i < C_FILENBR; i ++)                /* Set default stream pointers */
    C_fileTab[i].pntr = (C_fileTab[i].mode[0] == 'r') ? stdin : stdout;
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        C_fileTab[C_fileNum ++].name = argv[i];
      else {
        fprintf (stderr, "dummysizes: ERROR: main: too many file names given");
        exit    (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          printf ("Usage is:\ndummysizes [<input pattern header file> [<output header file>]]\n");
          exit       (0);
        case 'V' :
          fprintf (stderr, "dummysizes, version 4.0 - F. Pellegrini\n");
          fprintf (stderr, "Copyright 2004 INRIA, France\n");
          fprintf (stderr, "This software is libre/free software under LGPL -- see the user's manual for more information\n");
          return  (0);
        default :
          fprintf (stderr, "dummysizes: ERROR: main: unprocessed option (\"%s\")", argv[i]);
          exit    (1);
      }
    }
  }

  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      if ((C_fileTab[i].pntr = fopen (C_fileTab[i].name, C_fileTab[i].mode)) == NULL) { /* Open the file */
          fprintf (stderr, "dummysizes: ERROR: main: cannot open file (%d)", i);
          exit    (1);
      }
    }
  }

  substab[0][0] = "library.h";
  substab[0][1] = "scotch.h ";
  substab[1][0] = "libraryf.h";
  substab[1][1] = "scotchf.h ";
  substab[2][0] = "DUMMYINT";
  substab[2][1] = EXPAND(INT);
  substab[3][0] = "DUMMYMAXINT";
  substab[3][1] = EXPAND(INT_MAX);
  subsnbr = 4;
  subsFill (substab[subsnbr ++], "DUMMYSIZEARCH",  sizeof (Arch));
  subsFill (substab[subsnbr ++], "DUMMYSIZEGEOM",  sizeof (Geom));
  subsFill (substab[subsnbr ++], "DUMMYSIZEGRAPH", sizeof (Graph));
  subsFill (substab[subsnbr ++], "DUMMYSIZEMESH",  sizeof (Mesh));
  subsFill (substab[subsnbr ++], "DUMMYSIZEMAP",   sizeof (LibMapping));
  subsFill (substab[subsnbr ++], "DUMMYSIZEORDER", sizeof (LibOrder));
  subsFill (substab[subsnbr ++], "DUMMYSIZESTRAT", sizeof (Strat *));

  while (fgets (chartab, CHARMAX, C_filepntrhedinp) != NULL) { /* Infinite loop on file lines */
    int                 charnbr;
    int                 subsnum;

    if (((charnbr = strlen (chartab)) >= (CHARMAX - 1)) && /* If line read is at least as long as maximum size     */
        (chartab[CHARMAX - 1] != '\n')) {         /* And last character is not a newline, that is, some is missing */
      fprintf (stderr, "dummysizes: ERROR: line too long\n");
      exit    (1);
    }

    for (subsnum = 0; subsnum < subsnbr; subsnum ++) { /* Perform substitutions */
      char *              charptr;                /* Place where token found    */

      while ((charptr = strstr (chartab, substab[subsnum][0])) != NULL) { /* As long as substitution can be performed */
        int                 charnbr;
        int                 charnum;

        charnum = charptr - chartab;              /* Position where token found */
        charnbr = strlen (substab[subsnum][0]);   /* Length of token            */

        strcpy (chartmp, charptr + charnbr);      /* Save end of line */

        sprintf (charptr, "%s%s", substab[subsnum][1], chartmp); /* Replace end of line with substituted token */
      }
    }

    fputs (chartab, C_filepntrhedout);            /* Output possibly updated line */
  }

#ifdef SCOTCH_DEBUG_MAIN1
  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      fclose (C_fileTab[i].pntr);                 /* Close the stream */
    }
  }
#endif /* SCOTCH_DEBUG_MAIN1 */

  exit (0);
}
