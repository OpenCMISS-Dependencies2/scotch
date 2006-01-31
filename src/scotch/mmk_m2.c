/************************************************************/
/**                                                        **/
/**   NAME       : mmk_m2.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates source meshes for tridimen-     **/
/**                sional finite element grids.            **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 26 sep 2002     **/
/**                                 to   : 06 feb 2003     **/
/**                                                        **/
/**   NOTES      : # The nodes and elements of the         **/
/**                  (dX,dY) mesh are numbered so that     **/
/**                  t(0,0) = 0, t(1,0) = 1,               **/
/**                  t(dX - 1, 0) = dX - 1, t(0,1) =       **/
/**                  dX, and t(x,y) = y * dX + x.          **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MMK_M2

#include "common.h"
#include "scotch.h"
#include "mmk_m2.h"

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "-", NULL, "w" },
                              { "-", NULL, "w" } };

static const int            C_nghbTab[3] = { 4, 2, 1 };

static const char *         C_usageList[] = {
  "mmk_m2 <dimX> [<dimY> [<output mesh file>]] <options>",
  "  -g<file>  : Output mesh geometry to <file>",
  "  -h        : Display this help",
  "  -V        : Print program version and copyright",
  NULL };

/****************************************/
/*                                      */
/* The main routine, which computes the */
/* source mesh description.             */
/*                                      */
/****************************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Num          e[2] = { 1, 1 };            /* Mesh element dimensions */
  SCOTCH_Num          n[2];                       /* Mesh node dimensions    */
  SCOTCH_Num          c[2];                       /* Vertex coordinates      */
  SCOTCH_Num          velmnbr;                    /* First node number       */
  uint                flag;                       /* Process flags           */
  int                 i;

  errorProg ("mmk_m2");

  flag = C_FLAGDEFAULT;                           /* Set default flags */

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  for (i = 0; i < C_FILENBR; i ++)                /* Set default stream pointers */
    C_fileTab[i].pntr = (C_fileTab[i].mode[0] == 'r') ? stdin : stdout;
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_paraNum < 2) {                        /* If number of parameters not reached */
        if ((e[C_paraNum ++] = atoi (argv[i])) < 1) { /* Get the dimension               */
          errorPrint ("main: invalid dimension (\"%s\")", argv[i]);
          return     (1);
        }
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        C_fileTab[C_fileNum ++].name = argv[i];
      else {
        errorPrint ("main: too many file names given");
        return     (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'G' :                                /* Output mesh geometry */
        case 'g' :
          flag |= C_FLAGGEOOUT;
          if (argv[i][2] != '\0')
            C_filenamegeoout = &argv[i][2];
          break;
        case 'H' :                                /* Give program usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "mmk_m2, version %s - F. Pellegrini\n", SCOTCH_VERSION);
          fprintf (stderr, "Copyright 2004 INRIA, France\n");
          fprintf (stderr, "This software is libre/free software under LGPL -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option (\"%s\")", argv[i]);
          return     (1);
      }
    }
  }

  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      if ((C_fileTab[i].pntr = fopen (C_fileTab[i].name, C_fileTab[i].mode)) == NULL) { /* Open the file */
        errorPrint ("main: cannot open file (%d)", i);
        return     (1);
      }
    }
  }

  n[0] = e[0] + 1;
  n[1] = e[1] + 1;
  velmnbr = e[0] * e[1];

  fprintf (C_filepntrmshout, "1\n%ld\t%ld\t%ld\n0\t%ld\t000\n", /* Print mesh file header */
             (long) velmnbr,
             (long) (n[0] * n[1]),
             (long) (((velmnbr + n[0] * n[1]) - (e[0] + e[1] + 1)) * 4),
             (long) velmnbr);

  for (c[1] = 0; c[1] < e[1]; c[1] ++) {          /* Output element neighbor list */
    for (c[0] = 0; c[0] < e[0]; c[0] ++) {
      fprintf (C_filepntrmshout, "4\t%ld\t%ld\t%ld\t%ld\n", /* Output neighbors of element */
                (long) (c[1] * n[0] + c[0]),
                (long) (c[1] * n[0] + c[0] + 1),
                (long) ((c[1] + 1) * n[0] + c[0]),
                (long) ((c[1] + 1) * n[0] + c[0] + 1));
    }
  }
  for (c[1] = 0; c[1] < n[1]; c[1] ++) {          /* Output node neighbor list */
    for (c[0] = 0; c[0] < n[0]; c[0] ++) {
      fprintf (C_filepntrmshout, "%d",            /* Output number of neighboring elements */
                C_nghbTab[(((c[0] != 0) && (c[0] != e[0])) ? 0 : 1) +
                          (((c[1] != 0) && (c[1] != e[1])) ? 0 : 1)]);
      if (c[1] != 0) {                            /* Output neighbors of nodes */
        if (c[0] != 0)
          fprintf (C_filepntrmshout, "\t%ld",
                    (long) ((c[1] - 1) * e[0] + (c[0] - 1)));
        if (c[0] != e[0])
          fprintf (C_filepntrmshout, "\t%ld",
                    (long) ((c[1] - 1) * e[0] + c[0]));
      }
      if (c[1] != e[1]) {
        if (c[0] != 0)
          fprintf (C_filepntrmshout, "\t%ld",
                    (long) (c[1] * e[0] + (c[0] - 1)));
        if (c[0] != e[0])
          fprintf (C_filepntrmshout, "\t%ld",
                    (long) (c[1] * e[0] + c[0]));
      }
      fprintf (C_filepntrmshout, "\n");
    }
  }

  if (flag & C_FLAGGEOOUT) {                      /* If geometry is wanted       */
    fprintf (C_filepntrgeoout, "2\n%ld\n",        /* Output geometry file header */
              (long) (velmnbr + n[0] * n[1]));

    for (c[1] = 0; c[1] < e[1]; c[1] ++) {        /* Output element coordinates */
      for (c[0] = 0; c[0] < e[0]; c[0] ++)
        fprintf (C_filepntrgeoout, "%ld\t%ld.5\t%ld.5\n",
                  (long) (c[1] * e[0] + c[0]),
                  (long) c[0],
                  (long) (e[1] - 1 - c[1]));
    }
    for (c[1] = 0; c[1] < n[1]; c[1] ++) {        /* Output node coordinates */
      for (c[0] = 0; c[0] < n[0]; c[0] ++)
        fprintf (C_filepntrgeoout, "%ld\t%ld\t%ld\n",
                  (long) (velmnbr + c[1] * n[0] + c[0]),
                  (long) c[0],
                  (long) (e[1] - c[1]));
    }
  }

#ifdef SCOTCH_DEBUG_MAIN1
  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      fclose (C_fileTab[i].pntr);                 /* Close the stream */
    }
  }
#endif /* SCOTCH_DEBUG_MAIN1 */

  return (0);
}
