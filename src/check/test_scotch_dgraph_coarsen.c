/* Copyright 2014,2015,2018,2021 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_dgraph_coarsen.c            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_dgraphCoarsen() routine.     **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 28 sep 2014     **/
/**                                 to   : 22 may 2018     **/
/**                # Version 6.1  : from : 16 jun 2021     **/
/**                                 to   : 28 dec 2021     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <mpi.h>
#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <pthread.h>
#include <unistd.h>

#include "ptscotch.h"

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (

int                 argc,
char *              argv[])
{
  MPI_Comm            proccomm;
  int                 procglbnbr;                 /* Number of processes sharing graph data */
  int                 proclocnum;                 /* Number of this process                 */
  SCOTCH_Num          vertglbnbr;
  SCOTCH_Num          vertlocnbr;
  SCOTCH_Dgraph       grafdat;
  SCOTCH_Dgraph       coargrafdat;
  SCOTCH_Num          coarvertglbnbr;
  SCOTCH_Num          coarvertlocnbr;
  double              coarrat;
  FILE *              file;
#ifdef SCOTCH_PTHREAD
  int                 thrdlvlreqval;
  int                 thrdlvlproval;
#endif /* SCOTCH_PTHREAD */
  int                 i;

  SCOTCH_errorProg (argv[0]);

#ifdef SCOTCH_PTHREAD
  thrdlvlreqval = MPI_THREAD_MULTIPLE;
  if (MPI_Init_thread (&argc, &argv, thrdlvlreqval, &thrdlvlproval) != MPI_SUCCESS)
    SCOTCH_errorPrint ("main: Cannot initialize (1)");
  if (thrdlvlreqval > thrdlvlproval)
    SCOTCH_errorPrint ("main: MPI implementation is not thread-safe: recompile without SCOTCH_PTHREAD");
#else /* SCOTCH_PTHREAD */
  if (MPI_Init (&argc, &argv) != MPI_SUCCESS)
    SCOTCH_errorPrint ("main: Cannot initialize (2)");
#endif /* SCOTCH_PTHREAD */

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  proccomm = MPI_COMM_WORLD;
  MPI_Comm_size (proccomm, &procglbnbr);          /* Get communicator data */
  MPI_Comm_rank (proccomm, &proclocnum);

  fprintf (stderr, "Proc %2d of %2d, pid %d\n", proclocnum, procglbnbr, getpid ());

#ifdef SCOTCH_CHECK_NOAUTO
  if (proclocnum == 0) {                          /* Synchronize on keybord input */
    char           c;

    printf ("Waiting for key press...\n");
    scanf ("%c", &c);
  }
#endif /* SCOTCH_CHECK_NOAUTO */

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize for debug */
    SCOTCH_errorPrint ("main: cannot communicate (1)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_dgraphInit (&grafdat, proccomm) != 0) { /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph (1)");
    exit (EXIT_FAILURE);
  }

  file = NULL;
  if ((proclocnum == 0) &&
      ((file = fopen (argv[1], "r")) == NULL)) {
    SCOTCH_errorPrint ("main: cannot open graph file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_dgraphLoad (&grafdat, file, -1, 0) != 0) {
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  if (file != NULL)
    fclose (file);

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize for debug */
    SCOTCH_errorPrint ("main: cannot communicate (2)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_dgraphData (&grafdat, NULL, &vertglbnbr, &vertlocnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  coarrat = 0.8;                                  /* Lazy coarsening ratio */

  for (i = 0; i < 3; i ++) {                      /* For all test cases */
    SCOTCH_Num *        multloctab;
    SCOTCH_Num          multlocsiz;
    SCOTCH_Num          foldval;
    char *              foldstr;
    char *              coarstr;
    int                 procnum;
    int                 o;

    switch (i) {
      case 0 :
        foldval = SCOTCH_COARSENNONE;
        foldstr = "Plain coarsening";
        multlocsiz = vertlocnbr;
        break;
      case 1 :
        foldval = SCOTCH_COARSENFOLD;
        foldstr = "Folding";
        multlocsiz = (procglbnbr == 1) ? vertlocnbr : (((SCOTCH_Num) ((double) (vertglbnbr * 2) * coarrat)) / procglbnbr) + 1; /* Max ratio FOLD is 2 -> 1 */
        break;
      case 2 :
        foldval = SCOTCH_COARSENFOLDDUP;
        foldstr = "Folding with duplication";
	multlocsiz = (procglbnbr == 1) ? vertlocnbr : (((SCOTCH_Num) ((double) (vertglbnbr * 2) * coarrat)) / (procglbnbr - (procglbnbr % 2))) + 1; /* Max ratio FOLD-DUP is 3 -> 1 */
        break;
    }

    if (proclocnum == 0)
      printf ("%s\n", foldstr);

    if ((multloctab = malloc (vertglbnbr * 2 * sizeof (SCOTCH_Num))) == NULL) { /* Allocate maximum size for security */
      SCOTCH_errorPrint ("main: cannot allocate multinode array");
      exit (EXIT_FAILURE);
    }

    if (SCOTCH_dgraphInit (&coargrafdat, proccomm) != 0) { /* Initialize band graph */
      SCOTCH_errorPrint ("main: cannot initialize graph (2)");
      exit (EXIT_FAILURE);
    }

    o = SCOTCH_dgraphCoarsen (&grafdat, 0, coarrat, foldval, &coargrafdat, multloctab);

    SCOTCH_dgraphData (&coargrafdat, NULL, &coarvertglbnbr, &coarvertlocnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    for (procnum = 0; procnum < procglbnbr; procnum ++) {
      switch (o) {
        case 0 :
          coarstr = "coarse graph created";
          break;
        case 1 :
          coarstr = "graph could not be coarsened";
          break;
        case 2 :
          coarstr = "folded graph not created here";
          break;
        case 3 :
          coarstr = "cannot create coarse graph";
          break;
      }

      if (procnum == proclocnum)
        printf ("%d: %s (%ld / %ld / %ld / %f)\n", procnum, coarstr, (long) multlocsiz, (long) coarvertlocnbr, (long) vertlocnbr,
                (double) coarvertglbnbr / (double) vertglbnbr);

      MPI_Barrier (proccomm);
    }

    if (coarvertlocnbr > multlocsiz) {
      SCOTCH_errorPrint ("main: invalid local multinode array size");
      exit (EXIT_FAILURE);
    }

    SCOTCH_dgraphExit (&coargrafdat);
    free (multloctab);
  }

  SCOTCH_dgraphExit (&grafdat);

  MPI_Finalize ();
  exit (EXIT_SUCCESS);
}
