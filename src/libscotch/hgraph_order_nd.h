/* Copyright 2004,2007,2010,2018,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_nd.h                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the halo graph nested dissection    **/
/**                ordering algorithm.                     **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 17 oct 1996     **/
/**                                 to   : 18 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 13 mar 1999     **/
/**                # Version 4.0  : from : 03 jan 2002     **/
/**                                 to   : 24 apr 2004     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to   : 04 nov 2010     **/
/**                # Version 6.0  : from : 07 jun 2018     **/
/**                                 to   : 07 jun 2018     **/
/**                # Version 7.0  : from : 26 apr 2021     **/
/**                                 to   : 11 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct HgraphOrderNdParam_ {
  Strat *                   sepstrat;             /*+ Separation strategy         +*/
  Strat *                   ordstratlea;          /*+ Leaf ordering strategy      +*/
  Strat *                   ordstratsep;          /*+ Separator ordering strategy +*/
} HgraphOrderNdParam;

/*+ This structure holds the splitting parameters. +*/

typedef struct HgraphOrderNdSplit2_ {
  Gnum                      vnumnbr;              /*+ Number of induced vertices               +*/
  Gnum *                    vnumtab;              /*+ Pointer to list of induced vertices      +*/
  Gnum                      vhalmax;              /*+ Upper bound on induced halo vertices     +*/
  Gnum                      ordenum;              /*+ Local start index of inverse permutation +*/
  OrderCblk *               cblkptr;              /*+ Column block to process                  +*/
} HgraphOrderNdSplit2;

typedef struct HgraphOrderNdSplit_ {
  HgraphOrderNdSplit2         splttab[2];         /*+ Array of induced subgraph data +*/
  const Hgraph *              grafptr;            /*+ Original graph                 +*/
  Order *                     ordeptr;            /*+ Pointer to ordering            +*/
  const HgraphOrderNdParam *  paraptr;            /*+ Nested dissection parameters   +*/
  int *                       revaptr;            /*+ Pointer to return value        +*/
} HgraphOrderNdSplit;

/*
**  The function prototypes.
*/

#ifdef SCOTCH_HGRAPH_ORDER_ND
static void                 hgraphOrderNd2      (Context * restrict const, const int, const HgraphOrderNdSplit * const);
#endif /* SCOTCH_HGRAPH_ORDER_ND */

int                         hgraphOrderNd       (Hgraph * const, Order * const, const Gnum, OrderCblk * const, const HgraphOrderNdParam * const);
