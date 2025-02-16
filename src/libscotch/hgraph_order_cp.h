/* Copyright 2004,2007,2010,2014,2018,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_cp.h                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the graph compression         **/
/**                ordering routine.                       **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 29 aug 1998     **/
/**                                 to   : 09 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to   : 03 jan 1999     **/
/**                # Version 4.0  : from : 01 jan 2003     **/
/**                                 to   : 01 jan 2003     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to   : 04 nov 2010     **/
/**                # Version 6.0  : from : 09 nov 2014     **/
/**                                 to   : 07 jun 2018     **/
/**                # Version 7.0  : from : 19 jan 2023     **/
/**                                 to   : 11 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Prime number for hashing vertex numbers. +*/

#define HGRAPHORDERCPHASHPRIME      17            /* Prime number */

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct HgraphOrderCpParam_ {
  double                    comprat;              /*+ Compression ratio threshold             +*/
  Strat *                   stratcpr;             /*+ Compressed subgraph ordering strategy   +*/
  Strat *                   stratunc;             /*+ Uncompressed subgraph ordering strategy +*/
} HgraphOrderCpParam;

/*+ This structure holds fine neighbor hashing data. +*/

typedef struct HgraphOrderCpHash_ {
  Gnum                      vertnum;              /*+ Origin vertex (i.e. pass) number +*/
  Gnum                      vertend;              /*+ Adjacent end vertex number       +*/
} HgraphOrderCpHash;

/*+ This structure holds coarse neighbor mate data. +*/

typedef struct HgraphOrderCpMate_ {
  Gnum                      coarvertend;          /*+ Adjacent coarse end vertex number +*/
  Gnum                      finevertend;          /*+ Adjacent end vertex number        +*/
} HgraphOrderCpMate;

/*
**  The function prototypes.
*/

#ifdef SCOTCH_HGRAPH_ORDER_CP
static Gnum                 hgraphOrderCpTree   (const Gnum * const, const Gnum * const, OrderCblk * const, const Gnum);
#endif /* SCOTCH_HGRAPH_ORDER_CP */

int                         hgraphOrderCp       (Hgraph * const, Order * const, const Gnum, OrderCblk * const, const HgraphOrderCpParam * const);
