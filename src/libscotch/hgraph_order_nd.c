/* Copyright 2004,2007,2010,2012,2014,2016,2018,2019,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_nd.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders graphs using the     **/
/**                nested dissection algorithm.            **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 17 oct 1996     **/
/**                                 to   : 21 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 13 mar 1999     **/
/**                # Version 4.0  : from : 03 jan 2002     **/
/**                                 to   : 24 dec 2004     **/
/**                # Version 5.0  : from : 19 dec 2006     **/
/**                                 to   : 25 jul 2007     **/
/**                # Version 5.1  : from : 24 oct 2010     **/
/**                                 to   : 24 oct 2010     **/
/**                # Version 6.0  : from : 17 oct 2012     **/
/**                                 to   : 23 may 2018     **/
/**                # Version 6.1  : from : 01 nov 2021     **/
/**                                 to   : 21 nov 2021     **/
/**                # Version 7.0  : from : 05 may 2019     **/
/**                                 to   : 11 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SCOTCH_HGRAPH_ORDER_ND

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_nd.h"
#include "hgraph_order_st.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"

/***************************/
/*                         */
/* Multi-threaded routine. */
/*                         */
/***************************/

static
void
hgraphOrderNd2 (
Context * restrict const          contptr,        /*+ (Sub-)context                          +*/
const int                         spltnum,        /*+ Rank of sub-context in initial context +*/
const HgraphOrderNdSplit * const  spltptr)
{
  Hgraph              orggrafdat;
  Hgraph              indgrafdat;
  int                 o;

  orggrafdat = *spltptr->grafptr;                 /* Assign new context to work graph */
  orggrafdat.contptr = contptr;

  o = 0;
  if ((hgraphInduceList (&orggrafdat, spltptr->splttab[spltnum].vnumnbr, spltptr->splttab[spltnum].vnumtab,
                         spltptr->splttab[spltnum].vhalmax, &indgrafdat)) != 0) {
    errorPrint ("hgraphOrderNd2: cannot build induced subgraph");
    o = 1;
  }
  if (o == 0) {
    o = hgraphOrderNd (&indgrafdat, spltptr->ordeptr, spltptr->splttab[spltnum].ordenum, spltptr->splttab[spltnum].cblkptr, spltptr->paraptr);

    hgraphExit (&indgrafdat);
  }
  if (o != 0) {                                   /* Report any error */
#ifdef SCOTCH_PTHREAD
    pthread_mutex_lock (&spltptr->ordeptr->mutedat); /* Use ordering lock to avoid race condition */
#endif /* SCOTCH_PTHREAD */
    *spltptr->revaptr = 1;
#ifdef SCOTCH_PTHREAD
    pthread_mutex_unlock (&spltptr->ordeptr->mutedat);
#endif /* SCOTCH_PTHREAD */
  }
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the ordering.
** It returns:
** - 0   : if the ordering could be computed.
** - !0  : on error.
*/

int
hgraphOrderNd (
Hgraph * restrict const                   grafptr,
Order * restrict const                    ordeptr,
const Gnum                                ordenum,
OrderCblk * restrict const                cblkptr,
const HgraphOrderNdParam * restrict const paraptr)
{
  Gnum *                    vspvnumptr[3];        /* Pointers to vertex lists to fill   */
  VertList                  vsplisttab[3];        /* Array of separated part lists      */
  Vgraph                    vspgrafdat;           /* Vertex separation graph data       */
  Gnum                      vspvertnum;           /* Current vertex in separation graph */
  int                       o;

  hgraphUnhalo (grafptr, &vspgrafdat.s);          /* Keep only non-halo vertices for separation */
  if ((vspgrafdat.frontab = (Gnum *) memAlloc (vspgrafdat.s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("hgraphOrderNd: out of memory (1)");
    return (1);
  }
  if ((vspgrafdat.parttax = (GraphPart *) memAlloc (vspgrafdat.s.vertnbr * sizeof (GraphPart))) == NULL) {
    errorPrint ("hgraphOrderNd: out of memory (2)");
    memFree    (vspgrafdat.frontab);
    return (1);
  }
  memSet (vspgrafdat.parttax, 0, vspgrafdat.s.vertnbr * sizeof (GraphPart)); /* Set all vertices to part 0 */
  vspgrafdat.parttax    -= vspgrafdat.s.baseval;
  vspgrafdat.s.flagval  |= VGRAPHFREEFRON | VGRAPHFREEPART;
  vspgrafdat.dwgttab[0]  =                        /* Parts should have about the same weights */
  vspgrafdat.dwgttab[1]  = 1;
  vspgrafdat.compload[0] = vspgrafdat.s.velosum;
  vspgrafdat.compload[1] = 0;
  vspgrafdat.compload[2] = 0;
  vspgrafdat.comploaddlt = vspgrafdat.s.velosum;
  vspgrafdat.compsize[0] = vspgrafdat.s.vertnbr;
  vspgrafdat.compsize[1] = 0;
  vspgrafdat.fronnbr     = 0;
  vspgrafdat.levlnum     = grafptr->levlnum;      /* Set level of separation graph as that of halo graph */
  vspgrafdat.contptr     = grafptr->contptr;      /* Use same execution context                          */

  if (vgraphSeparateSt (&vspgrafdat, paraptr->sepstrat) != 0) { /* Separate vertex-separation graph */
    vgraphExit (&vspgrafdat);
    return (1);
  }

  if ((vspgrafdat.compsize[0] == 0) ||            /* If could not separate more */
      (vspgrafdat.compsize[1] == 0)) {
    vgraphExit    (&vspgrafdat);                  /* Free useless space */
    hgraphOrderSt (grafptr, ordeptr, ordenum, cblkptr, paraptr->ordstratlea); /* Order this leaf */
    return (0);                                   /* Leaf has been processed */
  }

  vsplisttab[0].vnumnbr = vspgrafdat.compsize[0]; /* Build vertex lists within frontier array */
  vsplisttab[0].vnumtab = vspgrafdat.frontab + vspgrafdat.fronnbr;
  vsplisttab[1].vnumnbr = vspgrafdat.compsize[1];
  vsplisttab[1].vnumtab = vsplisttab[0].vnumtab + vsplisttab[0].vnumnbr;
  vsplisttab[2].vnumnbr = vspgrafdat.fronnbr;
  vsplisttab[2].vnumtab = vspgrafdat.frontab;
  vspvnumptr[0] = vsplisttab[0].vnumtab;
  vspvnumptr[1] = vsplisttab[1].vnumtab;
  vspvnumptr[2] = vsplisttab[2].vnumtab;
  for (vspvertnum = vspgrafdat.s.baseval; vspvertnum < vspgrafdat.s.vertnnd; vspvertnum ++) { /* Fill lists */
    *vspvnumptr[vspgrafdat.parttax[vspvertnum]] ++ = vspvertnum;
#ifdef SCOTCH_DEBUG_HGRAPH2
    if (vspgrafdat.parttax[vspvertnum] != 2) {    /* If vertex does not separate */
      Gnum                vspedgenum;
      GraphPart           vsppartnum;

      vsppartnum = 1 - vspgrafdat.parttax[vspvertnum]; /* Get opposite part value */
      for (vspedgenum = vspgrafdat.s.verttax[vspvertnum];
           vspedgenum < vspgrafdat.s.vendtax[vspvertnum]; vspedgenum ++) {
        if (vspgrafdat.parttax[vspgrafdat.s.edgetax[vspedgenum]] == vsppartnum) { /* If an edge crosses the separator */
          errorPrint ("hgraphOrderNd: internal error (1)");
          vgraphExit (&vspgrafdat);
          return (1);
        }
      }
    }
#endif /* SCOTCH_DEBUG_HGRAPH2 */
  }
#ifdef SCOTCH_DEBUG_HGRAPH2
  if ((vspvnumptr[0] != vsplisttab[0].vnumtab + vsplisttab[0].vnumnbr) ||
      (vspvnumptr[1] != vsplisttab[1].vnumtab + vsplisttab[1].vnumnbr) ||
      (vspvnumptr[2] != vsplisttab[2].vnumtab + vsplisttab[2].vnumnbr)) {
    errorPrint ("hgraphOrderNd: internal error (2)");
    vgraphExit (&vspgrafdat);
    return (1);
  }
#endif /* SCOTCH_DEBUG_HGRAPH2 */

  memFree (vspgrafdat.parttax + vspgrafdat.s.baseval); /* Free useless space    */
  vspgrafdat.parttax = NULL;                      /* Prevent subsequent freeing */

  cblkptr->typeval = ORDERCBLKNEDI;               /* Node becomes a nested dissection node */
  if ((cblkptr->cblktab = (OrderCblk *) memAlloc (3 * sizeof (OrderCblk))) == NULL) {
    errorPrint ("hgraphOrderNd: out of memory (3)");
    vgraphExit (&vspgrafdat);
    return (1);
  }
  cblkptr->cblktab[0].typeval = ORDERCBLKLEAF;    /* Build column blocks */
  cblkptr->cblktab[0].vnodnbr = vsplisttab[0].vnumnbr;
  cblkptr->cblktab[0].cblknbr = 0;
  cblkptr->cblktab[0].cblktab = NULL;
  cblkptr->cblktab[1].typeval = ORDERCBLKLEAF;
  cblkptr->cblktab[1].vnodnbr = vsplisttab[1].vnumnbr;
  cblkptr->cblktab[1].cblknbr = 0;
  cblkptr->cblktab[1].cblktab = NULL;

  if (vsplisttab[2].vnumnbr != 0) {               /* If separator not empty */
    Hgraph              indgrafdat;

    cblkptr->cblknbr = 3;                         /* It is a three-cell tree node */
#ifdef SCOTCH_PTHREAD
    pthread_mutex_lock (&ordeptr->mutedat);
#endif /* SCOTCH_PTHREAD */
    ordeptr->cblknbr += 2;                        /* Two more column blocks created */
    ordeptr->treenbr += 3;                        /* Three more tree nodes created  */
#ifdef SCOTCH_PTHREAD
    pthread_mutex_unlock (&ordeptr->mutedat);
#endif /* SCOTCH_PTHREAD */

    cblkptr->cblktab[2].typeval = ORDERCBLKLEAF;
    cblkptr->cblktab[2].vnodnbr = vsplisttab[2].vnumnbr;
    cblkptr->cblktab[2].cblknbr = 0;
    cblkptr->cblktab[2].cblktab = NULL;

    if (graphInduceList (&grafptr->s, vsplisttab[2].vnumnbr, vsplisttab[2].vnumtab, &indgrafdat.s) != 0) { /* Perform non-halo induction for separator, as it will get highest numbers */
      errorPrint ("hgraphOrderNd: cannot build induced subgraph (1)");
      vgraphExit (&vspgrafdat);
      return (1);
    }
    indgrafdat.vnohnbr = indgrafdat.s.vertnbr;    /* Fill halo graph structure of non-halo graph */
    indgrafdat.vnohnnd = indgrafdat.s.vertnnd;
    indgrafdat.vnhdtax = indgrafdat.s.vendtax;
    indgrafdat.vnlosum = indgrafdat.s.velosum;
    indgrafdat.enohnbr = indgrafdat.s.edgenbr;
    indgrafdat.enlosum = indgrafdat.s.edlosum;
    indgrafdat.levlnum = grafptr->levlnum;        /* Separator graph is at level of original graph */
    indgrafdat.contptr = grafptr->contptr;        /* Use same execution context                    */

    o = hgraphOrderSt (&indgrafdat, ordeptr, ordenum + vsplisttab[0].vnumnbr + vsplisttab[1].vnumnbr,
                       cblkptr->cblktab + 2, paraptr->ordstratsep);
    hgraphExit (&indgrafdat);
  }
  else {                                          /* Separator is empty         */
    cblkptr->cblknbr = 2;                         /* It is a two-cell tree node */
#ifdef SCOTCH_PTHREAD
    pthread_mutex_lock (&ordeptr->mutedat);
#endif /* SCOTCH_PTHREAD */
    ordeptr->cblknbr ++;                          /* One more column block created */
    ordeptr->treenbr += 2;                        /* Two more tree nodes created   */
#ifdef SCOTCH_PTHREAD
    pthread_mutex_unlock (&ordeptr->mutedat);
#endif /* SCOTCH_PTHREAD */
    o = 0;                                        /* No separator ordering computed */
  }
  if (o == 0) {
    HgraphOrderNdSplit  spltdat;                  /* Parameters for context splitting */

    spltdat.splttab[0].vnumnbr = vsplisttab[0].vnumnbr;
    spltdat.splttab[0].vnumtab = vsplisttab[0].vnumtab;
    spltdat.splttab[0].vhalmax = vsplisttab[2].vnumnbr + grafptr->s.vertnbr - grafptr->vnohnbr;
    spltdat.splttab[0].ordenum = ordenum;
    spltdat.splttab[0].cblkptr = cblkptr->cblktab;
    spltdat.splttab[1].vnumnbr = vsplisttab[1].vnumnbr;
    spltdat.splttab[1].vnumtab = vsplisttab[1].vnumtab;
    spltdat.splttab[1].vhalmax = vsplisttab[2].vnumnbr + grafptr->s.vertnbr - grafptr->vnohnbr;
    spltdat.splttab[1].ordenum = ordenum + vsplisttab[0].vnumnbr;
    spltdat.splttab[1].cblkptr = cblkptr->cblktab + 1;
    spltdat.grafptr = grafptr;
    spltdat.ordeptr = ordeptr;
    spltdat.paraptr = paraptr;
    spltdat.revaptr = &o;

#ifndef HGRAPHORDERNDNOTHREAD
    if (contextThreadLaunchSplit (grafptr->contptr, (ContextSplitFunc) hgraphOrderNd2, &spltdat) != 0) /* If counld not split context to run concurrently */
#endif /* HGRAPHORDERNDNOTHREAD */
    {
      hgraphOrderNd2 (grafptr->contptr, 0, &spltdat); /* Run tasks in sequence */
      if (o == 0)
        hgraphOrderNd2 (grafptr->contptr, 1, &spltdat);
    }
  }

  vgraphExit (&vspgrafdat);

  return (o);
}
