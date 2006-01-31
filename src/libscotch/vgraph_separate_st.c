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
/**   NAME       : vgraph_separate_st.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the global         **/
/**                separation strategy and method tables.  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 25 oct 1996     **/
/**                                 to     14 nov 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     31 may 1999     **/
/**                # Version 4.0  : from : 06 jan 2002     **/
/**                                 to     31 aug 2005     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_ST

#include "module.h"
#include "common.h"
#include "gain.h"
#include "parser.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "vgraph.h"
#include "vgraph_separate_es.h"
#include "vgraph_separate_fm.h"
#include "vgraph_separate_gg.h"
#include "vgraph_separate_gp.h"
#include "vgraph_separate_ml.h"
#include "vgraph_separate_th.h"
#include "vgraph_separate_st.h"
#include "vgraph_separate_vw.h"
#include "vgraph_separate_zr.h"

/*
**  The static and global variables.
*/

static Vgraph               vgraphdummy;          /* Dummy separator graph for offset computations */

static union {
  VgraphSeparateEsParam     param;
  StratNodeMethodData       padding;
} vgraphseparatedefaultes = { { &stratdummy, VGRAPHSEPAESWIDTHTHIN } };

static union {
  VgraphSeparateFmParam     param;
  StratNodeMethodData       padding;
} vgraphseparatedefaultfm = { { 200, 1000, 0.1L } };

static union {
  VgraphSeparateGgParam     param;
  StratNodeMethodData       padding;
} vgraphseparatedefaultgg = { { 5 } };

static union {
  VgraphSeparateGpParam     param;
  StratNodeMethodData       padding;
} vgraphseparatedefaultgp = { { 9 } };

static union {
  VgraphSeparateMlParam     param;
  StratNodeMethodData       padding;
} vgraphseparatedefaultml = { { 100, 0.8L, GRAPHCOARHEM, &stratdummy, &stratdummy } };

static StratMethodTab       vgraphseparatestmethtab[] = { /* Graph separation methods array */
                              { VGRAPHSEPASTMETHES, "e",  vgraphSeparateEs, &vgraphseparatedefaultes },
                              { VGRAPHSEPASTMETHFM, "f",  vgraphSeparateFm, &vgraphseparatedefaultfm },
                              { VGRAPHSEPASTMETHGG, "h",  vgraphSeparateGg, &vgraphseparatedefaultgg },
                              { VGRAPHSEPASTMETHGP, "g",  vgraphSeparateGp, &vgraphseparatedefaultgp },
                              { VGRAPHSEPASTMETHML, "m",  vgraphSeparateMl, &vgraphseparatedefaultml },
                              { VGRAPHSEPASTMETHVW, "w",  vgraphSeparateVw, NULL },
                              { VGRAPHSEPASTMETHZR, "z",  vgraphSeparateZr, NULL },
                              { -1,                 NULL, NULL,             NULL } };

static StratParamTab        vgraphseparatestparatab[] = { /* Graph separation method parameter list */
                              { VGRAPHSEPASTMETHES,  STRATPARAMSTRAT,  "strat",
                                (byte *) &vgraphseparatedefaultes.param,
                                (byte *) &vgraphseparatedefaultes.param.strat,
                                (void *) &bgraphbipartststratab },
                              { VGRAPHSEPASTMETHES,  STRATPARAMCASE,   "type",
                                (byte *) &vgraphseparatedefaultes.param,
                                (byte *) &vgraphseparatedefaultes.param.widtval,
                                (void *) "tf" },
                              { VGRAPHSEPASTMETHFM,  STRATPARAMINT,    "move",
                                (byte *) &vgraphseparatedefaultfm.param,
                                (byte *) &vgraphseparatedefaultfm.param.movenbr,
                                NULL },
                              { VGRAPHSEPASTMETHFM,  STRATPARAMINT,    "pass",
                                (byte *) &vgraphseparatedefaultfm.param,
                                (byte *) &vgraphseparatedefaultfm.param.passnbr,
                                NULL },
                              { VGRAPHSEPASTMETHFM,  STRATPARAMDOUBLE, "bal",
                                (byte *) &vgraphseparatedefaultfm.param,
                                (byte *) &vgraphseparatedefaultfm.param.deltrat,
                                NULL },
                              { VGRAPHSEPASTMETHGG,  STRATPARAMINT,    "pass",
                                (byte *) &vgraphseparatedefaultgg.param,
                                (byte *) &vgraphseparatedefaultgg.param.passnbr,
                                NULL },
                              { VGRAPHSEPASTMETHGP,  STRATPARAMINT,    "pass",
                                (byte *) &vgraphseparatedefaultgp.param,
                                (byte *) &vgraphseparatedefaultgp.param.passnbr,
                                NULL },
                              { VGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "asc",
                                (byte *) &vgraphseparatedefaultml.param,
                                (byte *) &vgraphseparatedefaultml.param.stratasc,
                                (void *) &vgraphseparateststratab },
                              { VGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "low",
                                (byte *) &vgraphseparatedefaultml.param,
                                (byte *) &vgraphseparatedefaultml.param.stratlow,
                                (void *) &vgraphseparateststratab },
                              { VGRAPHSEPASTMETHML,  STRATPARAMCASE,   "type",
                                (byte *) &vgraphseparatedefaultml.param,
                                (byte *) &vgraphseparatedefaultml.param.coartype,
                                (void *) "hscd" },
                              { VGRAPHSEPASTMETHML,  STRATPARAMINT,    "vert",
                                (byte *) &vgraphseparatedefaultml.param,
                                (byte *) &vgraphseparatedefaultml.param.coarnbr,
                                NULL },
                              { VGRAPHSEPASTMETHML,  STRATPARAMDOUBLE, "rat",
                                (byte *) &vgraphseparatedefaultml.param,
                                (byte *) &vgraphseparatedefaultml.param.coarrat,
                                NULL },
                              { VGRAPHSEPASTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        vgraphseparatestcondtab[] = { /* Graph condition parameter table */
                              { STRATNODECOND,       STRATPARAMINT,    "levl",
                                (byte *) &vgraphdummy,
                                (byte *) &vgraphdummy.levlnum,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "load",
                                (byte *) &vgraphdummy,
                                (byte *) &vgraphdummy.s.velosum,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "vert",
                                (byte *) &vgraphdummy,
                                (byte *) &vgraphdummy.s.vertnbr,
                                NULL },
                              { STRATNODENBR,        STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    vgraphseparateststratab = { /* Strategy tables for vertex separation methods */
                              vgraphseparatestmethtab,
                              vgraphseparatestparatab,
                              vgraphseparatestcondtab };

/*******************************************/
/*                                         */
/* This is the generic separation routine. */
/*                                         */
/*******************************************/

/* This routine computes the separation of
** the given graph according to the given
** strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
vgraphSeparateSt (
Vgraph * restrict const       grafptr,            /*+ Separation graph    +*/
const Strat * restrict const  strat)              /*+ Separation strategy +*/
{
  StratTest           val;
  VgraphStore         savetab[2];                 /* Results of the two strategies */
  int                 o;

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("vgraphSeparateSt: invalid type specification for parser variables");
    return     (1);
  }
  if ((sizeof (VgraphSeparateFmParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (VgraphSeparateGgParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (VgraphSeparateGpParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (VgraphSeparateMlParam) > sizeof (StratNodeMethodData))) {
    errorPrint ("vgraphSeparateSt: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
#ifdef SCOTCH_DEBUG_VGRAPH1
  if ((strat->tabl != &vgraphseparateststratab) &&
      (strat       != &stratdummy)) {
    errorPrint ("vgraphSeparateSt: invalid parameter (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH1 */

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      o = vgraphSeparateSt (grafptr, strat->data.concat.strat[0]); /* Apply first strategy          */
      if (o == 0)                                 /* If it worked all right                         */
        o |= vgraphSeparateSt (grafptr, strat->data.concat.strat[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_VGRAPH2
        if ((val.typetest != STRATTESTVAL) ||
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("vgraphSeparateSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                      */
          o = vgraphSeparateSt (grafptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement          */
            o = vgraphSeparateSt (grafptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :
      if (((vgraphStoreInit (grafptr, &savetab[0])) != 0) || /* Allocate save areas */
          ((vgraphStoreInit (grafptr, &savetab[1])) != 0)) {
        errorPrint      ("vgraphSeparateSt: out of memory");
        vgraphStoreExit (&savetab[0]);
        return          (1);
      }

      vgraphStoreSave  (grafptr, &savetab[1]);    /* Save initial bipartition           */
      vgraphSeparateSt (grafptr, strat->data.select.strat[0]); /* Apply first strategy  */
      vgraphStoreSave  (grafptr, &savetab[0]);    /* Save its result                    */
      vgraphStoreUpdt  (grafptr, &savetab[1]);    /* Restore initial bipartition        */
      vgraphSeparateSt (grafptr, strat->data.select.strat[1]); /* Apply second strategy */

      if ( (savetab[0].fronnbr <  grafptr->fronnbr) || /* If first strategy is better */
          ((savetab[0].fronnbr == grafptr->fronnbr) &&
           (abs (savetab[0].comploaddlt) < abs (grafptr->comploaddlt))))
        vgraphStoreUpdt (grafptr, &savetab[0]);   /* Restore its result */

      vgraphStoreExit (&savetab[0]);              /* Free both save areas */
      vgraphStoreExit (&savetab[1]);
      break;
#ifdef SCOTCH_DEBUG_VGRAPH1
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_VGRAPH1 */
    default :
#endif /* SCOTCH_DEBUG_VGRAPH1 */
      return (strat->tabl->methtab[strat->data.method.meth].func (grafptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_VGRAPH1
    default :
      errorPrint ("vgraphSeparateSt: invalid parameter (2)");
      return     (1);
#endif /* SCOTCH_DEBUG_VGRAPH1 */
  }
  return (o);
}
