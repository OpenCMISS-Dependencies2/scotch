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
/**   NAME       : symbol.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic matrix.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 21 mar 2002     **/
/**                                 to     21 mar 2002     **/
/**                # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.3  : from : 10 apr 2003     **/
/**                                 to     10 jun 2003     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_H
#define SYMBOL_VERSION              1

#ifdef CXREF_DOC
#ifndef COMMON_H
#include "common.h"
#endif /* COMMON_H */
#ifndef GRAPH_H
#include "graph.h"
#endif /* GRAPH_H */
#ifndef DOF_H
#include "dof.h"
#endif /* DOF_H */
#endif /* CXREF_DOC */

/*
**  The type and structure definitions.
*/

/*+ The column block structure. +*/

typedef struct SymbolCblk_ {
  INT                       fcolnum;              /*+ First column index               +*/
  INT                       lcolnum;              /*+ Last column index (inclusive)    +*/
  INT                       bloknum;              /*+ First block in column (diagonal) +*/
} SymbolCblk;

/*+ The column block structure. +*/

typedef struct SymbolBlok_ {
  INT                       frownum;              /*+ First row index            +*/
  INT                       lrownum;              /*+ Last row index (inclusive) +*/
  INT                       cblknum;              /*+ Facing column block        +*/
  INT                       levfval;              /*+ Level-of-fill value        +*/
} SymbolBlok;

/*+ The symbolic block matrix. +*/

typedef struct SymbolMatrix_ {
  INT                       baseval;              /*+ Base value for numberings         +*/
  INT                       cblknbr;              /*+ Number of column blocks           +*/
  INT                       bloknbr;              /*+ Number of blocks                  +*/
  SymbolCblk * restrict     cblktab;              /*+ Array of column blocks [+1,based] +*/
  SymbolBlok * restrict     bloktab;              /*+ Array of blocks [based]           +*/
  INT                       nodenbr;              /*+ Number of nodes in matrix         +*/
} SymbolMatrix;

/*+ The type of cost computations. +*/

typedef enum SymbolCostType_ {
  SYMBOLCOSTLDLT                                  /*+ Crout (i.e. LDLt) cost function +*/
} SymbolCostType;

/*
**  The function prototypes.
*/

#ifndef SYMBOL
#define static
#endif

int                         symbolInit          (SymbolMatrix * const symbptr);
void                        symbolExit          (SymbolMatrix * const symbptr);
void                        symbolRealloc       (SymbolMatrix * const symbptr);
int                         symbolLoad          (SymbolMatrix * const symbptr, FILE * const stream);
int                         symbolSave          (const SymbolMatrix * const symbptr, FILE * const stream);
int                         symbolCheck         (const SymbolMatrix * const symbptr);
int                         symbolDraw          (const SymbolMatrix * const symbptr, FILE * const stream);
int                         symbolDrawFunc      (const SymbolMatrix * const symbptr, int (*) (const SymbolMatrix * const, const SymbolBlok * const, void * const, float * const), int (*) (const SymbolMatrix * const, const SymbolBlok * const, void * const, float * const), void * const, FILE * const stream);
void                        symbolDrawColor     (const INT labl, float * const coloptr);
#ifdef DOF_H
int                         symbolCost          (const SymbolMatrix * const symbptr, const Dof * const deofptr, const SymbolCostType typeval, double * const nnzptr, double * const opcptr);
int                         symbolCosti         (const SymbolMatrix * const symbptr, const Dof * const deofptr, const SymbolCostType typeval, const INT levfval, double * const nnzptr, double * const opcptr);
int                         symbolLevf          (const SymbolMatrix * const symbptr, INT * const levfmax, INT ** const levftab);
int                         symbolTree          (const SymbolMatrix * const symbptr, const Dof * const deofptr, INT * const leafnbr, INT * const heigmin, INT * const heigmax, double * const heigavg, double * const heigdlt);
int                         symbolNonzeros      (const SymbolMatrix * const symbptr, FILE * const stream);
#endif /* DOF_H */

#undef static
