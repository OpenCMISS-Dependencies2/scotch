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
/**   NAME       : arch_build.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the architecture building routine.  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 29 may 1997     **/
/**                                 to     01 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to     29 nov 2003     **/
/**                                                        **/
/**   NOTES      : # This file contains pieces of code     **/
/**                  extracted from release 3.1 of         **/
/**                  "amk_src.c".                          **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Job to process. +*/

typedef struct ArchBuildJob_ {
  struct ArchBuildJob_ *    joblink;              /*+ Link to job pool                        +*/
  ArchDomNum                domnum;               /*+ Mapping domain to which vertices belong +*/
  Graph                     grafdat;              /*+ Job graph data                          +*/
} ArchBuildJob;

/*+ Vertex distance information. +*/

typedef struct ArchBuildDistElem_ {
  int                       queued;               /*+ Flag set if vertex queued  +*/
  Anum                      distval;              /*+ Distance to initial vertex +*/
} ArchBuildDistElem;

/*+ Queue element. +*/

typedef struct ArchBuildQueuElem_ {
  Gnum                      vertnum;              /*+ Vertex number in source graph +*/
  Anum                      distval;              /*+ Distance reached              +*/
} ArchBuildQueuElem;

/*
**  The function prototypes.
*/

#ifndef ARCH_BUILD
#define static
#endif

static void                 archBuildJobExit    (ArchBuildJob *);

int                         archBuild           (Arch * const, const Graph * const, const VertList * const, const Strat * const);

#undef static
