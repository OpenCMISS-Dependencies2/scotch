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
/**   NAME       : vmesh_separate_gr.h                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the the graph separation-based node **/
/**                separation method.                      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 11 oct 2003     **/
/**                                 to     11 oct 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*
**  The type and structure definitions.
*/

/*+ Method parameters. +*/

typedef struct VmeshSeparateGrParam_ {
  Strat *                   stratptr;             /*+ Graph separation strategy +*/
} VmeshSeparateGrParam;

/*
**  The function prototypes.
*/

#ifndef VMESH_SEPARATE_GR
#define static
#endif

int                         vmeshSeparateGr     (Vmesh * restrict const, const VmeshSeparateGrParam * restrict const);

#undef static
