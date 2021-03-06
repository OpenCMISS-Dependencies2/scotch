/* This file is part of the Scotch distribution.
** It does not include the standard Scotch header because it is a
** slight adaptation of the qsort routine of glibc 2.4, taylored to
** match Scotch needs.
** Consequently, this file is distributed according to the terms of
** the GNU LGPL, see copyright notice below.
*/

/* Copyright (C) 1991,1992,1996,1997,1999,2004 Free Software Foundation, Inc.
   This file is extracted from the GNU C Library.
   Written by Douglas C. Schmidt (schmidt@ics.uci.edu).
   Modifications (C) 2019 IPB, Universite de Bordeaux, INRIA & CNRS
   Modified by Francois Pellegrini (francois.pellegrini@u-bordeaux.fr).

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/* If you consider tuning this algorithm, you should consult first:
   Engineering a sort function; Jon Bentley and M. Douglas McIlroy;
   Software - Practice and Experience; Vol. 23 (11), 1249-1265, 1993.  */

#ifndef STACK_NODE_L_DEFINED
#define STACK_NODE_L_DEFINED

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
  {
    char *lo;
    char *hi;
    int lv;
  } stack_node_l;

#endif /* STACK_NODE_L_DEFINED */

#ifndef MAX_THRESH

#define MAX_THRESH 6

#define max_thresh                  (MAX_THRESH * INTSORTSIZE) /* Variable turned into constant */

/* The next 4 #defines implement a very fast in-line stack abstraction. */
/* The stack needs log (total_elements) entries (we could even subtract
   log(MAX_THRESH)).  Since total_elements has type size_t, we get as
   upper bound for log (total_elements):
   bits per byte (CHAR_BIT) * sizeof(size_t).  */
#define STACK_SIZE	(CHAR_BIT * sizeof (INT))
#define PUSH(low, high, lvl)	((void) ((top->lo = (low)), (top->hi = (high)), (top->lv = (lvl)), ++top))
#define	POP(low, high, lvl)	((void) (--top, (low = top->lo), (high = top->hi), (lvl = top->lv)))
#define	STACK_NOT_EMPTY	(stack < top)

#endif /* MAX_THRESH */

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
      next array partition to sort.  To save time, this maximum amount
      of space required to store an array of SIZE_MAX is allocated on the
      stack.  Assuming a 32-bit (64 bit) integer for size_t, this needs
      only 32 * sizeof(stack_node) == 256 bytes (for 64 bit: 1024 bytes).
      Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
      insertion sort to order the MAX_THRESH items within each partition.
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (total_elems)
      stack size is needed (actually O(1) in this case)!  */

/* To be defined :
** INTSORTQUAL : whether the function is "static" or not
** INTSORTNAME : Name of function
** INTSORTSIZE : Size of elements to sort
** INTSORTSWAP : Swapping macro
** INTSORTCMP  : Comparison function
*/

#ifdef INTSORTQUAL
INTSORTQUAL
#endif /* INTSORTQUAL */
void
INTSORTNAME (
void * const                pbase,                /*+ Array to sort                  +*/
const INT                   total_elems,          /*+ Number of entries to sort      +*/
const int                   max_levels)           /*+ Maximum number of levels to go +*/
{
  register char *base_ptr = (char *) pbase;

  if (total_elems > MAX_THRESH)
    {
      char *lo = base_ptr;
      char *hi = &lo[INTSORTSIZE * (total_elems - 1)];
      int lvl = 0;
      stack_node_l stack[STACK_SIZE];
      stack_node_l *top = stack;

      PUSH (NULL, NULL, 0);

      while (STACK_NOT_EMPTY)
        {
          char *left_ptr;
          char *right_ptr;

	  /* Select median value from among LO, MID, and HI. Rearrange
	     LO and HI so the three values are sorted. This lowers the
	     probability of picking a pathological pivot value and
	     skips a comparison for both the LEFT_PTR and RIGHT_PTR in
	     the while loops. */

	  char *mid = lo + INTSORTSIZE * ((hi - lo) / INTSORTSIZE >> 1);

	  if (INTSORTCMP ((void *) mid, (void *) lo))
	    INTSORTSWAP (mid, lo);
	  if (INTSORTCMP ((void *) hi, (void *) mid))
	    {
	      INTSORTSWAP (mid, hi);
	      if (INTSORTCMP ((void *) mid, (void *) lo))
		INTSORTSWAP (mid, lo);
	    }

	  left_ptr  = lo + INTSORTSIZE;
	  right_ptr = hi - INTSORTSIZE;

	  /* Here's the famous ``collapse the walls'' section of quicksort.
	     Gotta like those tight inner loops!  They are the main reason
	     that this algorithm runs much faster than others. */
	  do
	    {
	      while (INTSORTCMP ((void *) left_ptr, (void *) mid))
		left_ptr += INTSORTSIZE;

	      while (INTSORTCMP ((void *) mid, (void *) right_ptr))
		right_ptr -= INTSORTSIZE;

	      if (left_ptr < right_ptr)
		{
		  INTSORTSWAP (left_ptr, right_ptr);
		  if (mid == left_ptr)
		    mid = right_ptr;
		  else if (mid == right_ptr)
		    mid = left_ptr;
		  left_ptr += INTSORTSIZE;
		  right_ptr -= INTSORTSIZE;
		}
	      else if (left_ptr == right_ptr)
		{
		  left_ptr += INTSORTSIZE;
		  right_ptr -= INTSORTSIZE;
		  break;
		}
	    }
	  while (left_ptr <= right_ptr);

          /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size.  If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

	  if (++lvl >= max_levels)
	    {
	      /* If maximum level achieved, prevent going on further on this branch */
	      lo = right_ptr;
	      hi = left_ptr;
	    }

          if ((size_t) (right_ptr - lo) <= max_thresh)
            {
              if ((size_t) (hi - left_ptr) <= max_thresh)
		/* Ignore both small partitions. */
                POP (lo, hi, lvl);
              else
		/* Ignore small left partition. */
                lo = left_ptr;
            }
          else if ((size_t) (hi - left_ptr) <= max_thresh)
	    /* Ignore small right partition. */
            hi = right_ptr;
          else if ((right_ptr - lo) > (hi - left_ptr))
            {
	      /* Push larger left partition indices. */
              PUSH (lo, right_ptr, lvl);
              lo = left_ptr;
            }
          else
            {
	      /* Push larger right partition indices. */
              PUSH (left_ptr, hi, lvl);
              hi = right_ptr;
            }
        }
    }
}

#undef MAX_THRESH
#undef max_thresh
#undef STACK_SIZE
#undef PUSH
#undef POP
#undef STACK_NOT_EMPTY
