/*
 * $Id: sort.c,v 1.4 2000/08/10 21:02:51 danny Exp $
 *
 * Copyright � 1990, 1992, 1993 Free Software Foundation, Inc.
 * 
 * This file is part of Oleo, the GNU Spreadsheet.
 * 
 * Oleo is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 * 
 * Oleo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Oleo; see the file COPYING.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Plug-compatible replacement for UNIX qsort.
   Copyright (C) 1989,1990, 1992, 1993 Free Software Foundation, Inc.
   Written by Douglas C. Schmidt (schmidt@ics.uci.edu)
   Modified by Jay Fenlason for use in the spreadsheet. */

/* The modifications for the spreadsheet consist of calling functions to
   compare, swap, and rotate array elements.  This allows it to work
   on sparse arrays, where some or all of the elements may not exist
 */

#include "sort.h"

/* Envoke the comparison function, returns either 0, < 0, or > 0. */
#define CMP(A,B) ((*cmp)((A),(B)))

/* Byte-wise swap two items of size SIZE. */
/* #define SWAP(A,B,SIZE) do {int sz = (SIZE); char *a = (A); char *b = (B); \
    do { char _temp = *a;*a++ = *b;*b++ = _temp;} while (--sz);} while (0) */
#define SWAP(A,B) ((*swap)((A),(B)))

#define ROT(a,b) ((*rot)((a),(b)))

/* Copy SIZE bytes from item B to item A. */
/* #define COPY(A,B,SIZE) {int sz = (SIZE); do { *(A)++ = *(B)++; } while (--sz); } */

/* This should be replaced by a standard ANSI macro. */
#define BYTES_PER_WORD 8

/* The next 4 #defines implement a very fast in-line stack abstraction. */
#define STACK_SIZE (BYTES_PER_WORD * sizeof (long))
#define PUSH(LOW,HIGH) do {top->lo = LOW;top++->hi = HIGH;} while (0)
#define POP(LOW,HIGH)  do {LOW = (--top)->lo;HIGH = top->hi;} while (0)
#define STACK_NOT_EMPTY (stack < top)

/* Discontinue quicksort algorithm when partition gets below this size.
   This particular magic number was chosen to work best on a Sun 4/260. */
#define MAX_THRESH 4

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
{
  int lo;
  int hi;
}

stack_node;

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
      next array partition to sort.  To save time, this maximum amount
      of space required to store an array of MAX_INT is allocated on the
      stack.  Assuming a 32-bit integer, this needs only 32 *
      sizeof (stack_node) == 136 bits.  Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
      insertion sort to order the MAX_THRESH items within each partition.
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segements.

   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (n)
      stack size is needed (actually O(1) in this case)! */

int
sort (total_elems, cmp, swap, rot)
     int total_elems;
     int (*cmp) ();
     void (*swap) ();
     void (*rot) ();
{
  /* Allocating SIZE bytes for a pivot buffer facilitates a better
     algorithm below since we can do comparisons directly on the pivot. */
  int max_thresh = MAX_THRESH;

  if (total_elems > MAX_THRESH)
    {
      int lo = 0;
      int hi = lo + (total_elems - 1);
      stack_node stack[STACK_SIZE];	/* Largest size needed for 32-bit int!!! */
      stack_node *top = stack + 1;

      while (STACK_NOT_EMPTY)
	{
	  int left_ptr;
	  int right_ptr;
	  {
	    int pivot;
	    {
	      /* Select median value from among LO, MID, and HI. Rearrange
                 LO and HI so the three values are sorted. This lowers the
                 probability of picking a pathological pivot value and
                 skips a comparison for both the LEFT_PTR and RIGHT_PTR. */

	      int mid = lo + ((hi - lo) >> 1);

	      if (CMP (mid, lo) < 0)
		SWAP (mid, lo);
	      if (CMP (hi, mid) < 0)
		SWAP (mid, hi);
	      else
		goto jump_over;
	      if (CMP (mid, lo) < 0)
		SWAP (mid, lo);
	    jump_over:
	      /* COPY (pivot, mid); */
	      pivot = mid;
	      /* pivot = pivot_buffer; */
	    }
	    left_ptr = lo + 1;
	    right_ptr = hi - 1;

	    /* Here's the famous ``collapse the walls'' section of quicksort.
               Gotta like those tight inner loops!  They are the main reason
               that this algorithm runs much faster than others. */
	    do
	      {
		while (CMP (left_ptr, pivot) < 0)
		  left_ptr += 1;

		while (CMP (pivot, right_ptr) < 0)
		  right_ptr -= 1;

		if (left_ptr < right_ptr)
		  {
		    SWAP (left_ptr, right_ptr);
		    left_ptr += 1;
		    right_ptr -= 1;
		  }
		else if (left_ptr == right_ptr)
		  {
		    left_ptr += 1;
		    right_ptr -= 1;
		    break;
		  }
	      }
	    while (left_ptr <= right_ptr);

	  }

	  /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size. If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

	  if ((right_ptr - lo) <= max_thresh)
	    {
	      if ((hi - left_ptr) <= max_thresh)	/* Ignore both small partitions. */
		POP (lo, hi);
	      else		/* Ignore small left partition. */
		lo = left_ptr;
	    }
	  else if ((hi - left_ptr) <= max_thresh)	/* Ignore small right partition. */
	    hi = right_ptr;
	  else if ((right_ptr - lo) > (hi - left_ptr))	/* Push larger left partition indices. */
	    {
	      PUSH (lo, right_ptr);
	      lo = left_ptr;
	    }
	  else
	    /* Push larger right partition indices. */
	    {
	      PUSH (left_ptr, hi);
	      hi = right_ptr;
	    }
	}
    }

  /* Once the BASE_PTR array is partially sorted by quicksort the rest
     is completely sorted using insertion sort, since this is efficient
     for partitions below MAX_THRESH size. BASE_PTR points to the beginning
     of the array to sort, and END_PTR points at the very last element in
     the array (*not* one beyond it!). */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

  {
    int end_ptr = 0 + (total_elems - 1);
    int run_ptr;
    int tmp_ptr = 0;
    int thresh = MIN (end_ptr, 0 + max_thresh);

    /* Find smallest element in first threshold and place it at the
       array's beginning.  This is the smallest array element,
       and the operation speeds up insertion sort's inner loop. */

    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr += 1)
      if (CMP (run_ptr, tmp_ptr) < 0)
	tmp_ptr = run_ptr;

    if (tmp_ptr != 0)
      SWAP (tmp_ptr, 0);

    /* Insertion sort, running from left-hand-side up to `right-hand-side.'
       Pretty much straight out of the original GNU qsort routine. */

    for (run_ptr = 0 + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;)
      {

	while (CMP (run_ptr, tmp_ptr -= 1) < 0)
	  ;

	if ((tmp_ptr += 1) != run_ptr)
	  {
	    ROT (tmp_ptr, run_ptr);
	    /* int trav;

	    for (trav = run_ptr + 1; --trav >= run_ptr;)
              {
                char c = *trav;
                char *hi, *lo;

                for (hi = lo = trav; (lo -= 1) >= tmp_ptr; hi = lo)
                  *hi = *lo;
                *hi = c;
              } */
	  }

      }
  }
  return 1;
}


#ifdef TEST_ME
int buf[25] =
{
  1, 15, 37, 9, 100, 3, 14, 2, 88,
  6, 97, 12, 34, 8, 92, 11, 15, 38,
  16, 6, 93, 42, 45, 55, 64};
main ()
{
  int com ();
  void swa ();
  void rot ();
  int n;

  sort (25, com, swa, rot);
  for (n = 0; n < 25; n++)
    printf ("%d ", buf[n]);
}

com (n1, n2)
{
  printf ("Cmp %d,%d\n", n1, n2);
  return buf[n1] - buf[n2];
}

void
swa (n1, n2)
{
  int t;

  printf ("Swap %d,%d\n", n1, n2);
  t = buf[n1];
  buf[n1] = buf[n2];
  buf[n2] = t;
}

void
rot (n1, n2)
{
  int t;

  printf ("Rot %d,%d\n", n1, n2);
  t = buf[n2];
  while (n2 > n1)
    {
      buf[n2] = buf[n2 - 1];
      --n2;
    }
  buf[n2] = t;
}

#endif
