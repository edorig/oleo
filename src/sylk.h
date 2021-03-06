/*
 * $Id: sylk.h,v 1.5 2000/08/10 21:02:51 danny Exp $
 *
 * Copyright � 1992, 1993, 1999 Free Software Foundation, Inc.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef SYLKH
#define SYLKH

extern void sylk_read_file (FILE *,int);
extern void sylk_write_file (FILE *,struct rng *);
extern  int sylk_set_options (int, char *);
extern void sylk_show_options (void);

#endif

