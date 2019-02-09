/*
 * $Id: xrdb.c,v 1.8 2000/08/10 21:02:51 danny Exp $
 *
 * Deal with the X Resource Manager.
 * Copyright � 1990, 1992, 1993 Free Software Foundation.
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
 * along with this program; see the file COPYING.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Coopted for oleo -t 12/92 */
/* Written by jla, 4/90 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef	WITH_DMALLOC
#include <dmalloc.h>
#endif

#ifndef X_DISPLAY_MISSING

#ifndef alloca
#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* __GNUC__ not defined.  */
#if HAVE_ALLOCA_H
#include <alloca.h>
#else /* not HAVE_ALLOCA_H */
#if defined (MSDOS) && !defined (__TURBOC__)
#include <malloc.h>
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
#include <malloc.h>
 #pragma alloca
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not HAVE_ALLOCA_H */
#endif /* __GNUC__ not defined.  */
#endif /* alloca not defined.  */

#include <stdio.h>
#include "sysdef.h"

/* This should be included before the X include files; otherwise, we get
   warnings about redefining NULL under BSD 4.3.  */
#include <sys/param.h>
#define NeedFunctionPrototypes 1 
#include <X11/Xfuncproto.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/Xresource.h>
#ifdef VMS
#include "vms-pwd.h"
#else
#include <pwd.h>
#endif

/*
 * Already included <sys/stat.h> via sysdef.h
 */

#if !defined(S_ISDIR) && defined(S_IFDIR)
#define       S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#endif


#ifndef MAXPATHLEN
#define MAXPATHLEN	256
#endif

extern char *getenv ();
#ifndef __sgi
extern struct passwd *getpwuid ();
#endif /* !__sgi */
extern struct passwd *getpwnam ();

static char *
gethomedir (dirname)
     char *dirname;
{
  int uid;
  struct passwd *pw;
  char *ptr;

  if ((ptr = getenv ("HOME")) == NULL)
    {
      if ((ptr = getenv ("USER")) != NULL)
	pw = getpwnam (ptr);
      else
	{
	  uid = getuid ();
	  pw = getpwuid (uid);
	}
      if (pw)
	ptr = pw->pw_dir;
      else
	{
	  ptr = NULL;
	  *dirname = '\0';
	}
    }

  if (ptr != NULL) 
    strcpy (dirname, ptr);

  dirname += strlen (dirname);
  *dirname = '/';
  dirname++;
  *dirname = '\0';

  return dirname;
}

static int
file_p (path)
     char *path;
{
  struct stat status;

  return (access (path, 4) == 0			/* exists and is readable */
	  && stat (path, &status) == 0		/* get the status */
	  && !S_ISDIR(status.st_mode));
}


/* Isn't this just disgusting? */

#define X_DEFAULT_SEARCH_PATH "/usr/lib/X11/%L/%T/%N%S:/usr/lib/X11/%l/%T/%N%S:/usr/lib/X11/%T/%N%S"

static int
decode_magic (string, file, return_path)
     char *string, *file, *return_path;
{
  char *p = string;
  char *t = return_path;

  while (*p)
    {
      if (*p == '%')
	switch (*++p)
	  {
	  case '%':
	    *t++ = '%';
	    p++;
	    break;

	  case 'N':
	  case 'T':
	  case 'S':
	  case 'L':
	  case 'l':
	  case 't':
	  case 'c':
	  default:
	    p++;
	    if (*t == '/' && *p == '/')
	      p++;
	    break;
	  }
      else
	*t++ = *p++;
    }
  *t = '\0';
  strcat (return_path, file);

  if (file_p (return_path))
    return 1;

  return_path[0] = '\0';
  return 0;
}

static int
magic_searchpath_decoder (incantation_string, file, return_path)
     char *incantation_string, *return_path, *file;
{
  register char *s = incantation_string;
  register char *p;

  /* Must be big enough for "%N%S".  */
  register int string_size = MAXPATHLEN;
  register char *string = (char *) alloca (string_size * sizeof (*string));

  while (*s)
    {
      p = s;

      while (*p && *p != ':')
	p++;

      if (p == s)
	{
	  /* We know string is big enough for this.  */
	  bcopy ("%N%S", string, 5);
	  if (decode_magic (string, file, return_path))
	    return 1;
	}
      else
	{
	  int len = p - s;

	  if (string_size < len+1)
	    {
	      string_size = 2 * len;
	      string = (char *) alloca (string_size * sizeof (*string));
	    }
	  bcopy (s, string, len);
	  string[len + 1] = '\0';
	  if (decode_magic (string, file, return_path))
	    return 1;
	}

      if (*p)
	s = p + 1;
      else
	return 0;
    }

  return 0;
}

static XrmDatabase
get_system_app (class)
     char *class;
{
  XrmDatabase db;
  char path[MAXPATHLEN];
  char *p;

  if ((p = getenv ("XFILESEARCHPATH")) == NULL)
    p = X_DEFAULT_SEARCH_PATH;

  if (! magic_searchpath_decoder (p, class, path))
    return NULL;

  db = XrmGetFileDatabase (path);
  return db;
}

static XrmDatabase
get_fallback (display)
     Display *display;
{
  return NULL;
}

static XrmDatabase
get_user_app (class)
     char *class;
{
  XrmDatabase db;
  char *magic_path;
  char path[MAXPATHLEN];

  if ((magic_path = getenv ("XUSERFILESEARCHPATH")) == NULL)
    {
      char homedir[MAXPATHLEN];
      char *default_magic;
      char *p;

      gethomedir (homedir);

      if ((p = getenv ("XAPPLRESDIR")) == NULL)
	{
	  default_magic = "%s/%%L/%%N:%s/%%l/%%N:%s/%%N";
	  magic_path = (char *) alloca ((3 * strlen (homedir))
					+ strlen (default_magic));
	  sprintf (magic_path, default_magic, homedir, homedir, homedir);
	}
      else
	{
	  default_magic = "%s/%%L/%%N:%s/%%l/%%N:%s/%%N:%s/%%N";
	  magic_path = (char *) alloca ((3 * strlen (p))
					+ strlen (default_magic)
					+ strlen (homedir));
	  sprintf (magic_path, default_magic, p, p, p, homedir);
	}
    }

  if (! magic_searchpath_decoder (magic_path, class, path))
    return NULL;
  
  db = XrmGetFileDatabase (path);
  return db;
}

static XrmDatabase
get_user_db (display)
     Display *display;
{
  XrmDatabase db;
  char *xdefs;

  xdefs = XResourceManagerString (display);

  if (xdefs != NULL)
    db = XrmGetStringDatabase (xdefs);
  else
    {
      char xdefault[MAXPATHLEN];

      gethomedir (xdefault);
      strcat (xdefault, ".Xdefaults");
      db = XrmGetFileDatabase (xdefault);
    }

  return db;
}

static XrmDatabase
get_environ_db ()
{
  XrmDatabase db;
  char *p;
  char path[MAXPATHLEN];

  if ((p = getenv ("XENVIRONMENT")) == NULL)
    {
      gethomedir (path);
      strcat (path, ".Xdefaults-");
      gethostname (path + strlen (path), MAXPATHLEN - strlen (path));
      p = path;
    }

  db = XrmGetFileDatabase (p);
  return db;
}

/* Types of values that we can find in a database */

#define XrmStringType "String"	/* String representation */
XrmRepresentation x_rm_string = 0;	/* Quark representation */

/* Load X resources based on the display and a possible -xrm option. */

XrmDatabase
x_load_resources (display, xrm_string, myclass)
     Display *display;
     char *xrm_string, *myclass;
{
  XrmDatabase rdb;
  XrmDatabase db;

  if (!x_rm_string)
    x_rm_string = XrmStringToQuark (XrmStringType);

  XrmInitialize ();
  rdb = XrmGetStringDatabase ("");

  /* Get application system defaults */
  db = get_system_app (myclass);
  if (db != NULL)
    XrmMergeDatabases (db, &rdb);

  /* Get Fallback resources */
  db = get_fallback (display);
  if (db != NULL)
    XrmMergeDatabases (db, &rdb);

  /* Get application user defaults */
  db = get_user_app (myclass);
  if (db != NULL)
    XrmMergeDatabases (db, &rdb);

  /* get User defaults */
  db = get_user_db (display);
  if (db != NULL)
    XrmMergeDatabases (db, &rdb);

  /* Get Environment defaults. */
  db = get_environ_db ();
  if (db != NULL)
    XrmMergeDatabases (db, &rdb);
  
  /* Last, merge in any specification from the command line. */
  if (xrm_string != NULL)
    {
      db = XrmGetStringDatabase (xrm_string);
      if (db != NULL)
	XrmMergeDatabases (db, &rdb);
    }

  return rdb;
}

/* Retrieve the value of the resource specified by NAME with class CLASS
   and of type TYPE from database RDB.  The value is returned in RET_VALUE. */

int
x_get_resource (rdb, name, class, expected_type, ret_value)
     XrmDatabase rdb;
     char *name, *class;
     XrmRepresentation expected_type;
     XrmValue *ret_value;
{
  XrmValue value;
  XrmName namelist[100];
  XrmClass classlist[100];
  XrmRepresentation type;

  XrmStringToNameList(name, namelist);
  XrmStringToClassList(class, classlist);

  if (XrmQGetResource (rdb, namelist, classlist, &type, &value) == True
      && (type == expected_type))
    {
      if (type == x_rm_string)
	ret_value->addr = (char *) value.addr;
      else
	bcopy (value.addr, ret_value->addr, ret_value->size);

      return value.size;
    }

  return 0;
}

/* Retrieve the string resource specified by NAME with CLASS from
   database RDB. */

char *
x_get_string_resource (rdb, class, name)
     XrmDatabase rdb;
     char *name, *class;
{
  XrmValue value;

  if (!x_rm_string)
    x_rm_string = XrmStringToQuark (XrmStringType);

  if (x_get_resource (rdb, name, class, x_rm_string, &value))
    return (char *) value.addr;

  return (char *) 0;
}

#ifdef TESTRM
#include <stdio.h>
#include "arg-list.h"

static void
fatal (msg, prog, x1, x2, x3, x4, x5)
    char *msg, *prog;
    int x1, x2, x3, x4, x5;
{
    extern int errno;

    if (errno)
      perror (prog);

    (void) fprintf (stderr, msg, prog, x1, x2, x3, x4, x5);
    exit (1);
}

main (argc, argv)
    int argc;
    char **argv;
{
  Display *display;
  char *displayname, *resource_string, *class;
  XrmDatabase xdb;
  List *arg_list, *lp;

  arg_list = arg_listify (argc, argv);

  lp = member ("-d", arg_list);
  if (!NIL (lp))
    displayname = car (cdr (lp));
  else
    displayname = "localhost:0.0";

  lp = member ("-xrm", arg_list);
  if (! NIL (lp))
    resource_string = car (cdr (lp));
  else
    resource_string = (char *) 0;

  lp = member ("-c", arg_list);
  if (! NIL (lp))
    class = car (cdr (lp));
  else
    class = "Emacs";

  free_arglist (arg_list);



  if (!(display = XOpenDisplay (displayname)))
    fatal ("Can't open display '%s'\n", XDisplayName (displayname));

  xdb = x_load_resources (display, resource_string, class);

#if 0
  /* In a real program, you'd want to also do this: */
  display->db = xdb;
#endif

  while (1)
    {
      char line[90];

      printf ("String: ");
      gets (line);
      if (strlen (line))
	{
	  char *value = x_get_string_resource (xdb, line, class);

	  if (value != NULL)
	    printf ("\t%s:  %s\n\n", line, value);
	  else
	    printf ("\tNo Value.\n\n");
	}
      else
	break;
    }
  printf ("\tExit.\n\n");

  XCloseDisplay (display);
}
#endif /* TESTRM */
#endif /* X_DISPLAY_MISSING */

