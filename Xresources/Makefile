# Generated automatically from Makefile.in by configure.
# Makefile.in generated automatically by automake 1.4 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

#
# This directory contains files to be installed as application defaults
# files for oleo. One of the purposes is to internationalise the Motif
# interface, which is not easily done with the GNU Gettext facilities in
# the po/ directory.
#
# $Header: /cvs/oleo/Xresources/Makefile.am,v 1.3 2001/02/08 00:21:23 pw Exp $
#


SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /usr/local
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/oleo
pkglibdir = $(libdir)/oleo
pkgincludedir = $(includedir)/oleo

top_builddir = ..

ACLOCAL = aclocal
AUTOCONF = autoconf
AUTOMAKE = automake
AUTOHEADER = autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
BUILD_PREFIX = /usr/local
CATALOGS = 
CATOBJEXT = 
CC = gcc
CPP = gcc -E
CXX = c++
CXXCPP = @CXXCPP@
DATADIRNAME = share
EXEEXT = 
GENCAT = 
GMOFILES =  en.gmo nl.gmo fr.gmo
GMSGFMT = 
GT_NO = 
GT_YES = #YES#
HAVE_PERL = /usr/bin/perl
HAVE_TEXI2HTML = 
INCLUDE_LOCALE_H = #include <locale.h>
INSTOBJEXT = 
INTLDEPS = 
INTLLIBS = 
INTLOBJS = 
LIBOBJS =  mktime.o
LN_S = ln -s
MAINT = #
MAKEINFO = /home/edmond/software/oleo/missing makeinfo
MKINSTALLDIRS = ./mkinstalldirs
MOTIF_CFLAGS = 
MOTIF_LIBS = -lXm
MSGFMT = 
PACKAGE = oleo
POFILES =  en.po nl.po fr.po
POSUB = 
RANLIB = ranlib
U = 
USE_INCLUDED_LIBINTL = no
USE_NLS = no
VERSION = 1.99.16
XBAE_CFLAGS = 
XBAE_LIBS = -lXbae
XLT_CFLAGS = 
XLT_LIBS = 
XPM_CFLAGS = 
XPM_LIBS = -lXpm
YACC = bison -y
include_motif = 
include_sciplot = -I/usr/include
include_xbae = 
include_xmhtml = 
l = 
link_motif = -lXm
link_sciplot = -L/usr/lib -lsciplot
link_xbae = 
link_xmhtml = -L/usr/lib/x86_64-linux-gnu -lXmHTML -lXext -ljpeg -lpng -lz

xdefaultsdir = $(prefix)/share/oleo

noinst_DATA = Oleo.nl Oleo.fr Oleo.en
EXTRA_DIST = Oleo.nl Oleo.fr Oleo.en

nldir = $(xdefaultsdir)/nl
frdir = $(xdefaultsdir)/fr
endir = $(xdefaultsdir)

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES = 
DATA =  $(noinst_DATA)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
all: all-redirect
.SUFFIXES:
$(srcdir)/Makefile.in: # Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps Xresources/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status

tags: TAGS
TAGS:


distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = Xresources

distdir: $(DISTFILES)
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pr $$/$$file $(distdir)/$$file; \
	  else \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file || :; \
	  fi; \
	done
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am:
install-exec: install-exec-am

install-data-am:
	@$(NORMAL_INSTALL)
	$(MAKE) $(AM_MAKEFLAGS) install-data-hook
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am:
uninstall: uninstall-am
all-am: Makefile $(DATA)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:


mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
	-test -z "$(MAINTAINERCLEANFILES)" || rm -f $(MAINTAINERCLEANFILES)
mostlyclean-am:  mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-generic mostlyclean-am

clean: clean-am

distclean-am:  distclean-generic clean-am

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-generic distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: tags distdir info-am info dvi-am dvi check check-am \
installcheck-am installcheck install-exec-am install-exec \
install-data-am install-data install-am install uninstall-am uninstall \
all-redirect all-am all installdirs mostlyclean-generic \
distclean-generic clean-generic maintainer-clean-generic clean \
mostlyclean distclean maintainer-clean


install-data-hook:
		$(mkinstalldirs) $(nldir) $(frdir) $(endir)
		$(INSTALL_DATA) $(srcdir)/Oleo.nl $(nldir)/Oleo
		$(INSTALL_DATA) $(srcdir)/Oleo.fr $(frdir)/Oleo
		$(INSTALL_DATA) $(srcdir)/Oleo.en $(endir)/Oleo

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
