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
# $Header: /cvs/oleo/afm/Makefile.am,v 1.5 2001/02/08 00:21:23 pw Exp $
#
# Process this file with automake to produce Makefile.in
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

afmdir = $(prefix)/share/afm

afm_DATA = Courier-Bold.afm Courier-BoldOblique.afm 		Courier-Oblique.afm Courier.afm 		Helvetica-Bold.afm Helvetica-BoldOblique.afm 		Helvetica-Oblique.afm Helvetica.afm 		Symbol.afm 		Times-Bold.afm Times-BoldItalic.afm 		Times-Italic.afm Times-Roman.afm


EXTRA_DIST = $(afm_DATA)

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES = 
DATA =  $(afm_DATA)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
all: all-redirect
.SUFFIXES:
$(srcdir)/Makefile.in: # Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps afm/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


install-afmDATA: $(afm_DATA)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(afmdir)
	@list='$(afm_DATA)'; for p in $$list; do \
	  if test -f $(srcdir)/$$p; then \
	    echo " $(INSTALL_DATA) $(srcdir)/$$p $(DESTDIR)$(afmdir)/$$p"; \
	    $(INSTALL_DATA) $(srcdir)/$$p $(DESTDIR)$(afmdir)/$$p; \
	  else if test -f $$p; then \
	    echo " $(INSTALL_DATA) $$p $(DESTDIR)$(afmdir)/$$p"; \
	    $(INSTALL_DATA) $$p $(DESTDIR)$(afmdir)/$$p; \
	  fi; fi; \
	done

uninstall-afmDATA:
	@$(NORMAL_UNINSTALL)
	list='$(afm_DATA)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(afmdir)/$$p; \
	done
tags: TAGS
TAGS:


distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = afm

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

install-data-am: install-afmDATA
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-afmDATA
uninstall: uninstall-am
all-am: Makefile $(DATA)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(afmdir)


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

.PHONY: uninstall-afmDATA install-afmDATA tags distdir info-am info \
dvi-am dvi check check-am installcheck-am installcheck install-exec-am \
install-exec install-data-am install-data install-am install \
uninstall-am uninstall all-redirect all-am all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
