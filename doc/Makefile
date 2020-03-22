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
# $Header: /cvs/oleo/doc/Makefile.am,v 1.15 2001/02/14 20:54:49 danny Exp $
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
MAKEINFO = /home/edmond/software/oleo-1.99.16/missing makeinfo
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
include_sciplot = 
include_xbae = 
include_xmhtml = 
l = 
link_motif = -lXm
link_sciplot = 
link_xbae = 
link_xmhtml = -L/usr/lib/x86_64-linux-gnu -lXmHTML -lXext -ljpeg -lpng -lz

info_TEXINFOS = oleo.texi

TEXIS = basics.texi bugs.texi cell-refs.texi entering-data.texi 		expressions.texi functions.texi interface.texi keymaps.texi 		macros.texi making-changes.texi options.texi print-plot.texi 		programs.texi ranlib.texi style.texi texinfo.tex windows.texi 		motif.texi dbms.texi


EXTRA_DIST = $(TEXIS) oleo_web.html texi2html oleo1.png
noinst_DATA = oleo_web.html $(TEXIS)

htmldir = $(prefix)/share/oleo
html_DATA = oleo.html

DISTCLEANFILES = $(html_DATA)
MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES =  texi2html
TEXI2DVI = texi2dvi
INFO_DEPS = oleo.info
DVIS = oleo.dvi
TEXINFOS = oleo.texi
DATA =  $(html_DATA) $(noinst_DATA)

DIST_COMMON =  Makefile.am Makefile.in texi2html.in texinfo.tex


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
all: all-redirect
.SUFFIXES:
.SUFFIXES: .dvi .info .ps .texi .texinfo .txi
$(srcdir)/Makefile.in: # Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps doc/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status

texi2html: $(top_builddir)/config.status texi2html.in
	cd $(top_builddir) && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status

oleo.info: oleo.texi
oleo.dvi: oleo.texi


DVIPS = dvips

.texi.info:
	@cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
	cd $(srcdir) \
	  && $(MAKEINFO) `echo $< | sed 's,.*/,,'`

.texi.dvi:
	TEXINPUTS=.:$$TEXINPUTS \
	  MAKEINFO='$(MAKEINFO) -I $(srcdir)' $(TEXI2DVI) $<

.texi:
	@cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
	cd $(srcdir) \
	  && $(MAKEINFO) `echo $< | sed 's,.*/,,'`

.texinfo.info:
	@cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
	cd $(srcdir) \
	  && $(MAKEINFO) `echo $< | sed 's,.*/,,'`

.texinfo:
	@cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
	cd $(srcdir) \
	  && $(MAKEINFO) `echo $< | sed 's,.*/,,'`

.texinfo.dvi:
	TEXINPUTS=.:$$TEXINPUTS \
	  MAKEINFO='$(MAKEINFO) -I $(srcdir)' $(TEXI2DVI) $<

.txi.info:
	@cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
	cd $(srcdir) \
	  && $(MAKEINFO) `echo $< | sed 's,.*/,,'`

.txi.dvi:
	TEXINPUTS=.:$$TEXINPUTS \
	  MAKEINFO='$(MAKEINFO) -I $(srcdir)' $(TEXI2DVI) $<

.txi:
	@cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
	cd $(srcdir) \
	  && $(MAKEINFO) `echo $< | sed 's,.*/,,'`
.dvi.ps:
	$(DVIPS) $< -o $@

install-info-am: $(INFO_DEPS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(infodir)
	@list='$(INFO_DEPS)'; \
	for file in $$list; do \
	  d=$(srcdir); \
	  for ifile in `cd $$d && echo $$file $$file-[0-9] $$file-[0-9][0-9]`; do \
	    if test -f $$d/$$ifile; then \
	      echo " $(INSTALL_DATA) $$d/$$ifile $(DESTDIR)$(infodir)/$$ifile"; \
	      $(INSTALL_DATA) $$d/$$ifile $(DESTDIR)$(infodir)/$$ifile; \
	    else : ; fi; \
	  done; \
	done
	@$(POST_INSTALL)
	@if $(SHELL) -c 'install-info --version | sed 1q | fgrep -s -v -i debian' >/dev/null 2>&1; then \
	  list='$(INFO_DEPS)'; \
	  for file in $$list; do \
	    echo " install-info --info-dir=$(DESTDIR)$(infodir) $(DESTDIR)$(infodir)/$$file";\
	    install-info --info-dir=$(DESTDIR)$(infodir) $(DESTDIR)$(infodir)/$$file || :;\
	  done; \
	else : ; fi

uninstall-info:
	$(PRE_UNINSTALL)
	@if $(SHELL) -c 'install-info --version | sed 1q | fgrep -s -v -i debian' >/dev/null 2>&1; then \
	  ii=yes; \
	else ii=; fi; \
	list='$(INFO_DEPS)'; \
	for file in $$list; do \
	  test -z "$ii" \
	    || install-info --info-dir=$(DESTDIR)$(infodir) --remove $$file; \
	done
	@$(NORMAL_UNINSTALL)
	list='$(INFO_DEPS)'; \
	for file in $$list; do \
	  (cd $(DESTDIR)$(infodir) && rm -f $$file $$file-[0-9] $$file-[0-9][0-9]); \
	done

dist-info: $(INFO_DEPS)
	list='$(INFO_DEPS)'; \
	for base in $$list; do \
	  d=$(srcdir); \
	  for file in `cd $$d && eval echo $$base*`; do \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file; \
	  done; \
	done

mostlyclean-aminfo:
	-rm -f oleo.aux oleo.cp oleo.cps oleo.dvi oleo.fn oleo.fns oleo.ky \
	  oleo.kys oleo.ps oleo.log oleo.pg oleo.toc oleo.tp oleo.tps \
	  oleo.vr oleo.vrs oleo.op oleo.tr oleo.cv oleo.cn

clean-aminfo:

distclean-aminfo:

maintainer-clean-aminfo:
	cd $(srcdir) && for i in $(INFO_DEPS); do \
	  rm -f $$i; \
	  if test "`echo $$i-[0-9]*`" != "$$i-[0-9]*"; then \
	    rm -f $$i-[0-9]*; \
	  fi; \
	done

install-htmlDATA: $(html_DATA)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(htmldir)
	@list='$(html_DATA)'; for p in $$list; do \
	  if test -f $(srcdir)/$$p; then \
	    echo " $(INSTALL_DATA) $(srcdir)/$$p $(DESTDIR)$(htmldir)/$$p"; \
	    $(INSTALL_DATA) $(srcdir)/$$p $(DESTDIR)$(htmldir)/$$p; \
	  else if test -f $$p; then \
	    echo " $(INSTALL_DATA) $$p $(DESTDIR)$(htmldir)/$$p"; \
	    $(INSTALL_DATA) $$p $(DESTDIR)$(htmldir)/$$p; \
	  fi; fi; \
	done

uninstall-htmlDATA:
	@$(NORMAL_UNINSTALL)
	list='$(html_DATA)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(htmldir)/$$p; \
	done
tags: TAGS
TAGS:


distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = doc

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
	$(MAKE) $(AM_MAKEFLAGS) top_distdir="$(top_distdir)" distdir="$(distdir)" dist-info
info-am: $(INFO_DEPS)
info: info-am
dvi-am: $(DVIS)
dvi: dvi-am
check-am: all-am
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am:
install-exec: install-exec-am

install-data-am: install-info-am install-htmlDATA
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-info uninstall-htmlDATA
uninstall: uninstall-am
all-am: Makefile $(INFO_DEPS) $(DATA)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(infodir) $(DESTDIR)$(htmldir)


mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*
	-test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)

maintainer-clean-generic:
	-test -z "$(MAINTAINERCLEANFILES)" || rm -f $(MAINTAINERCLEANFILES)
mostlyclean-am:  mostlyclean-aminfo mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-aminfo clean-generic mostlyclean-am

clean: clean-am

distclean-am:  distclean-aminfo distclean-generic clean-am

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-aminfo maintainer-clean-generic \
		distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: install-info-am uninstall-info mostlyclean-aminfo \
distclean-aminfo clean-aminfo maintainer-clean-aminfo \
uninstall-htmlDATA install-htmlDATA tags distdir info-am info dvi-am \
dvi check check-am installcheck-am installcheck install-exec-am \
install-exec install-data-am install-data install-am install \
uninstall-am uninstall all-redirect all-am all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


oleo.info:	oleo.texi $(TEXIS)

oleo.html:	$(srcdir)/oleo.texi ${TEXIS} $(top_builddir)/doc/texi2html
		/usr/bin/perl $(top_builddir)/doc/texi2html $(srcdir)/oleo.texi

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
