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
# $Header: /cvs/oleo/src/Makefile.am,v 1.56 2001/02/08 00:26:36 pw Exp $
#
# CPPFLAGS = -Wall


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

localedir = $(prefix)/share/locale

AUTOMAKE_OPTIONS = ../lib/ansi2knr dejagnu

EXTRA_DIST = xbase.cpp oleo_xb.h mysql.c parse.h gpl.sed
CLEANFILES = parse.c parse.h posixtm.c posixtm.h
DISTCLEANFILES = tags
MAINTAINERCLEANFILES = Makefile.in

YFLAGS = -d

bin_PROGRAMS = oleo
GUI = io-motif.c appres.c fallback.c oleo_icon.xpm
#GUI = 

#
# Temporary presence of GNU Plotutils Sciplot here.
# (until it gets released as part of plotutils)
#
SCIPLOT_SRC = plotter.c sciplot.h sciplotI.h legend.c

MYSQL = mysql.c oleosql.h mysql.h
#XBSRC = xbase.cpp oleo_xb.h
XBSRC = 
#CXXFLAGS = $(CFLAGS)
#oleo_LINK = $(CXXLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
oleo_LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@

oleo_SOURCES = parse.y byte-compile.c eval.c ref.c  decompile.c sort.c regions.c utils.c lists.c  io-term.c io-utils.c io-x11.c window.c io-edit.c userpref.h  line.c key.c format.c funcs.c  hash.c ir.c io-curses.c font.c display.c init.c  sylk.c oleofile.c sc.c list.c busi.c string.c cells.c  cmd.c basic.c args.c  stub.c input.c info.c help.c  graph.c plot.c  forminfo.c posixtm.y date.c xrdb.c  $(GUI) $(XBSRC) $(MYSQL)  print.c postscript.c pcl.c epson.c prtext.c afm.c  mdi.c gsl.c  transc_funs.c ${SCIPLOT_SRC}


# MYLD = @GTK_LIBS@
# MYINC =
# MYCFL = @GTK_CFLAGS@

#
# Allow non-Motif X version to get compiled.
#
MYLD = -lXm -L/usr/lib/x86_64-linux-gnu -lXmHTML -lXext -ljpeg -lpng -lz -lXbae -L/usr/lib -lsciplot 
#MYLD = ${X_LIBS}
MYINC =    -I/usr/include 
#MYINC = 
MYCFL = ${X_CFLAGS}
#MYCFL = ${X_CFLAGS}

oleo_LDADD = ${X_PRE_LIBS} ${MYLD}  ../lib/liboleo.a ../intl/libintl.a ${X_LIBS}
INCLUDES = -I.. -I$(srcdir) -I$(top_srcdir)/lib -I../intl ${MYINC}
CFLAGS =   -g -O2 -DLOCALEDIR=\"$(localedir)\" ${MYCFL} 	-DBUILD_PREFIX='"'${BUILD_PREFIX}'"'


noinst_HEADERS = args.h basic.h byte-compile.h cell.h cmd.h decompile.h  defun.h defuns.h display.h errors.h eval.h font.h format.h forminfo.h  funcdef.h funcs.h global.h hash.h help.h info.h init.h  input.h  io-abstract.h io-curses.h io-edit.h io-generic.h io-term.h  io-utils.h io-x11.h  ir.h key.h line.h lists.h node.h oleofile.h  graph.h oleo_plot.h  proto.h ref.h regions.h  sc.h stub.h stubs.h sylk.h  sysdef.h utils.h window.h list.h ir.h sort.h  io-motif.h appres.h  print.h pcl.h postscript.h epson.h prtext.h afm.h  mdi.h

mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES = 
bin_PROGRAMS =  oleo$(EXEEXT)
PROGRAMS =  $(bin_PROGRAMS)


DEFS = -DHAVE_CONFIG_H -I. -I$(srcdir) -I..
CPPFLAGS = 
LDFLAGS = 
LIBS = -lXpm -lXm   -lSM -lICE -lXt -lX11  -lm -lncurses -L/usr/local/lib   -lplot
X_CFLAGS = 
X_LIBS = 
X_EXTRA_LIBS = 
X_PRE_LIBS =  -lSM -lICE
ANSI2KNR = ../lib/ansi2knr
#oleo_OBJECTS =  parse$U.o \
#byte-compile$U.o eval$U.o ref$U.o \
#decompile$U.o sort$U.o regions$U.o \
#utils$U.o lists$U.o io-term$U.o \
#io-utils$U.o io-x11$U.o window$U.o \
#io-edit$U.o line$U.o key$U.o format$U.o \
#funcs$U.o hash$U.o ir$U.o io-curses$U.o \
#font$U.o display$U.o init$U.o sylk$U.o \
#oleofile$U.o sc$U.o list$U.o busi$U.o \
#string$U.o cells$U.o cmd$U.o basic$U.o \
#args$U.o stub$U.o input$U.o info$U.o \
#help$U.o graph$U.o plot$U.o forminfo$U.o \
#posixtm$U.o date$U.o xrdb$U.o \
#io-motif$U.o appres$U.o fallback$U.o \
#xbase$U.o mysql$U.o print$U.o \
#postscript$U.o pcl$U.o epson$U.o \
#prtext$U.o afm$U.o mdi$U.o gsl$U.o \
#plotter$U.o legend$U.o transc_funs$U.o
oleo_OBJECTS =  parse$U.o \
byte-compile$U.o eval$U.o ref$U.o \
decompile$U.o sort$U.o regions$U.o \
utils$U.o lists$U.o io-term$U.o \
io-utils$U.o io-x11$U.o window$U.o \
io-edit$U.o line$U.o key$U.o format$U.o \
funcs$U.o hash$U.o ir$U.o io-curses$U.o \
font$U.o display$U.o init$U.o sylk$U.o \
oleofile$U.o sc$U.o list$U.o busi$U.o \
string$U.o cells$U.o cmd$U.o basic$U.o \
args$U.o stub$U.o input$U.o info$U.o \
help$U.o graph$U.o plot$U.o \
forminfo$U.o posixtm$U.o date$U.o \
xrdb$U.o io-motif$U.o appres$U.o \
fallback$U.o mysql$U.o print$U.o \
postscript$U.o pcl$U.o epson$U.o \
prtext$U.o afm$U.o mdi$U.o gsl$U.o \
plotter$U.o legend$U.o transc_funs$U.o
##oleo_OBJECTS =  parse$U.o \
##byte-compile$U.o eval$U.o ref$U.o \
##decompile$U.o sort$U.o regions$U.o \
##utils$U.o lists$U.o io-term$U.o \
##io-utils$U.o io-x11$U.o window$U.o \
##io-edit$U.o line$U.o key$U.o format$U.o \
##funcs$U.o hash$U.o ir$U.o io-curses$U.o \
##font$U.o display$U.o init$U.o sylk$U.o \
##oleofile$U.o sc$U.o list$U.o busi$U.o \
##string$U.o cells$U.o cmd$U.o basic$U.o \
##args$U.o stub$U.o input$U.o info$U.o \
##help$U.o graph$U.o plot$U.o \
##forminfo$U.o posixtm$U.o date$U.o \
##xrdb$U.o xbase$U.o mysql$U.o print$U.o \
##postscript$U.o pcl$U.o epson$U.o \
##prtext$U.o afm$U.o mdi$U.o gsl$U.o \
##plotter$U.o legend$U.o transc_funs$U.o
#oleo_OBJECTS =  parse$U.o \
#byte-compile$U.o eval$U.o ref$U.o \
#decompile$U.o sort$U.o regions$U.o \
#utils$U.o lists$U.o io-term$U.o \
#io-utils$U.o io-x11$U.o window$U.o \
#io-edit$U.o line$U.o key$U.o \
#format$U.o funcs$U.o hash$U.o ir$U.o \
#io-curses$U.o font$U.o display$U.o \
#init$U.o sylk$U.o oleofile$U.o sc$U.o \
#list$U.o busi$U.o string$U.o cells$U.o \
#cmd$U.o basic$U.o args$U.o stub$U.o \
#input$U.o info$U.o help$U.o graph$U.o \
#plot$U.o forminfo$U.o posixtm$U.o \
#date$U.o xrdb$U.o mysql$U.o print$U.o \
#postscript$U.o pcl$U.o epson$U.o \
#prtext$U.o afm$U.o mdi$U.o gsl$U.o \
#plotter$U.o legend$U.o transc_funs$U.o
oleo_DEPENDENCIES =  ../lib/liboleo.a ../intl/libintl.a
#oleo_DEPENDENCIES =  ../lib/liboleo.a ../intl/libintl.a
oleo_LDFLAGS = 
CXXCOMPILE = $(CXX) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
CXXLD = $(CXX)
CXXLINK = $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@
YLWRAP = $(srcdir)/ylwrap
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
HEADERS =  $(noinst_HEADERS)

DIST_COMMON =  Makefile.am Makefile.in parse.c posixtm.c ylwrap


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
EXPECT = expect
RUNTEST = runtest
SOURCES = $(oleo_SOURCES)
OBJECTS = $(oleo_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .cpp .o .s .y
$(srcdir)/Makefile.in: # Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps src/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


mostlyclean-binPROGRAMS:

clean-binPROGRAMS:
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)

distclean-binPROGRAMS:

maintainer-clean-binPROGRAMS:

install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo "  $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`"; \
	     $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	  else :; fi; \
	done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	list='$(bin_PROGRAMS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	done

.c.o:
	$(COMPILE) -c $<

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

mostlyclean-compile:
	-rm -f *.o core *.core

clean-compile:

distclean-compile:
	-rm -f *.tab.c

maintainer-clean-compile:
../lib/ansi2knr: ../lib/ansi2knr.o
	cd ../lib && $(MAKE) $(AM_MAKEFLAGS) ansi2knr

../lib/ansi2knr.o:
	cd ../lib && $(MAKE) $(AM_MAKEFLAGS) ansi2knr.o


mostlyclean-kr:
	-rm -f *_.c

clean-kr:

distclean-kr:

maintainer-clean-kr:

oleo$(EXEEXT): $(oleo_OBJECTS) $(oleo_DEPENDENCIES)
	@rm -f oleo$(EXEEXT)
	$(oleo_LINK) $(oleo_LDFLAGS) $(oleo_OBJECTS) $(oleo_LDADD) $(LIBS)
afm_.c: afm.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/afm.c; then echo $(srcdir)/afm.c; else echo afm.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > afm_.c
appres_.c: appres.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/appres.c; then echo $(srcdir)/appres.c; else echo appres.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > appres_.c
args_.c: args.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/args.c; then echo $(srcdir)/args.c; else echo args.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > args_.c
basic_.c: basic.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/basic.c; then echo $(srcdir)/basic.c; else echo basic.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > basic_.c
busi_.c: busi.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/busi.c; then echo $(srcdir)/busi.c; else echo busi.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > busi_.c
byte-compile_.c: byte-compile.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/byte-compile.c; then echo $(srcdir)/byte-compile.c; else echo byte-compile.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > byte-compile_.c
cells_.c: cells.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/cells.c; then echo $(srcdir)/cells.c; else echo cells.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > cells_.c
cmd_.c: cmd.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/cmd.c; then echo $(srcdir)/cmd.c; else echo cmd.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > cmd_.c
date_.c: date.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/date.c; then echo $(srcdir)/date.c; else echo date.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > date_.c
decompile_.c: decompile.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/decompile.c; then echo $(srcdir)/decompile.c; else echo decompile.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > decompile_.c
display_.c: display.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/display.c; then echo $(srcdir)/display.c; else echo display.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > display_.c
epson_.c: epson.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/epson.c; then echo $(srcdir)/epson.c; else echo epson.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > epson_.c
eval_.c: eval.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/eval.c; then echo $(srcdir)/eval.c; else echo eval.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > eval_.c
fallback_.c: fallback.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/fallback.c; then echo $(srcdir)/fallback.c; else echo fallback.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > fallback_.c
font_.c: font.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/font.c; then echo $(srcdir)/font.c; else echo font.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > font_.c
format_.c: format.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/format.c; then echo $(srcdir)/format.c; else echo format.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > format_.c
forminfo_.c: forminfo.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/forminfo.c; then echo $(srcdir)/forminfo.c; else echo forminfo.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > forminfo_.c
funcs_.c: funcs.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/funcs.c; then echo $(srcdir)/funcs.c; else echo funcs.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > funcs_.c
graph_.c: graph.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/graph.c; then echo $(srcdir)/graph.c; else echo graph.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > graph_.c
gsl_.c: gsl.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/gsl.c; then echo $(srcdir)/gsl.c; else echo gsl.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > gsl_.c
hash_.c: hash.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/hash.c; then echo $(srcdir)/hash.c; else echo hash.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > hash_.c
help_.c: help.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/help.c; then echo $(srcdir)/help.c; else echo help.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > help_.c
info_.c: info.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/info.c; then echo $(srcdir)/info.c; else echo info.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > info_.c
init_.c: init.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/init.c; then echo $(srcdir)/init.c; else echo init.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > init_.c
input_.c: input.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/input.c; then echo $(srcdir)/input.c; else echo input.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > input_.c
io-curses_.c: io-curses.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/io-curses.c; then echo $(srcdir)/io-curses.c; else echo io-curses.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > io-curses_.c
io-edit_.c: io-edit.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/io-edit.c; then echo $(srcdir)/io-edit.c; else echo io-edit.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > io-edit_.c
io-motif_.c: io-motif.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/io-motif.c; then echo $(srcdir)/io-motif.c; else echo io-motif.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > io-motif_.c
io-term_.c: io-term.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/io-term.c; then echo $(srcdir)/io-term.c; else echo io-term.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > io-term_.c
io-utils_.c: io-utils.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/io-utils.c; then echo $(srcdir)/io-utils.c; else echo io-utils.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > io-utils_.c
io-x11_.c: io-x11.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/io-x11.c; then echo $(srcdir)/io-x11.c; else echo io-x11.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > io-x11_.c
ir_.c: ir.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/ir.c; then echo $(srcdir)/ir.c; else echo ir.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > ir_.c
key_.c: key.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/key.c; then echo $(srcdir)/key.c; else echo key.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > key_.c
legend_.c: legend.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/legend.c; then echo $(srcdir)/legend.c; else echo legend.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > legend_.c
line_.c: line.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/line.c; then echo $(srcdir)/line.c; else echo line.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > line_.c
list_.c: list.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/list.c; then echo $(srcdir)/list.c; else echo list.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > list_.c
lists_.c: lists.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/lists.c; then echo $(srcdir)/lists.c; else echo lists.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > lists_.c
mdi_.c: mdi.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/mdi.c; then echo $(srcdir)/mdi.c; else echo mdi.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > mdi_.c
mysql_.c: mysql.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/mysql.c; then echo $(srcdir)/mysql.c; else echo mysql.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > mysql_.c
oleofile_.c: oleofile.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/oleofile.c; then echo $(srcdir)/oleofile.c; else echo oleofile.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > oleofile_.c
parse_.c: parse.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/parse.c; then echo $(srcdir)/parse.c; else echo parse.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > parse_.c
pcl_.c: pcl.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/pcl.c; then echo $(srcdir)/pcl.c; else echo pcl.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > pcl_.c
plot_.c: plot.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/plot.c; then echo $(srcdir)/plot.c; else echo plot.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > plot_.c
plotter_.c: plotter.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/plotter.c; then echo $(srcdir)/plotter.c; else echo plotter.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > plotter_.c
posixtm_.c: posixtm.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/posixtm.c; then echo $(srcdir)/posixtm.c; else echo posixtm.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > posixtm_.c
postscript_.c: postscript.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/postscript.c; then echo $(srcdir)/postscript.c; else echo postscript.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > postscript_.c
print_.c: print.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/print.c; then echo $(srcdir)/print.c; else echo print.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > print_.c
prtext_.c: prtext.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/prtext.c; then echo $(srcdir)/prtext.c; else echo prtext.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > prtext_.c
ref_.c: ref.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/ref.c; then echo $(srcdir)/ref.c; else echo ref.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > ref_.c
regions_.c: regions.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/regions.c; then echo $(srcdir)/regions.c; else echo regions.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > regions_.c
sc_.c: sc.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/sc.c; then echo $(srcdir)/sc.c; else echo sc.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > sc_.c
sort_.c: sort.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/sort.c; then echo $(srcdir)/sort.c; else echo sort.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > sort_.c
string_.c: string.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/string.c; then echo $(srcdir)/string.c; else echo string.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > string_.c
stub_.c: stub.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/stub.c; then echo $(srcdir)/stub.c; else echo stub.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > stub_.c
sylk_.c: sylk.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/sylk.c; then echo $(srcdir)/sylk.c; else echo sylk.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > sylk_.c
utils_.c: utils.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/utils.c; then echo $(srcdir)/utils.c; else echo utils.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > utils_.c
window_.c: window.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/window.c; then echo $(srcdir)/window.c; else echo window.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > window_.c
xrdb_.c: xrdb.c $(ANSI2KNR)
	$(CPP) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) `if test -f $(srcdir)/xrdb.c; then echo $(srcdir)/xrdb.c; else echo xrdb.c; fi` | sed 's/^# \([0-9]\)/#line \1/' | $(ANSI2KNR) > xrdb_.c
afm_.o appres_.o args_.o basic_.o busi_.o byte-compile_.o cells_.o \
cmd_.o date_.o decompile_.o display_.o epson_.o eval_.o fallback_.o \
font_.o format_.o forminfo_.o funcs_.o graph_.o gsl_.o hash_.o help_.o \
info_.o init_.o input_.o io-curses_.o io-edit_.o io-motif_.o io-term_.o \
io-utils_.o io-x11_.o ir_.o key_.o legend_.o line_.o list_.o lists_.o \
mdi_.o mysql_.o oleofile_.o parse_.o pcl_.o plot_.o plotter_.o \
posixtm_.o postscript_.o print_.o prtext_.o ref_.o regions_.o sc_.o \
sort_.o string_.o stub_.o sylk_.o utils_.o window_.o xrdb_.o : \
$(ANSI2KNR)
.cpp.o:
	$(CXXCOMPILE) -c $<
.y.c:
	$(SHELL) $(YLWRAP) "$(YACC)" $< y.tab.c $*.c y.tab.h $*.h -- $(AM_YFLAGS) $(YFLAGS)
parse.h: parse.c
posixtm.h: posixtm.c


tags: TAGS

ID: $(HEADERS) $(SOURCES) $(LISP)
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	here=`pwd` && cd $(srcdir) \
	  && mkid -f$$here/ID $$unique $(LISP)

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(ETAGS_ARGS)$$unique$(LISP)$$tags" \
	  || (cd $(srcdir) && etags $(ETAGS_ARGS) $$tags  $$unique $(LISP) -o $$here/TAGS)

mostlyclean-tags:

clean-tags:

distclean-tags:
	-rm -f TAGS ID

maintainer-clean-tags:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = src

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

RUNTESTFLAGS =

DEJATOOL = $(PACKAGE)

RUNTESTDEFAULTFLAGS = --tool $(DEJATOOL) --srcdir $$srcdir

check-DEJAGNU: site.exp
	srcdir=`cd $(srcdir) && pwd`; export srcdir; \
	EXPECT=$(EXPECT); export EXPECT; \
	runtest=$(RUNTEST); \
	if $(SHELL) -c "$$runtest --version" > /dev/null 2>&1; then \
	  $$runtest $(RUNTESTDEFAULTFLAGS) $(RUNTESTFLAGS); \
	else echo "WARNING: could not find \`runtest'" 1>&2; :;\
	fi
site.exp: Makefile
	@echo 'Making a new site.exp file...'
	@test ! -f site.bak || rm -f site.bak
	@echo '## these variables are automatically generated by make ##' > $@-t
	@echo '# Do not edit here.  If you wish to override these values' >> $@-t
	@echo '# edit the last section' >> $@-t
	@echo 'set tool $(DEJATOOL)' >> $@-t
	@echo 'set srcdir $(srcdir)' >> $@-t
	@echo 'set objdir' `pwd` >> $@-t
	@echo '## All variables above are generated by configure. Do Not Edit ##' >> $@-t
	@test ! -f site.exp || sed '1,/^## All variables above are.*##/ d' site.exp >> $@-t
	@test ! -f site.exp || mv site.exp site.bak
	@mv $@-t site.exp
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
	$(MAKE) $(AM_MAKEFLAGS) check-DEJAGNU
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am: install-binPROGRAMS
install-exec: install-exec-am

install-data-am:
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-binPROGRAMS
uninstall: uninstall-am
all-am: Makefile $(PROGRAMS) $(HEADERS)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(bindir)


mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*
	-test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)

maintainer-clean-generic:
	-test -z "parsehparsecposixtmhposixtmc$(MAINTAINERCLEANFILES)" || rm -f parseh parsec posixtmh posixtmc $(MAINTAINERCLEANFILES)
mostlyclean-am:  mostlyclean-binPROGRAMS mostlyclean-compile \
		mostlyclean-kr mostlyclean-tags mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-binPROGRAMS clean-compile clean-kr clean-tags \
		clean-generic mostlyclean-am

clean: clean-am

distclean-am:  distclean-binPROGRAMS distclean-compile distclean-kr \
		distclean-tags distclean-generic clean-am

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-binPROGRAMS \
		maintainer-clean-compile maintainer-clean-kr \
		maintainer-clean-tags maintainer-clean-generic \
		distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: mostlyclean-binPROGRAMS distclean-binPROGRAMS clean-binPROGRAMS \
maintainer-clean-binPROGRAMS uninstall-binPROGRAMS install-binPROGRAMS \
mostlyclean-compile distclean-compile clean-compile \
maintainer-clean-compile mostlyclean-kr distclean-kr clean-kr \
maintainer-clean-kr tags mostlyclean-tags distclean-tags clean-tags \
maintainer-clean-tags distdir check-DEJAGNU info-am info dvi-am dvi \
check check-am installcheck-am installcheck install-exec-am \
install-exec install-data-am install-data install-am install \
uninstall-am uninstall all-redirect all-am all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


gpl.h:	$(top_srcdir)/COPYING $(srcdir)/gpl.sed
	sed -f $(srcdir)/gpl.sed <$(top_srcdir)/COPYING >gpl.h

io-motif.o:	gpl.h $(top_srcdir)/configure $(srcdir)/oleo_icon.xpm ../config.h \
		oleo_plot.h

io-term.o:	../config.h
utils.o:	../config.h
io-utils.o:	../config.h
plot.o:		../config.h
io-term.o:	../config.h
cmd.o:		../config.h
io-x11.o:	../config.h
lists.o:	../config.h
afm.o:		../config.h $(srcdir)/afm.h
cells.o:	../config.h
print.o:	$(srcdir)/afm.h

ref.o:		parse.h
sc.o:		parse.h

plot.o:		oleo_plot.h ../config.h
graph.o:	oleo_plot.h

cmd.o funcs.o:	defun.h defuns.h

cmd.o stub.o:	stubs.h

io-term.o io-motif.o:	mdi.h
mdi.o:		../config.h

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
