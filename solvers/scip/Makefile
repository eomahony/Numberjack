#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.4 2008/09/29 20:33:54 bzfviger Exp $

#@file    Makefile
#@brief   Makefile for the SCIP tutorial example by Cornelius Schwarz
#@author  Marc Pfetsch


#MAINDIR = ../..
include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

ZIBDIR = $(ZIBPATH)

SCIPDIR = $(shell find $(ZIBPATH) -name "scip-?.?.?")
SOPLEXDIR = $(shell find $(ZIBPATH) -name "soplex-?.?.?")

PYTHON = ./python

ifeq ($(OSTYPE),darwin)
	ifeq ($(ARCH),x86_64)
		# Fixes an issue with gcc4.2 creating two 64bit binarys and trying to link them with `lipo`, even though it is given -arch i386 -arch x86_64 by python-config.
		export ARCHFLAGS = -arch x86_64
	endif
endif

ifeq ($(PYVERSION),python2.6)
PYLIB = -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall  -DENABLE_DTRACE  -arch x86_64 
PYOBJ = `python2.6-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.5)
PYLIB = `python2.5-config --cflags` 
PYOBJ = `python2.5-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.4)
PYLIB = `/usr/bin/python2.4-config --cflags` 
PYOBJ = `/usr/bin/python2.4-config --prefix`
SWIG = swig -classic
else 
PYLIB = `python-config --cflags`
PYOBJ = `python-config --prefix`
SWIG = swig
endif
endif
endif
#SWIG = swig
#PYLIB = `python-config --cflags` $(HACK)
ZIMPL=false

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
include ../../make.project  # We need to include this again because SCIP's make.project clears CFLAGS and we need to re-add -flat_namespace for compiling on Mac.
ZIMPL=false 
ZLIB=false 
READLINE=false

#include $(SCIPDIR)/Makefile

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	SCIP
MAINOBJ		=	SCIP.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))

MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))


#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE		=	$(PYTHON)/$(MAINNAME).i
PYFILE		=	$(PYTHON)/$(MAINNAME).py
WXXFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.o
SOFILE		= 	$(PYTHON)/_$(MAINNAME).so

MIPDIR          =       $(MAINDIR)/solvers/mip
MIPFLAGS        =       -I$(MIPDIR)/src/ 


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

NFLAGS = $(FLAGS) #-D_DEBUGWRAP

wrapper: $(SCIPDIR)/done.txt $(IFILE) $(WXXFILE) $(SOFILE)


install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install	

$(OBJDIR):
	@-mkdir -p $(OBJDIR)

clean:
	-rm -rf $(OBJDIR)/* $(PYTHON)/*o $(PYTHON)/*~ $(PYTHON)/build *~ $(SRCDIR)/*~ $(SCIPDIR)/done.txt

clean_swig: 
	-rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~ 

clean_all: clean clean_swig
	cd ./ziboptsuite/; make clean


$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(NFLAGS) $(MIPFLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -fPIC -c $< $(CXX_o)$@

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
		@echo "-> compiling $@"
		$(CXX) $(MIPFLAGS) $(NFLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -I$(PYLIB) -fPIC -c $< $(CXX_o)$@

#$(PYTHON)/%.py: $(PYTHON)/%.i
#	python tools/finalize.py $(<:.i=.py)

#$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
#	$(SWIG) -c++ -python $<
#	python $(MAINDIR)/tools/finalize.py $(<:.i=.py)

#$(PYTHON)/%.i: $(SRCDIR)/%.hpp
#	python $(MAINDIR)/tools/mk_interface.py $<

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"	
	$(SWIG) -c++ -python $<
	python $(MAINDIR)/tools/finalize.py $(<:.i=.py) MipWrapper

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	python $(MAINDIR)/tools/mk_subinterface.py mip MipWrapper $<


$(PYTHON)/_%.so: $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(MAINOBJFILES) $(OBJDIR)/%.o  $(MIPDIR)/obj/MipWrapper.o $(PYTHON)/%_wrap.o
	@echo "-> linking $@"
	$(LINKCXX) $(BUNDLE) -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX)  $(OFLAGS) $(LPSLDFLAGS) $(MIPDIR)/obj/MipWrapper.o python/SCIP_wrap.o $(LINKCXX_o)$@

#scip_lib: $(SCIPDIR)/Makefile.backup 
#	cp  $(SCIPDIR)/Makefile $(SCIPDIR)/Makefile.backup
#	cp  $(MAINDIR)/tools/Makefile.scip $(SCIPDIR)/Makefile
#	cp  $(SOPLEXDIR)/Makefile $(SOPLEXDIR)/Makefile.backup
#	cp  $(MAINDIR)/tools/Makefile.soplex $(SOPLEXDIR)/Makefile
#	cd  $(SCIPDIR)/..; make ZIMPL=false ZLIB=false READLINE=false

$(SCIPDIR)/done.txt: 
	cp $(MAINDIR)/tools/Makefile.Numberjack $(SCIPDIR)/../
	cp $(MAINDIR)/tools/Makefile.doit.Numberjack $(SCIPDIR)/../
	cp $(MAINDIR)/tools/Makefile.scip.Numberjack $(SCIPDIR)/
	cp $(MAINDIR)/tools/Makefile.soplex.Numberjack $(SOPLEXDIR)/
	cd $(SCIPDIR)/..; make -f Makefile.Numberjack ZIMPL=false ZLIB=false READLINE=false	
	@echo finished | cat > $(SCIPDIR)/done.txt

#---- EOF --------------------------------------------------------------------
