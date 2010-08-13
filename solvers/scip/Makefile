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

SCIPDIR         =       /Users/eomahony/ziboptsuite-1.2.0/scip-1.2.0
#SCIPDIR         =       /Users/ehebrard/work/dev/ziboptsuite-1.2.0/scip-1.2.0
#SCIPDIR         =       /Users/eomahony/Documents/workspace/ziboptsuite-1.1.0/scip-1.1.0 
#SCIPDIR         =       /Users/eoinomahony/ProgressiveParty/ziboptsuite-1.1.0/scip-1.1.0

PYTHON = ./python

ifeq ($(PYVERSION),python2.6)
PYLIB = `python2.6-config --cflags` $(HACK) 
PYOBJ = `python2.6-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.5)
PYLIB = `python2.5-config --cflags` 
PYOBJ = `python2.5-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.4)
PYLIB = `/usr/bin/python2.4-config --cflags` $(HACK) 
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
		$(CXX) $(MIPFLAGS) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -I$(PYLIB) -fPIC -c $< $(CXX_o)$@

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
	$(LINKCXX) $(BUNDLE) -flat_namespace -Wno-long-double $(MAINOBJFILES) $(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX)  $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(MIPDIR)/obj/MipWrapper.o python/SCIP_wrap.o $(LINKCXX_o)$@

scip_lib: $(SCIPDIR)/done.txt
	cd  $(SCIPDIR)/..; make ZIMPL=false READLINE=false

$(SCIPDIR)/done.txt: 
	cd $(SCIPDIR)/..; make ZIMPL=false READLINE=false	
	@echo finished | cat > $(SCIPDIR)/done.txt

#---- EOF --------------------------------------------------------------------
