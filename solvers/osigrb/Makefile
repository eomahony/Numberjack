
include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

OBJDIR = $(CURDIR)/obj
SRCDIR = $(CURDIR)/src
PYTHON = $(CURDIR)/python
SWIG = swig
PYLIB = `python-config --cflags`

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	OsiGrb
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

MIPDIR = $(MAINDIR)/solvers/mip
MIPFLAGS = -I$(MIPDIR)/src/

OSIDIR = $(MAINDIR)/solvers/osi
OSIFLAGS = -I$(OSIDIR)/src

# Be careful with spaces here.
OSIVER     = 0.105.2
OSISRC     = $(CURDIR)/Osi-$(OSIVER)
OSIINCS    = -I$(OSISRC)/Osi/src/Osi
OSIGRBINCS = -I$(OSISRC)/Osi/src/OsiGrb
COININCS   = -I$(OSISRC)/CoinUtils/src/

OSILIBS     = $(OSISRC)/Osi/src/Osi/.libs/*.o
OSIGRBLIBS  = $(OSISRC)/Osi/src/OsiGrb/.libs/*.o
OSICOMPILED = $(OSISRC)/Osi/src/OsiGrb/.libs/libOsiGrb.so
COINLIBS    = $(OSISRC)/CoinUtils/src/.libs/*.o

DYNLINKS = -llapack -lz -lbz2

GRBMAJVER = 46
GRBMINVER = 1
GRBHOST = linux64
GRBDIR = $(MAINDIR)/solvers/osigrb/gurobi$(GRBMAJVER)$(GRBMINVER)/$(GRBHOST)/
GRBINCS = -I$(GRBDIR)/include/
GRBLIBS = -L$(GRBDIR)/lib/ -lgurobi$(GRBMAJVER) -lm -pthread

DEPSCOMPILED = $(OSIDIR)/obj/Osi.o $(MIPDIR)/obj/MipWrapper.o
INCLUDES     = $(OSIFLAGS) $(MIPFLAGS) $(OSIINCS) $(OSIGRBINCS) $(COININCS) $(GRBINCS)
LIBRARIES    = $(OSILIBS) $(OSIGRBLIBS) $(COINLIBS) $(DEPSCOMPILED)
#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE	= $(PYTHON)/$(MAINNAME).i
PYFILE	= $(PYTHON)/$(MAINNAME).py
WXXFILE	= $(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE	= $(PYTHON)/$(MAINNAME)_wrap.o
SOFILE	= $(PYTHON)/_$(MAINNAME).so

CFLAGS += -Wall -ffloat-store -lm $(DYNLINKS) $(INCLUDES) -O3 #-g
LFLAGS = $(LIBRARIES) $(GRBLIBS)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++

wrapper: $(IFILE) $(WXXFILE) $(DEPSCOMPILED) $(OBJDIR)/OsiGrb.o $(OSICOMPILED) $(SOFILE)

$(OSICOMPILED):
	cd $(OSISRC) && (./configure --with-gurobi-incdir=$(GRBDIR)/include --with-gurobi-lib=$(GRBLIBS) && make)

$(MIPDIR)/obj/%.o: $(MIPDIR)/src/MipWrapper.cpp
	@echo "-> compiling $@"
	$(CXX) $(CFLAGS) -fPIC -c -o $@ $<

$(OSIDIR)/obj/%.o: $(OSIDIR)/src/Osi.cpp
	@echo "-> compiling $@"
	$(CXX) $(CFLAGS) -fPIC -c -o $@ $<

install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install

clean_swig:
	-rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~

clean_all: clean clean_swig
	cd $(CORE)/; make clean

clean:
	-rm -rf $(OBJDIR)/* $(PYTHON)/*o $(PYTHON)/*~ $(PYTHON)/build *~ $(SRCDIR)/*~

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(CFLAGS) -fPIC -c -o $@ $<

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
		@echo "-> compiling $@"
		$(CXX) $(CFLAGS) -I$(SRCDIR) -I$(PYLIB) -fPIC -c $< -o $@

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"
	$(SWIG) -c++ -python $<
	python $(MAINDIR)/tools/finalize.py $(<:.i=.py) MipWrapper

$(PYTHON)/%.i: $(SRCDIR)/%.hpp $(OSIDIR)/src/Osi.hpp $(MIPDIR)/src/MipWrapper.hpp
	-python $(MAINDIR)/tools/mk_subinterface.py osi Osi $< mip MipWrapper

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(DEPSCOMPILED) $(OSICOMPILED)
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(LFLAGS) $(PYTHON)/OsiGrb_wrap.o -o $@

#---- EOF --------------------------------------------------------------------
