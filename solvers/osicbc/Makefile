
include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

MINDIR         =       .

OBJDIR = ./obj
SRCDIR = ./src
PYTHON = ./python
SWIG = swig
PYLIB = `python-config --cflags` 

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	OsiCbc
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

MIPDIR = $(MAINDIR)/solvers/mip
MIPFLAGS = -I$(MIPDIR)/src/

OSIDIR = $(MAINDIR)/solvers/osi
OSISRC = $(OSIDIR)/src
OSIFLAGS = -I$(OSISRC)

#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE	= $(PYTHON)/$(MAINNAME).i
PYFILE	= $(PYTHON)/$(MAINNAME).py
WXXFILE	= $(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE	= $(PYTHON)/$(MAINNAME)_wrap.o
SOFILE	= $(PYTHON)/_$(MAINNAME).so

COINLIBS = -lCoinUtils -llapack -lOsi -lz -lbz2
CLPLIBS = -lClp -lOsiClp
CBCLIBS = -lCbc -lOsiCbc -lCgl $(CLPLIBS)
CFLAGS = -g -Wall -ffloat-store -O3 -lm $(COINLIBS) $(MIPFLAGS) $(CBCLIBS) $(OSIFLAGS)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ 

wrapper: $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/OsiCbc.o 

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

$(PYTHON)/%.i: $(SRCDIR)/%.hpp ../osi/src/Osi.hpp ../mip/src/MipWrapper.hpp
	-python $(MAINDIR)/tools/mk_subinterface.py osi Osi $< mip MipWrapper

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o ../osi/obj/Osi.o ../mip/obj/MipWrapper.o
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(MAINOBJFILES) ../mip/obj/MipWrapper.o ../osi/obj/Osi.o $(PYTHON)/OsiCbc_wrap.o -o $@

#---- EOF --------------------------------------------------------------------
