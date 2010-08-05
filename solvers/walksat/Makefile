
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
PYLIB = `python-config --cflags` $(HACK) 

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	Walksat
MAINOBJ		=	$(MAINNAME).o
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





#-----------------------------------------------------------------------------
# Walksat
#-----------------------------------------------------------------------------


SATDIR    = $(MAINDIR)/solvers/sat
MCORE     = $(MAINDIR)/solvers/sat/minisat_src/core
MTL       = $(MAINDIR)/solvers/sat/minisat_src/mtl
WALKSRC   = $(MAINDIR)/solvers/sat/walksat_src
CFLAGS    = -I$(MCORE) -I$(MTL) -I$(SATDIR)/src/ -I$(WALKSRC) $(HACK) -Wall -ffloat-store -O3 -D NDEBUG #-D_DEBUGWRAP


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ #-g -D_DEBUGWRAP


wrapper: $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/Walksat.o 

$(WALKSRC)/%.o: $(WALKSRC)/%.cpp
	@echo Compiling: "$@ ( $< )"
	gcc $(CFLAGS) -fPIC -c -o $@ $<

$(SATDIR)/obj/%.o: $(SATDIR)/src/%.cpp
	@echo Compiling: "$@ ( $< )"
	gcc $(CFLAGS) -fPIC -c -o $@ $<

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
	python $(MAINDIR)/tools/finalize.py $(<:.i=.py) SatWrapper

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	python $(MAINDIR)/tools/mk_subinterface.py sat SatWrapper $<


$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(WALKSRC)/cpp_walksat.o $(SATDIR)/obj/SatWrapper.o
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(PYTHON)/Walksat_wrap.o $(WALKSRC)/cpp_walksat.o $(SATDIR)/obj/SatWrapper.o -o $@

#---- EOF --------------------------------------------------------------------
