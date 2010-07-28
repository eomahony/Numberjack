
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

MAINNAME	=	MiniSat
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

AUXNAME		=	SimpSolver
AUXOBJ		=	$(AUXNAME).o
AUXSRC		=	$(addprefix $(SRCDIR)/,$(AUXOBJ:.o=.cpp))
AUXOBJFILES	=	$(addprefix $(OBJDIR)/,$(AUXOBJ))


#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE		=	$(PYTHON)/$(MAINNAME).i
PYFILE		=	$(PYTHON)/$(MAINNAME).py
WXXFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.o
SOFILE		= 	$(PYTHON)/_$(MAINNAME).so





#-----------------------------------------------------------------------------
# MiniSat
#-----------------------------------------------------------------------------


MTL       = $(MINDIR)/minisat_src/mtl
CORE      = $(MINDIR)/minisat_src/core
CFLAGS    = -I$(MTL) $(HACK) -I$(CORE) -Wall -ffloat-store -O3 -D NDEBUG #-D_DEBUGWRAP
LFLAGS    = -lz


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ #-g -D_DEBUGWRAP


wrapper: $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/MiniSat.o 

$(CORE)/%.o: $(CORE)/%.C 
	@echo Compiling: "$@ ( $< )"
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
	python $(MAINDIR)/tools/finalize.py $(<:.i=.py)

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	python $(MAINDIR)/tools/mk_interface.py $<


$(PYTHON)/_%.so: $(AUXOBJFILES) $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(CORE)/Solver.o 
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(LFLAGS) $(AUXOBJFILES) $(MAINOBJFILES) $(CORE)/Solver.o python/MiniSat_wrap.o -o $@



#---- EOF --------------------------------------------------------------------
