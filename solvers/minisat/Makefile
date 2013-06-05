
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


SATDIR    = $(MAINDIR)/solvers/sat
MTL       = $(SATDIR)/minisat_src/mtl
CORE      = $(SATDIR)/minisat_src/core
CFLAGS    += -I$(SATDIR)/src/ -I$(MTL) $(HACK) -I$(CORE) -Wall -ffloat-store -O3 # -DNDEBUG
#LFLAGS    = -lz


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ #-g -D_DEBUGWRAP


wrapper: $(OBJDIR) $(PYTHON) $(SATDIR)/$(OBJDIR) $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/MiniSat.o

debug: CFLAGS += -D_DEBUGWRAP -g
debug: wrapper

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(PYTHON):
	@mkdir -p $(PYTHON)

$(SATDIR)/$(OBJDIR):
	@mkdir -p $(SATDIR)/$(OBJDIR)

$(CORE)/%.o: $(CORE)/%.C 
	@echo "-> compiling: $@ ( $< )"
	@$(CXX) $(CFLAGS) -fPIC -c -o $@ $<

$(SATDIR)/obj/%.o: $(SATDIR)/src/%.cpp
	@echo "-> compiling: $@ ( $< )"
	@gcc $(CFLAGS) -fPIC -c -o $@ $<

install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install	

clean_swig:
	@echo "-> cleaning swig interface"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~

clean_all: clean clean_swig
	@cd $(CORE)/; make clean


clean:
	@echo "-> cleaning"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/*o $(PYTHON)/*~ $(PYTHON)/build *~ $(SRCDIR)/*~


$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
	@echo "-> compiling $@"
	@$(CXX) $(CFLAGS) -fPIC -c -o $@ $< 

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
	@echo "-> compiling $@"
	@$(CXX) $(CFLAGS) -I$(SRCDIR) -I$(PYLIB) -fPIC -c $< -o $@

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"	
	@$(SWIG) -c++ -python  $<
	@python $(MAINDIR)/tools/finalize.py $(<:.i=.py) SatWrapper


$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	@echo "-> creating interface $@"
	@python $(MAINDIR)/tools/mk_subinterface.py sat SatWrapper $<


$(PYTHON)/_%.so: $(AUXOBJFILES) $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(CORE)/Solver.o $(SATDIR)/obj/SatWrapper.o
	@echo "-> linking $@"
	@$(CXX) $(BUNDLE) -Wno-long-double $(CFLAGS) $(AUXOBJFILES) $(MAINOBJFILES) $(CORE)/Solver.o $(SATDIR)/obj/SatWrapper.o python/MiniSat_wrap.o -o $@



#---- EOF --------------------------------------------------------------------
