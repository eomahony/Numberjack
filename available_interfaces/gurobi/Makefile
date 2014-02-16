include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

# Numberjack will try to find the Gurobi, install location. It will first check
# for the GUROBI_HOME environment variable which should be set accodring to the
# Gurobi install guide. If not, it will try to find the path based on the
# location of the gurobi_cl executable.

OBJDIR = ./obj
SRCDIR = ./src
PYTHON = ./python
SWIG = swig
PYLIB = `python-config --cflags`
PYOBJ = `python-config --prefix`

GRBCCFLAGS = `python $(MAINDIR)/tools/solverconfig.py Gurobi --cflags`
GRBLDFLAGS = `python $(MAINDIR)/tools/solverconfig.py Gurobi --lflags`


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	Gurobi
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

MIPDIR = $(MAINDIR)/solvers/mip
CFLAGS += -I$(MIPDIR)/src/ $(GRBCCFLAGS)
LFLAGS += $(GRBLDFLAGS)


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX ?= g++


wrapper: $(OBJDIR) $(PYTHON) $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/Gurobi.o $(PYTHON)/Gurobi_wrap.o

debug: CFLAGS += -D_DEBUGWRAP -g
debug: wrapper

install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install

$(OBJDIR):
	@-mkdir -p $(OBJDIR)

$(PYTHON):
	@mkdir -p $(PYTHON)

clean:
	@echo "-> cleaning"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/*o $(PYTHON)/*~ $(PYTHON)/build *~ $(SRCDIR)/*~

clean_swig:
	@echo "-> cleaning swig interface"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~ 

clean_all: clean clean_swig

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
	@echo "-> compiling $@"
	@$(CXX) $(CFLAGS) -fPIC -c $< -o $@

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
	@echo "-> compiling $@"
	@$(CXX) $(PYLIB) $(CFLAGS) -fPIC -c $< -o $@

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"	
	@$(SWIG) -c++ -python $<
	@python $(MAINDIR)/tools/finalize.py $(<:.i=.py) MipWrapper

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	@echo "-> creating interface $@"
	@python $(MAINDIR)/tools/mk_subinterface.py mip MipWrapper $<

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(MIPDIR)/obj/MipWrapper.o $(PYTHON)/%_wrap.o 
	@echo "-> linking $@"
	@$(CXX) $(BUNDLE) -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(LFLAGS) $(MIPDIR)/obj/MipWrapper.o $(PYTHON)/Gurobi_wrap.o -o $@
