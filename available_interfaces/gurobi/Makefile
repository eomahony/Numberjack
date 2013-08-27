include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

# GUROBI_HOME environment variable should be set according to install instructions:
# http://www.gurobi.com/documentation/current/quick-start-guide/

OBJDIR = ./obj
SRCDIR = ./src
PYTHON = ./python
SWIG = swig
PYLIB = `python-config --cflags`
PYOBJ = `python-config --prefix`


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	Gurobi
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

GRBLIBFILES = $(shell find $(GUROBI_HOME)/lib -name "libgurobi??.so")
GUROBILIB=$(strip $(patsubst $(GUROBI_HOME)/lib/lib%.so, %, $(GRBLIBFILES)))


#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE		=	$(PYTHON)/$(MAINNAME).i
PYFILE		=	$(PYTHON)/$(MAINNAME).py
WXXFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.o
SOFILE		= 	$(PYTHON)/_$(MAINNAME).so

MIPDIR = $(MAINDIR)/solvers/mip
CFLAGS += -I$(MIPDIR)/src/ -I$(GUROBI_HOME)/include
LFLAGS = -L$(GUROBI_HOME)/lib -lgurobi_c++ -l$(GUROBILIB)


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++


wrapper: check-env $(OBJDIR) $(PYTHON) $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/Gurobi.o $(PYTHON)/Gurobi_wrap.o

debug: CFLAGS += -D_DEBUGWRAP -g
debug: wrapper

check-env:
	@if [ "${GUROBI_HOME}" == "" ]; then \
		echo "The environment variable GUROBI_HOME is not set. Please set it according to:"; \
		echo "http://www.gurobi.com/documentation/current/quick-start-guide/"; \
		exit 1; \
	fi

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
	@$(CXX) $(CFLAGS) $(PYLIB) -fPIC -c $< -o $@

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
