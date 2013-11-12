include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

OBJDIR = ./obj
SRCDIR = ./src
PYTHON = ./python
MIPDIR = $(MAINDIR)/solvers/mip
SWIG = swig
PYLIB = `python-config --cflags`
PYOBJ = `python-config --prefix`


CPLEXCCFLAGS = `python $(MAINDIR)/tools/solverconfig.py CPLEX cflags`
CPLEXLDFLAGS = `python $(MAINDIR)/tools/solverconfig.py CPLEX lflags`


OPTFLAGS = -O3
COMPILFLAGS = -Wall -D_UNIX -fPIC -I$(MIPDIR)/src/ $(CPLEXCCFLAGS)

CCC ?= g++


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	CPLEX
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

#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------


wrapper: $(OBJDIR) $(PYTHON) $(IFILE) $(WXXFILE) $(SOFILE)

debug: CFLAGS += -D_DEBUGWRAP -g
debug: wrapper $(MAINOBJFILES) $(WOBFILE)

test:
	echo CPLEXDIR = $(CPLEXDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(PYTHON):
	@mkdir -p $(PYTHON)

clean:
	@echo "-> cleaning"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/*o $(PYTHON)/*~ $(PYTHON)/build *~ $(SRCDIR)/*~

clean_swig:
	@echo "-> cleaning swig interface"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~ 

install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install	

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
	@echo "-> compiling $@"
	@$(CCC) $(OPTFLAGS) $(COMPILFLAGS) $(CFLAGS) $(CCFLAGS) -c -o $@ $< 

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
	@echo "-> compiling $@"
	@$(CCC) $(OPTFLAGS) $(COMPILFLAGS) $(CFLAGS) $(PYLIB) -I$(SRCDIR) -c $< -o $@

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"	
	@$(SWIG) -c++ -python $<
	@python $(MAINDIR)/tools/finalize.py $(<:.i=.py) MipWrapper

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	@echo "-> creating interface $@"
	@python $(MAINDIR)/tools/mk_subinterface.py mip MipWrapper $<

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(MIPDIR)/obj/MipWrapper.o $(PYTHON)/%_wrap.o
	@echo "-> linking $@"
	@$(CCC) $(OPTFLAGS) $(COMPILFLAGS) $(BUNDLE) -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(LFLAGS) $(MIPDIR)/obj/MipWrapper.o python/$(MAINNAME)_wrap.o -o $@ $(CPLEXLDFLAGS)
