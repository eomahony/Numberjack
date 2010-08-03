
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

MAINNAME	=	SatWrapper
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





MCORE     = minisat_src/core
MTL       = minisat_src/mtl
CFLAGS    = -I$(MCORE) -I$(MTL) $(HACK) -Wall -ffloat-store -O3 -D NDEBUG #-D_DEBUGWRAP

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ 


wrapper: $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/SatWrapper.o 

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


$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o 
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(LFLAGS) $(MAINOBJFILES) $(PYTHON)/SatWrapper_wrap.o -o $@

#---- EOF --------------------------------------------------------------------
