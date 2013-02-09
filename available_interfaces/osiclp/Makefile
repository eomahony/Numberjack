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
MAINNAME	=	OsiClp
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

MIPDIR   = $(MAINDIR)/solvers/mip
MIPFLAGS = -I$(MIPDIR)/src/
MIPLIB   = $(MIPDIR)/obj/MipWrapper.o

OSIDIR = $(MAINDIR)/solvers/osi
OSIFLAGS = -I$(OSIDIR)/src/
OSILIB = $(OSIDIR)/obj/Osi.o

INCLUDES  = $(OSIFLAGS) $(MIPFLAGS) `pkg-config osi-clp --cflags`
LIBRARIES = $(OSILIB) $(MIPLIB) 

#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------
IFILE	= $(PYTHON)/$(MAINNAME).i
PYFILE	= $(PYTHON)/$(MAINNAME).py
WXXFILE	= $(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE	= $(PYTHON)/$(MAINNAME)_wrap.o
SOFILE	= $(PYTHON)/_$(MAINNAME).so

CFLAGS += -O3 -Wall -ffloat-store -lm $(INCLUDES) `pkg-config osi-clp --libs` #-g

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------
CXX = g++ 

wrapper: $(IFILE) $(WXXFILE) $(LIBRARIES) $(OBJDIR)/OsiClp.o $(SOFILE) 

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

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(LIBRARIES)
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(LIBRARIES) $(PYTHON)/OsiClp_wrap.o -o $@

#---- EOF --------------------------------------------------------------------
