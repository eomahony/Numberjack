
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

MAINNAME	=	Osi
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

MIPDIR = $(MAINDIR)/solvers/mip
MIPFLAGS = -I$(MIPDIR)/src/

#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE	= $(PYTHON)/$(MAINNAME).i
PYFILE	= $(PYTHON)/$(MAINNAME).py
WXXFILE	= $(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE	= $(PYTHON)/$(MAINNAME)_wrap.o
SOFILE	= $(PYTHON)/_$(MAINNAME).so

COINDEPS = -llapack -lz -lbz2
OSIVER = 0.105.2
OSIDIR = Osi-$(OSIVER)
OSIINCS = -I$(OSIDIR)/Osi/src/Osi
OSILIBS = $(OSIDIR)/Osi/src/Osi/.libs/*.o
COINDIR = $(OSIDIR)/CoinUtils
COININCS = -I$(COINDIR)/src
COINLIBS = $(COINDIR)/src/.libs/*.o
CFLAGS = -g -Wall -ffloat-store -lm $(COININCS) $(OSIINCS) $(MIPFLAGS) $(COINDEPS)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ 

wrapper: $(IFILE) $(WXXFILE) $(SOFILE) $(OBJDIR)/Osi.o $(OSIDIR)/Osi/src/Osi/.libs/libOsi.so
$(OSIDIR)/Osi/src/Osi/.libs/libOsi.so:
	cd $(OSIDIR) && make || (./configure && make)

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

$(PYTHON)/%.i: $(SRCDIR)/%.hpp $(MIPDIR)/src/MipWrapper.hpp
	-python $(MAINDIR)/tools/mk_subinterface.py mip MipWrapper $<

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(MIPDIR)/obj/MipWrapper.o
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(LFLAGS) $(MAINOBJFILES) $(MIPDIR)/obj/MipWrapper.o $(OSILIBS) $(COINLIBS) $(PYTHON)/Osi_wrap.o -o $@

#---- EOF --------------------------------------------------------------------
