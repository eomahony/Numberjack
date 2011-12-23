
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

MAINNAME	=	OsiVol
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

COINLIBS = -llapack -lz -lbz2
VOLVER = 1.3.2
VOLDIR = Vol-$(VOLVER)
VOLINCS = -I$(VOLDIR)/Vol/src/ -I$(VOLDIR)/CoinUtils/src/ -I$(VOLDIR)/Osi/src/Osi/ -I$(VOLDIR)/Vol/src/OsiVol/
VOLLIBS = $(VOLDIR)/Vol/src/.libs/*.o $(VOLDIR)/CoinUtils/src/.libs/*.o $(VOLDIR)/Osi/src/Osi/.libs/*.o $(VOLDIR)/Vol/src/OsiVol/.libs/*.o
CFLAGS = -g -Wall -ffloat-store -lm $(COINLIBS) $(MIPFLAGS) $(VOLINCS) $(OSIFLAGS)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

CXX = g++ 

wrapper: $(IFILE) $(WXXFILE) $(OSIDIR)/obj/Osi.o $(MIPDIR)/obj/MipWrapper.o $(OBJDIR)/OsiVol.o $(VOLDIR)/Vol/src/.libs/libVol.so $(SOFILE)

$(VOLDIR)/Vol/src/.libs/libVol.so:
	cd $(VOLDIR) && make || (./configure && make)

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

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(OSIDIR)/obj/Osi.o $(MIPDIR)/obj/MipWrapper.o
	@echo "-> linking $@"
	$(CXX) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(MIPDIR)/obj/MipWrapper.o $(OSIDIR)/obj/Osi.o $(VOLLIBS) $(PYTHON)/OsiVol_wrap.o -o $@

#---- EOF --------------------------------------------------------------------
