include ../../make.project
MAINDIR = $(CURDIR)/../..

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

ZIBDIR = $(ZIBPATH)

SCIPDIR = $(shell find $(ZIBPATH) -name "scip-?.?.?")
SOPLEXDIR = $(shell find $(ZIBPATH) -name "soplex-?.?.?")

PYTHON = ./python

ifeq ($(OSTYPE),darwin)
	ifeq ($(ARCH),x86_64)
		# Fixes an issue with gcc4.2 creating two 64bit binarys and trying to link them with `lipo`, even though it is given -arch i386 -arch x86_64 by python-config.
		export ARCHFLAGS = -arch x86_64
	endif
endif

ifeq ($(PYVERSION),python2.6)
PYLIB = `python2.6-config --cflags`
PYOBJ = `python2.6-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.5)
PYLIB = `python2.5-config --cflags` 
PYOBJ = `python2.5-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.4)
PYLIB = `/usr/bin/python2.4-config --cflags` 
PYOBJ = `/usr/bin/python2.4-config --prefix`
SWIG = swig -classic
else 
PYLIB = `python-config --cflags`
PYOBJ = `python-config --prefix`
SWIG = swig
endif
endif
endif

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
include ../../make.project  # We need to include this again because SCIP's make.project clears CFLAGS and we need to re-add -flat_namespace for compiling on Mac.
ZIMPL=false 
ZLIB=false 
READLINE=false

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	SCIP
MAINOBJ		=	SCIP.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))

MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

OBJDIR		=   ./obj

#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE		=	$(PYTHON)/$(MAINNAME).i
PYFILE		=	$(PYTHON)/$(MAINNAME).py
WXXFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.cxx
WOBFILE		= 	$(PYTHON)/$(MAINNAME)_wrap.o
SOFILE		= 	$(PYTHON)/_$(MAINNAME).so

MIPDIR          =       $(MAINDIR)/solvers/mip
MIPFLAGS        =       -I$(MIPDIR)/src/ 
CXXFLAGS += $(CFLAGS)


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

NFLAGS = $(FLAGS)

wrapper: check-env $(OBJDIR) $(PYTHON) $(SCIPDIR)/done.txt $(IFILE) $(WXXFILE) $(SOFILE)

debug: NFLAGS += -D_DEBUGWRAP -g
debug: wrapper

check-env:
	@if test "$(ZIBPATH)" = ""; then \
		echo ZIBPATH is undefined. Please set the path to scipoptsuite.; \
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
	@-rm -rf $(OBJDIR)/* $(PYTHON)/*o $(PYTHON)/*~ $(PYTHON)/build *~ $(SRCDIR)/*~ $(SCIPDIR)/done.txt

clean_swig:
	@echo "-> cleaning swig interface"
	@-rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~ 

clean_all: clean clean_swig
	cd $(SCIPDIR); make clean


$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
	@echo "-> compiling $@"
	@$(CXX) $(NFLAGS) $(MIPFLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -fPIC -c $< $(CXX_o)$@

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
	@echo "-> compiling $@"
	@$(CXX) $(PYLIB) $(MIPFLAGS) $(NFLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -fPIC -c $< $(CXX_o)$@

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"	
	@$(SWIG) -c++ -python $<
	@python $(MAINDIR)/tools/finalize.py $(<:.i=.py) MipWrapper

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	@echo "-> creating interface $@"
	@python $(MAINDIR)/tools/mk_subinterface.py mip MipWrapper $<


$(PYTHON)/_%.so: $(SCIPLIBFILE) $(NLPILIBFILE) $(LPILIBFILE) $(MAINOBJFILES) $(OBJDIR)/%.o  $(MIPDIR)/obj/MipWrapper.o $(PYTHON)/%_wrap.o
	@echo "-> linking $@"
	@$(LINKCXX) $(BUNDLE) -Wno-long-double $(CFLAGS) $(MAINOBJFILES) $(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX)  $(OFLAGS) $(LPSLDFLAGS) $(MIPDIR)/obj/MipWrapper.o python/SCIP_wrap.o $(LINKCXX_o)$@

#scip_lib: $(SCIPDIR)/Makefile.backup 
#	cp  $(SCIPDIR)/Makefile $(SCIPDIR)/Makefile.backup
#	cp  $(MAINDIR)/tools/Makefile.scip $(SCIPDIR)/Makefile
#	cp  $(SOPLEXDIR)/Makefile $(SOPLEXDIR)/Makefile.backup
#	cp  $(MAINDIR)/tools/Makefile.soplex $(SOPLEXDIR)/Makefile
#	cd  $(SCIPDIR)/..; make ZIMPL=false ZLIB=false READLINE=false

$(SCIPDIR)/done.txt: 
	@cp $(MAINDIR)/tools/Makefile.Numberjack $(SCIPDIR)/../
	@cp $(MAINDIR)/tools/Makefile.doit.Numberjack $(SCIPDIR)/../
	@cp $(MAINDIR)/tools/Makefile.scip.Numberjack $(SCIPDIR)/
	@cp $(MAINDIR)/tools/Makefile.soplex.Numberjack $(SOPLEXDIR)/
	@cd $(SCIPDIR)/..; make -f Makefile.Numberjack ZIMPL=false ZLIB=false READLINE=false	
	@echo finished | cat > $(SCIPDIR)/done.txt

#---- EOF --------------------------------------------------------------------
