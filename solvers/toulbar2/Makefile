include ../../make.project

MAINDIR = $(shell pwd)/../..
TOULBAR2DIR = .

OSTYPE          :=      $(shell uname -s | \
                        tr '[:upper:]' '[:lower:]' | \
                        sed \
			-e s/cygwin.*/cygwin/ \
                        -e s/irix../irix/ )

PYVERSION	:=	$(shell python -V 2>&1 | \
                        tr '[:upper:]' '[:lower:]' | \
                        sed \
			-e s/" 2.4."*/2.4/ \
			-e s/" 2.5."*/2.5/ \
			-e s/" 2.6."*/2.6/ \
			-e s/" 2.7."*/2.7/ )


ifeq ($(OSTYPE),darwin)
BUNDLE = -bundle -undefined suppress 
HACK = -I/usr/lib/gcc/i686-apple-darwin9/4.2.1/include
else
BUNDLE = -shared 
endif


ifeq ($(PYVERSION),python2.7)
PYLIB = `python2.7-config --cflags` $(HACK) 
PYOBJ = `python2.7-config --prefix`
SWIG = swig
else  
ifeq ($(PYVERSION),python2.6)
PYLIB = `python2.6-config --cflags` $(HACK) 
PYOBJ = `python2.6-config --prefix`
SWIG = swig 
else 
ifeq ($(PYVERSION),python2.5)
PYLIB = `python2.5-config --cflags` 
PYOBJ = `python2.5-config --prefix`
SWIG = swig
else 
ifeq ($(PYVERSION),python2.4)
PYLIB = `/usr/bin/python2.4-config --cflags` $(HACK)  
PYOBJ = `/usr/bin/python2.4-config --prefix`
SWIG = swig -classic
else 
PYLIB = `python-config --cflags` $(HACK) 
PYOBJ = `python-config --prefix`
SWIG = swig
endif
endif
endif
endif

OPTFLAGS = -g -O3 -DNDEBUG
# -O3 -DNDEBUG
# -D_DEBUGWRAP

COMPILFLAGS = -Wall \
-D_UNIX \
-fPIC \
-DLINUX -DLONGLONG_COST -DWIDE_STRING -DLONGDOUBLE_PROB -DNARYCHAR \

CCC = g++ $(OPTFLAGS) $(COMPILFLAGS)

LIBSRC=$(TOULBAR2DIR)/lib/src
LIBOBJ=$(TOULBAR2DIR)/lib/obj
LIBINC=$(TOULBAR2DIR)/include

OBJDIR = $(TOULBAR2DIR)/obj
SRCDIR = $(TOULBAR2DIR)/src
PYTHON = $(TOULBAR2DIR)/python


CFLAGS += -I$(LIBINC) -I$(SRCDIR)\
`xml2-config --cflags`

LFLAGS = -L$(LIBOBJ) \
`xml2-config --libs` \
-lgmp \


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	Toulbar2
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


PINCSRC = $(wildcard $(LIBINC)/*.hpp)
PLIBSRC = $(wildcard $(LIBSRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(LIBSRC)/%, $(LIBOBJ)/%, $(PLIBAUX))


wrapper: $(OBJDIR) $(PYTHON) $(LIBOBJ) $(IFILE) $(WXXFILE) $(SOFILE) 

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(PYTHON):
	@mkdir -p $(PYTHON)

$(LIBOBJ):
	@mkdir -p $(LIBOBJ)

# ------------------------------------------------------------
clean_swig: 
	rm -rf $(OBJDIR)/* $(PYTHON)/* *~ $(SRCDIR)/*~

clean : 
	rm -rf $(LIBOBJ)/*.o $(LIBOBJ)/*.a $(OBJDIR)/* $(LIBSRC)/*~ $(PYTHON)/*.so $(PYTHON)/*.cxx $(PYTHON)/*.py $(LIBSRC)/*~ $(LIBINC)/*~ *~

install_python: wrapper
	cd $(PYTHON); python $(MAINDIR)/tools/setup.py install	

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CCC) $(CFLAGS) -c -o $@ $< 

$(PYTHON)/%_wrap.o: $(PYTHON)/%_wrap.cxx
		@echo "-> compiling $@"
		$(CCC) $(CFLAGS) -I$(SRCDIR) $(PYLIB)  -c $< -o $@

$(PYTHON)/%_wrap.cxx: $(PYTHON)/%.i
	@echo "-> compiling $@"	
	$(SWIG) -c++ -python $<
	python $(MAINDIR)/tools/finalize.py $(<:.i=.py)

$(PYTHON)/%.i: $(SRCDIR)/%.hpp
	python $(MAINDIR)/tools/mk_interface.py $<

$(PYTHON)/_%.so: $(OBJDIR)/%.o $(PYTHON)/%_wrap.o $(PLIBOBJ)
	@echo "-> linking $@"
	$(CCC) $(BUNDLE) $(CFLAGS) $(AUXOBJFILES) $(MAINOBJFILES) $(PLIBOBJ) python/Toulbar2_wrap.o -o $@ $(LFLAGS)

$(LIBOBJ)/%.o: $(LIBSRC)/%.cpp
	$(CCC) $(CFLAGS) -c   $< -o $@ 

.PRECIOUS: %.c %.o
.SECONDARY:
