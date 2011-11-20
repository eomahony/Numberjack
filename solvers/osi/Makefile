
MAINDIR = $(shell pwd)/../..

OSTYPE          :=      $(shell uname -s | \
                        tr '[:upper:]' '[:lower:]' | \
                        sed \
			-e s/cygwin.*/cygwin/ \
                        -e s/irix../irix/ )

PYVERSION	:=	$(shell python2 -V 2>&1 | \
                        tr '[:upper:]' '[:lower:]' | \
                        sed \
			-e s/" 2.4."*/2.4/ \
			-e s/" 2.5."*/2.5/ \
			-e s/" 2.6."*/2.6/ \
			-e s/" 2.7."*/2.7/ )

PYTHON = python2

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

OPTFLAGS = #-m32
COMPILFLAGS = -Wall -D_UNIX -fPIC -g #-D_DEBUG #-D_DEBUGSEARCH #-D_DEBUGPROPAG 

CCC = g++ $(OPTFLAGS) $(COMPILFLAGS)

COINLIBS= -lCoinUtils -llapack -lOsi -lz -lbz2
CLPLIBS = -lClp -lOsiClp
CBCLIBS = $(CLPLIBS) -lCbc -lCgl -lOsiCbc

MIPDIR = $(MAINDIR)/solvers/mip
MIPFLAGS = -I$(MIPDIR)/src/

BIN=./bin
SRC=./lib/src
MOD=./models
OBJ=./lib/obj
INC=./include

OBJDIR = ./obj
SRCDIR = ./src
PYTHONDIR = ./python

XMLLIBS= `xml2-config --libs`

CFLAGS = -I$(INC) $(MIPFLAGS) $(COINLIBS) $(CLPLIBS) 

LFLAGS = -L$(OBJ)

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	Osi
MAINOBJ		=	$(MAINNAME).o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

#-----------------------------------------------------------------------------
# Interface
#-----------------------------------------------------------------------------

IFILE		=	$(PYTHONDIR)/$(MAINNAME).i
PYFILE		=	$(PYTHONDIR)/$(MAINNAME).py
WXXFILE		= 	$(PYTHONDIR)/$(MAINNAME)_wrap.cxx
WOBFILE		= 	$(PYTHONDIR)/$(MAINNAME)_wrap.o
SOFILE		= 	$(PYTHONDIR)/_$(MAINNAME).so

#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------

PMODELS = $(wildcard $(MOD)/src/*.cpp)
MODAUX = $(PMODELS:.cpp=)
MODELS = $(patsubst $(MOD)/src/%, %, $(MODAUX))
BINS = $(patsubst %, $(BIN)/%, $(MODELS))
EXE = $(BINS:=;)

PINCSRC = $(wildcard $(INC)/*.h)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))


wrapper: $(IFILE) $(WXXFILE) $(SOFILE) 

# ------------------------------------------------------------
clean_swig: 
	rm -rf $(OBJDIR)/* $(PYTHONDIR)/* *~ $(SRCDIR)/*~

clean: 
	rm -rf $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~ *~ $(BIN)/* 

install_python: wrapper
	cd $(PYTHONDIR); $(PYTHON) $(MAINDIR)/tools/setup.py install	

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CCC) $(CFLAGS) -c -o $@ $< 

$(PYTHONDIR)/%_wrap.o: $(PYTHONDIR)/%_wrap.cxx
		@echo "-> compiling $@"
		$(CCC) $(CFLAGS) -I$(SRCDIR) -I$(PYLIB)  -c $< -o $@

$(PYTHONDIR)/%_wrap.cxx: $(PYTHONDIR)/%.i
	@echo "-> compiling $@"	
	$(SWIG) -c++ -python $<
	$(python) $(MAINDIR)/tools/finalize.py $(<:.i=.py) MipWrapper

$(PYTHONDIR)/%.i: $(SRCDIR)/%.hpp
	$(PYTHON) $(MAINDIR)/tools/mk_subinterface.py mip MipWrapper $<

$(PYTHONDIR)/_%.so: $(OBJDIR)/%.o $(PYTHONDIR)/%_wrap.o $(PLIBOBJ)
	@echo "-> linking $@"
	$(CCC) $(BUNDLE) -flat_namespace -Wno-long-double $(CFLAGS) $(LFLAGS) $(AUXOBJFILES) $(MAINOBJFILES) $(PLIBOBJ) $(MIPDIR)/obj/MipWrapper.o  python/$(MAINNAME)_wrap.o -o $@

# The library
lib: $(PLIBOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.h
	$(CCC) $(CFLAGS) -c   $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	$(CCC) $(CFLAGS) $< $(PLIBOBJ) -lm -o $@ $(LFLAGS)

%: $(MOD)/obj/%.o $(PLIBOBJ)
	$(CCC) $(CFLAGS) $< $(PLIBOBJ) -lm -o $(BIN)/$@ $(LFLAGS)

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	$(CCC) -c $(CFLAGS)  -o $@ $<
