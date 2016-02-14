
MAINDIR ?= .

COPTIMIZE ?= -O3

CCC = g++ 

BIN=$(MAINDIR)/bin
SRC=$(MAINDIR)/src/lib
MOD=$(MAINDIR)/examples
OBJ=$(MAINDIR)/src/obj
INC=$(MAINDIR)/src/include
DOC=$(MAINDIR)/doc
TCL=$(MAINDIR)/tools/tclap/include

CFLAGS = -D_DEBUG_AC -D_DEBUG_PROPAGATE -D_DEBUG_SEARCH -I$(INC) -I$(TCL)  #-Wall -ffloat-store 
LFLAGS = -L$(OBJ)


MODELS = $(wildcard $(MOD)/src/*.cpp)
BINS = $(patsubst $(MOD)/src/%, $(BIN)/%, $(MODELS:.cpp=))


PINCSRC = $(wildcard $(INC)/*.hpp)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))


## Compile options
%.o:			CFLAGS +=$(COPTIMIZE)  $(COMPILFLAGS) #-ggdb -D DEBUG
%.op:			CFLAGS +=$(COPTIMIZE) -pg -ggdb -D NDEBUG
%.od:			CFLAGS +=-O0 -ggdb -D DEBUG -D INVARIANTS #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
%.or:			CFLAGS +=$(COPTIMIZE) -D NDEBUG
%.oc:                   CFLAGS +=-O0 -fprofile-arcs -ftest-coverage -ggdb -D DEBUG


#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------


default: flatzinc

flatzinc:
	cd fz; make
	cp fz/mistral-fz ./bin/fzn-mistral

#fz/mistral-fz: fz/mistral-fzn
#	cp fz/mistral-fzn ./bin

#fz/mistral-fzn: 
#	cd fz; make

all: lib $(BINS) flatzinc

# The library
lib: $(PLIBOBJ) $(PUTIOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(PLIBOBJ) $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<	
	$(CCC) $(CFLAGS) $(PLIBOBJ) $< -lm -o $(BIN)/$@

