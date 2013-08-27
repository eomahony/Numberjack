
MAINDIR ?= .

COPTIMIZE ?= -O3 --param inline-unit-growth=60
#OPTFLAGS ?= -O3 #-m
#OPTFLAGS = -O3 -g
#OPTFLAGS = -g 

#COMPILFLAGS = -Wall -D_UNIX -D_BIT32 -D_DEBUG_SEARCH #-D_DEBUG_AC -D_DEBUG_PROPAG #-D_DELTA 
COMPILFLAGS ?= -D_UNIX -D_BIT32 -DNDEBUG #-D_DEBUG_SEARCH -D_DEBUG_NOGOOD -D_DEBUG_UNITPROP -D_DEBUG_WATCH #-D_DEBUG_PROPAG #-D_DEBUG_REWRITE #-D_DEBUG_AC  #-D_CHRONOLOGICAL #-D_DEBUG_AC 

CCC = g++ $(COPTIMIZE) $(COMPILFLAGS)

BIN=$(MAINDIR)/bin
SRC=$(MAINDIR)/src/lib
MOD=$(MAINDIR)/examples
OBJ=$(MAINDIR)/src/obj
INC=$(MAINDIR)/src/include
DOC=$(MAINDIR)/doc
TCL=$(MAINDIR)/tools/tclap/include

CFLAGS = -I$(INC) -I$(TCL) #-Wall -ffloat-store 
LFLAGS = -L$(OBJ)


MODELS = $(wildcard $(MOD)/src/*.cpp)
BINS = $(patsubst $(MOD)/src/%, $(BIN)/%, $(MODELS:.cpp=))


PINCSRC = $(wildcard $(INC)/*.hpp)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))



#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------


default: flatzinc

flatzinc: fz/mistral-fzn
	cp fz/mistral-fzn ./bin/fzn-mistral

fz/mistral-fz: 
	cd fz; make

all: lib $(BINS) flatzinc

# The library
lib: $(PLIBOBJ) $(PUTIOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	@$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS)   $(PLIBOBJ) $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<	
	@$(CCC) $(CFLAGS)   $(PLIBOBJ) $< -lm -o $(BIN)/$@ 

