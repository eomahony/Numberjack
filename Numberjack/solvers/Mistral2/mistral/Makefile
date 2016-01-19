


COMPILFLAGS=-D_UNIX -D_BIT32 -D_DEBUG_AC -D_DEBUG_PROPAGATE -D_DEBUG_SEARCH #-DNDEBUG 
#COPTIMIZE = -O3
#COMPILFLAGS = -Wint-to-pointer-cast -D_UNIX -D_BIT32 -DNDEBUG -D_DEBUG_SEARCH -D_DEBUG_AC #-D_DEBUG_NOGOOD -D_DEBUG_FAIL # -D_DEBUG_AC # #-D_DEBUG_GENPROPAG #-D_MONITOR -D_DEBUG_SEARCH  #-D_DEBUG_NOGOOD # #-D_DEBUG_PROPAGATE #-D_DEBUG_CARRAY  #-D_VARNCONQUEUE  -D_DEBUG_PRUNING #-D_DEBUG_AC -D_DEBUG_QUEUE #-D_VARNCONQUEUE #-D_PROFILING #-D_DEBUG_SEARCH #-D_DEBUG_PROPAG # -D_DEBUG_PRUNING -D_DEBUG_AC -D_DEBUG_QUEUE #-D_DEBUG_REWRITE #-D_DEBUG_CGRAPH #-D_DEBUG_PROPAG -D_DEBUG_MUL #-D_DEBUG_AC # -D_DEBUG_UNITPROP -D_DEBUG_WATCH #-D_DEBUG_PROPAG #-D_DEBUG_REWRITE #-D_DEBUG_AC  #-D_CHRONOLOGICAL #-D_DEBUG_AC 


include ./template.mk

DATE := $(shell date '+%y-%m-%d')

clean : 
	rm -rf $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~ $(UTI)/*~  *~ $(BIN)/* $(DOC)/*~ ./fzn-mistral fz/mistral-fzn

archive: 
	@echo Export Mistral version 2.0.$(DATE)
	rm -rf Mistral-2.0.$(DATE).tar.gz
	rm -rf Mistral-2.0.$(DATE)
	mkdir Mistral-2.0.$(DATE)
	git archive master --format=tar | tar -x -C ./Mistral-2.0.$(DATE)
	mkdir ./Mistral-2.0.$(DATE)/globals
	#cp -rf ./fz/mznlib/*.mzn ./Mistral-2.0.$(DATE)/globals 
	tar -cvzf Mistral-2.0.$(DATE).tar.gz Mistral-2.0.$(DATE)

test: bin/unit_test
	cd fz; make; 
	bin/unit_test
	cd fz; ./test_fcts.sh passed
