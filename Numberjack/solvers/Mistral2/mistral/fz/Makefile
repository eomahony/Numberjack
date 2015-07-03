
MAINDIR   = ..
CHDRS     = $(wildcard *.hpp) $(wildcard $(MTL)/*.h)
EXEC      = mistral-fzn


main: $(EXEC)



COPTIMIZE=-O3
COMPILFLAGS=-D_UNIX -D_BIT32 #-DNDEBUG #-D_DEBUG_SEARCH -D_DEBUG_NOGOOD #-D_DEBUG_AC #-D_DEBUG_REWRITE #-D_DEBUG_SEARCH #-D_DEBUG_AC
#-D_FLATZINC_OUTPUT #-D_DEBUG_SEARCH #-D_MONITOR #-D_DEBUG_AC  #-D_DEBUG_PROPAG # # -D_DEBUG_REWRITE -D_DEBUG_MUL #-D_DEBUG_PROPAG #-D_DEBUG_NOGOOD -D_DEBUG_UNITPROP -D_DEBUG_WATCH #-D_DEBUG_PROPAG #-D_DEBUG_REWRITE #-D_DEBUG_AC  #-D_CHRONOLOGICAL #-D_DEBUG_AC 


include $(MAINDIR)/template.mk

COREOBJS  = $(PLIBOBJ)
CSRCS     = $(wildcard *.cpp) lexer.yy.cpp parser.tab.cpp
COBJS     = $(addsuffix .o, $(basename $(CSRCS))) $(COREOBJS) 



include ./template.mk




lexer.yy.cpp: lexer.lxx
	flex -o lexer.yy.cpp lexer.lxx

parser.tab.cpp parser.tab.hh: parser.yxx
	bison -t -o parser.tab.cpp -d $<

