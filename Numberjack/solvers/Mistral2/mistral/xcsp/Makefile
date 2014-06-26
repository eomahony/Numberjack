
MAINDIR   = ..
EXEC      = ../bin/mistral-xcsp


main: $(EXEC)



COPTIMIZE = -O3
COMPILFLAGS = -Wall -D_UNIX -D_BIT32 -DNDEBUG #-D_DEBUG_SEARCH #-D_MONITOR #-D_DEBUG_AC  #-D_DEBUG_PROPAG # # -D_DEBUG_REWRITE -D_DEBUG_MUL #-D_DEBUG_PROPAG #-D_DEBUG_NOGOOD -D_DEBUG_UNITPROP -D_DEBUG_WATCH #-D_DEBUG_PROPAG #-D_DEBUG_REWRITE #-D_DEBUG_AC  #-D_CHRONOLOGICAL #-D_DEBUG_AC 


include $(MAINDIR)/template.mk

COREOBJS  = $(PLIBOBJ)
CSRCS     = xsolve.cpp
COBJS     = $(addsuffix .o, $(basename $(CSRCS))) $(COREOBJS) 


xsolve.o: xsolve.cpp xml/MistralCallback.hh $(PINCSRC) $(PLIBSRC)
	@echo 'compile '$<
	$(CCC) $(CFLAGS)`xml2-config --cflags` -c $< -o $@ 

$(EXEC): $(COBJS)
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(COBJS) `xml2-config --libs` -lm -o $(EXEC)

clean:
	rm -f *o *~ $(EXEC)

