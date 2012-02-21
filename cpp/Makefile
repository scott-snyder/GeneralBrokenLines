# #################################################################
# Makefile for GBLpp (test general broken lines)
# #################################################################
#
# ROOT stuff
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --libs)  -lMathCore -lThread
#
CCOMP = g++
C_FLAGS = -Wall -O3 $(ROOTCFLAGS)
# profiling: -pg
C_LIBS = $(ROOTLIBS)
DEBUG =          # e.g. -g
#
LOADER = g++
L_FLAGS = -Wall -O3 $(ROOTLDFLAGS)
#
# objects for this project
#
USER_OBJ_GBL = GBLpp.o BorderedBandMatrix.o GblData.o GblPoint.o GblTrajectory.o MilleBinary.o VMatrix.o example1.o
#
EXECUTABLES = GBLpp
#

all:	$(EXECUTABLES)

# The single special one:
GBLpp: ${USER_OBJ_GBL} Makefile
	$(LOADER) $(L_FLAGS) $(C_LIBS) \
		-o $@ ${USER_OBJ_GBL} 
#  
clean:
	rm -f *.o *~ */*.o */*~
#
clobber: clean 
	rm -f $(EXECUTABLES)

install: $(EXECUTABLES) #clean
	mkdir -p bin
	mv $(EXECUTABLES) bin

# Make the object files - depend on source and include file 
#
%.o : %.cpp %.h Makefile 
	$(CCOMP) -c $(C_FLAGS) $(DEFINES) $(C_INCLUDEDIRS) $(DEBUG) -o $@ $<
#
# ##################################################################
# END
# ##################################################################
