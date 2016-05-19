#
# comments...
#

# c++ compiler
#CXX=		g++
#CXX=		icpc

# mpi...
CXX=		mpicxx
CC =		mpicc

CXXFLAGS=	-c -DHAVE_CONFIG_H -DHYPRE_TIMING 
SOURCES=	Transport.cpp Friendly_Vector.cpp  TableReading.cpp  Hypre_Poisson_Solver.cpp
OBJECTS=	$(SOURCES:.cpp=.o)
LDFLAGS=	
EXE=	final
HYPRE_DIR = /home/shima86/hypre-2.9.0b/src/hypre
HYPRE_INC = -I$(HYPRE_DIR)/include
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE -limf -lm

COPTS     = -g -Wall
LINKOPTS  = $(COPTS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -limf -lm
LFLAGS    = $(LINKOPTS) -lstdc++

all : $(SOURCES) $(EXE)  

$(EXE) : $(OBJECTS)
	$(CXX) -o  $(EXE) $(OBJECTS) $(HYPRE_LIB)

.cpp.o:
	$(CXX)  $(CXXFLAGS) $< -o $@ $(HYPRE_INC)

clean : 
	rm -rf *.o *.dat *.txt  $(EXE) *.out.* 
