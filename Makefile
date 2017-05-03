
# Specifying some environmental variables for Linux, note this can be done in the shell prompt

COMP	= GCC
TECIO	= NO
CODE	= RELEASE
OPENMP	= NO

# Specifing Standard Variables:
CXX	= g++ -std=gnu++11 #-pedantic-errors # c++ gcc compiler
CXXFLAGS=       # C++ compiler flags
LDLFLAGS=	# linker flags
CPPFLAGS=	# c/c++ preprocessor flags

OPTS	= 	# optimization flags and other options

# Includes

OPTS	+= -I include
#OPTS    += -I /home/mhawwary/Libraries/Eigen3.3.2/Eigen
#OPTS += -I /home/mhawwary/work/hpMusic/contrib/eigen/Eigen

ifeq ($(TECIO),YES)
	OPTS += -I $(TECIO_DIR)/tecsrc
endif

ifeq ($(CODE),RELEASE)
	ifeq ($(COMP),GCC)
		OPTS	+= -O3 
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif

ifeq ($(OPENMP),YES)	
	OPTS	+=  -lgomp -fopenmp 
endif

ifeq ($(CODE),DEBUG)
	ifeq ($(COMP),GCC)
		OPTS	+= -fbuiltin -g -Wall #-Werror
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif


# Source

SRC	= src/
OBJ	= obj/
BIN	= bin/
INC	= include/      

vpath %.cpp src
vpath %.c src
vpath %.o   obj
vpath %.h include src
vpath %.hpp include src

# Objects
OBJS	= $(OBJ)FVMFlow.o  $(OBJ)SimCase.o $(OBJ)general_tools.o $(OBJ)Mesh.o $(OBJ)SimData.o $(OBJ)Euler2DSolver.o $(OBJ)NS2DSolver.o $(OBJ)solver_tools.o $(OBJ)ExplicitTimeSolver.o  # objects 
INCLS	= 

# Compile

.PHONY: default help clean


default: all
help:	
	@echo 'help'

all: FVMFlow.exe

FVMFlow.exe: $(OBJS)
	$(CXX) $(OPTS) -o $(BIN)$@ $+


$(OBJ)%.o : %.cpp 
	$(CXX) $(OPTS) -c -o $@ $<

$(OBJ)%.o : %.c
	$(CXX) $(OPTS) -c -o $@ $<


$(OBJ)FVMFlow.o:   FVMFlow.cpp 
$(OBJ)SimCase.o:   SimCase.hpp SimCase.cpp
$(OBJ)SimData.o:   SimData.hpp SimData.cpp
$(OBJ)general_tools.o:   general_tools.h general_tools.cpp
$(OBJ)Mesh.o:   Mesh.hpp Mesh.cpp
$(OBJ)Euler2DSolver.o: Euler2DSolver.hpp Euler2DSolver.cpp 
$(OBJ)NS2DSolver.o: NS2DSolver.hpp NS2DSolver.cpp
$(OBJ)ExplicitTimeSolver.o: ExplicitTimeSolver.hpp ExplicitTimeSolver.cpp
$(OBJ)solver_tools.o: solver_tools.h solver_tools.c

clean:
	rm -f ./$(OBJ)*.o ./$(BIN)*.exe 
	#rm -rf ./post_process/ 
	@echo  removing all object and executable files
clean2:
	rm -rf ./post_process/
 

