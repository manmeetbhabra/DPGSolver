# Makefile
# See the GNU Make manual for guidelines.

# Instructions:
# 	Place "Makefile" in the top directory
# 	Place source code "*.c" and "*.h" files in "src" folder
#
# Additional Options:
# 	make clean


# C standard and compiler
CSTD = -std=c99
# CC   = mpicc -fopenmp -m64 

# Options
# OPTS = -O3
OPTS = -g -Wall -Wextra

# Standard libraries (Math)
# STD_LIB = -lm


# Includes
KERNEL  := $(shell uname -s)
MACHINE := $(shell uname -m)

#all:
#	@echo $(KERNEL)
#	@echo ${MACHINE}


# OSX
ifeq ($(KERNEL),Darwin)
  PROG_PATH = /Users/philipzwanenburg/Desktop/Research_Codes/Downloaded

  STD_LIB = -lm

  CC = ${PROG_PATH}/petsc/petsc-3.6.3/arch-darwin-mpich-c-debug/bin/
#  CC = ${PROG_PATH}/petsc/petsc-3.6.3/arch-darwin-mpich-c-opt/bin/
#  CC = mpicc -fopenmp -m64

#  PETSC_DIR = ${PROG_PATH}/petsc/petsc-3.2-p7 (ToBeDeleted)
  PETSC_DIR = ${PROG_PATH}/petsc/petsc-3.6.3
  PETSC_ARCH = arch-darwin-mpich-c-debug
#  PETSC_ARCH = arch-darwin-mpich-c-opt
#  PETSC_ARCH = arch-darwin-c-opt

  # METIS_DIR = ${PROG_PATH}/parmetis/parmetis-4.0.3
  METIS_DIR = ${PROG_PATH}/parmetis_mpich/parmetis-4.0.3/build/debug
#  METIS_DIR = ${PROG_PATH}/parmetis_mpich/parmetis-4.0.3/build/opt

  METIS_INC      = -I${METIS_DIR}/metis/include
  METIS_LDINC    = -L${METIS_DIR}/libmetis -lmetis
  PARMETIS_INC   = -I${METIS_DIR}/include
  PARMETIS_LDINC = -L${METIS_DIR}/libparmetis -lparmetis

  MKL_DIR   = ${PROG_PATH}/intel/mkl
  MKL_INC = -I${MKL_DIR}/include
  # MKL statically linked on OSX as the -Wl,--no-as-needed option is not supported by the OSX linker
  MKL_LDINC = ${MKL_DIR}/lib/libmkl_intel_lp64.a ${MKL_DIR}/lib/libmkl_core.a	${MKL_DIR}/lib/libmkl_sequential.a -lpthread

  OP_SYS    = Darwin-x86_64

  PETSC_INC = -I PETSC_DIR/include ${PETSC_CC_INCLUDES}
  # Run PETSC's 'variables' makefile
	include ${PETSC_DIR}/lib/petsc/conf/variables
endif
													
# LINUX
ifeq ($(KERNEL),Linux)
  CC   = mpicc -fopenmp -m64 
  STD_LIB =

  PETSC_DIR = /software/CentOS-6/libraries/petsc-3.5.3-openmpi-1.6.3-intel
  # Specify PROG_PATH in .bashrc
  METIS_DIR = ${PROG_PATH}/parmetis-4.0.2
  MKL_INC   = -I$(MKLROOT)/include
  MKL_LDINC = -Wl,--no-as-needed -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core	-lmkl_gnu_thread -ldl -lpthread -lgomp

  METIS_INC      = -I${METIS_DIR}/metis/include
  METIS_LDINC    = -L${METIS_DIR}/build/${OP_SYS}/libmetis -lmetis
  PARMETIS_INC   = -I${METIS_DIR}/include
  PARMETIS_LDINC = -L${METIS_DIR}/build/${OP_SYS}/libparmetis -lparmetis

  OP_SYS    = Linux-x86_64

	PETSC_INC = -I PETSC_DIR/include ${PETSC_CC_INCLUDES}
  # Run PETSC's 'variables' makefile
	include ${PETSC_DIR}/conf/variables
endif

#PETSC_INC = -I PETSC_DIR/include ${PETSC_CC_INCLUDES}
# Run PETSC's 'variables' makefile
#ifneq (,$(findstring 3.2-p7,$(PETSC_DIR)))
#  include ${PETSC_DIR}/conf/variables
#else
#	include ${PETSC_DIR}/lib/petsc/conf/variables
#endif


# missing LIBPATH, DEFINES
LIBS = $(STD_LIB) $(PETSC_INC) $(PETSC_LIB) $(METIS_INC) $(METIS_LDINC) $(PARMETIS_INC) $(PARMETIS_LDINC) $(MKL_INC) $(MKL_LDINC)

EXECUTABLE = DPGSolver.exe

SRCDIR  = src
OBJDIR  = obj
EXECDIR = bin

# OBJECTS    := $(addprefix $(OBJDIR)/,$(OBJECTS))
EXECUTABLE := $(addprefix $(EXECDIR)/,$(EXECUTABLE))

SOURCES = $(wildcard $(SRCDIR)/*.c)
HEADERS = $(wildcard $(SRCDIR)/*.h)
#OBJECTS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)


# Formatting for "rules" in makefiles is as follows:
# 	target ... : prerequisites ...
#   	recipe
#   	...
#
# 	Note that a "tab" character must be placed before each line of the recipe
#
# Automatic variables:
#
# 	$@: target
# 	$^: prerequisites
# 	$<: first prerequisite only (desired if excluding headers for example) 
#
# Compile commands:
# 	-c:       generate object files 
# 	-o <arg>: output to <arg>

### Default goal + Additional required rules ###

# Compile executable file (Default goal)
$(EXECUTABLE) : $(OBJECTS)
	$(CC) -o $@ $(OPTS) $^ $(LIBS)

# Create objects
# Still need to figure out how to include header dependencies.
$(OBJECTS) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(OPTS) $(CSTD) -c -o $@ $< $(LIBS) 

# Create directories if not present
$(OBJECTS): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)

$(EXECUTABLE): | $(EXECDIR)
$(EXECDIR):
	mkdir $(EXECDIR)


### Additional Rules ###

# Make clean ($ make clean)
# ".PHONY" is included here to specify that clean is a phony target

.PHONY : clean
clean:
	rm $(EXECUTABLE) $(OBJECTS)
