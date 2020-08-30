# Please make copies of this file called `parameters_[debug/release].sh` and use parameters specific to your system.

export CMAKE_PREFIX_PATH=/home/philip/software/mpich/mpich-3.2.1/build
export MKLROOT=/home/philip/software/mkl/mkl_2018_2/mkl
export PETSC_DIR=/home/philip/software/petsc/petsc-3.9.2
export PETSC_ARCH=arch-macOS-mkl-mpich-debug
export SLEPC_DIR=/home/pzwan/software/slepc/slepc-3.9.1

# Specific paths for libraries needed for optimization. 
# This will be needed to run any optimization cases on the IGA_Manmeet branch.
# Here, specify the path to the gfortran library file (should be in the gcc directory)
# and the NLPQLP directory holding the library with the NLPQLP wrapper function.
# Note that the NLPQLP library files may need to be rebuilt, in which case simply
# enter that directory and use the Makefile.
export LIBGFORTRAN_DIR="/usr/local/Cellar/gcc/10.2.0/lib/gcc/10"
export LIBNLPQLP_DIR="/Users/jm-034232/Documents/McGill/Research/DPGSolver/NLPQLP_Optimizer"  
