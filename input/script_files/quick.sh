#!/bin/bash
#PBS -A rck-371-aa
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=12
#PBS -q hb
#PBS -o outputfile
#PBS -e errorfile
#PBS -V
#PBS -N pz

# Executable and command line arguments
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_convergence
ARGV="p2-3/RESULTS_euler_joukowski_full_dg_2D_ar4_super_p_le_1 petsc_options_gmres_r240 euler_joukowski_full_dg_2D_ar4_super_p_le_1"
ARGV="RESULTS_euler_joukowski_full_dg_2D petsc_options_gmres_r120"
#EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_fe_init_refined
#ARGV="extern_mesh/TEST_h_refined_2d_tri_h1+uniform"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""


# DO NOT MODIFY ANYTHING BELOW THIS LINE

@MPIEXEC@ -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
