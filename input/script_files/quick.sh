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
ARGV="RESULTS_navier_stokes_joukowski_full_dg_2D petsc_options_gmres_tol_1e-4"
#ARGV="RESULTS_euler_joukowski_full_dg_2D petsc_options_gmres_tol_1e-2"
#ARGV="RESULTS_advection_vortex_dg_2d petsc_options_gmres_tol_1e-15"
ARGV="euler/free_stream/TEST_Euler_FreeStream_ParametricMixed2D petsc_options_gmres_tol_1e-4"
ARGV="advection/free_stream/TEST_advection_free_stream_mixed2d_non_conforming petsc_options_gmres_tol_1e-15"
ARGV="advection/vortex/TEST_advection_vortex_slipwall_mixed2d petsc_options_gmres_tol_1e-15"

ARGV="euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D_test petsc_options_gmres_tol_1e-15"
#ARGV="advection/vortex/TEST_advection_vortex_n_cube_mixed2d petsc_options_gmres_tol_1e-15"

ARGV="euler/gaussian_bump/TEST_Euler_GaussianBump_ParametricMixed2D petsc_options_gmres_tol_1e-15"

EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_restart
ARGV="integration/restart/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D"

EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_unit_inverse_mapping
ARGV="quad"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""


# DO NOT MODIFY ANYTHING BELOW THIS LINE

@MPIEXEC@ -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
