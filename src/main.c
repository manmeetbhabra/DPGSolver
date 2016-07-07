// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef TEST

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

#include "test.h"
struct S_TEST TestDB;

/*
 *	Purpose:
 *		Solve the (N)avier-(S)tokes equations (or a subset of the NS equations) using the (D)iscontinuous
 *		(P)etrov-(G)alerkin method.
 *
 *	Comments:
 *
 *	Notation:
 *
 *		DB : (D)ata(B)ase
 *
 *	References:
 *		Demkowicz(2010)_A class of discontinuous Petrov–Galerkin methods. Part I - The transport equation
 *		ToBeModified: Add significant references.
 *
 */

struct S_DB DB;

int main(int nargc, char **argv)
{
	printf("\n\n\n*** Test to see when unrolled mv multiplications break even with BLAS on Guillimin before running"
	       " any large jobs. ***\n\n\n");

	char *fNameOut, *string;
	int  MPIrank, MPIsize;

	struct S_TIME {
		clock_t ts, te;
		double  tt;
	} total, preproc, solving, postproc;

	total.ts = clock();

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // free
	string   = malloc(STRLEN_MIN * sizeof *string);   // free

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), exit(1);

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	preproc.ts = clock();
	initialization(nargc,argv);

	// Preprocessing
	if (!DB.MPIrank)
		printf("Preprocessing:\n\n");

	if (!DB.MPIrank)
		printf("  Set up Parameters\n");
	setup_parameters();

	if (!DB.MPIrank)
		printf("  Set up Mesh\n");
	setup_mesh();

	if (!DB.MPIrank)
		printf("  Set up Operators\n");
	setup_operators();

	if (!DB.MPIrank)
		printf("  Set up Structures\n");
	setup_structures();

	if (!DB.MPIrank)
		printf("  Set up Geometry\n");
	setup_geometry();

	preproc.te = clock();

	// Solving
	solving.ts = clock();
	if (!DB.MPIrank)
		printf("\n\nSolving:\n\n");

	if (!DB.MPIrank)
		printf("  Initializing\n");
	initialize_test_case();

	// Output initial solution to paraview
	strcpy(fNameOut,"SolInitial_");
	                             strcat(fNameOut,DB.TestCase);
	sprintf(string,"%dD",DB.d); strcat(fNameOut,string);
	output_to_paraview(fNameOut);

	if (DB.Restart >= 0) {
		if (!DB.MPIrank)
			printf("  Initializing restarted solution if enabled.\n");
		// Need to ensure that the same proc distribution is used as for lower order solutions.
//		restart_read();
	}

	if (!DB.MPIrank)
		printf("  Nonlinear Iterative Solve\n\n");

	if (strstr(DB.SolverType,"Explicit") != NULL) {
		solver_explicit();
	} else if (strstr(DB.SolverType,"Implicit") != NULL) {
		; //solver_implicit();
	} else {
		printf("Error: Unsupported SolverType in dpg_solver.\n"), exit(1);
	}
	solving.te = clock();

	// Postprocessing
	postproc.ts = clock();
	if (!DB.MPIrank)
		printf("\n\nPostprocessing:\n\n");

	// Output final solution to paraview
	printf("  Output final solution to paraview\n");
	strcpy(fNameOut,"SolFinal_");
	sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
	                               strcat(fNameOut,DB.MeshType);
	sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
	if (DB.Adapt == ADAPT_0)
		sprintf(string,"P%d_",DB.PGlobal), strcat(fNameOut,string);
	output_to_paraview(fNameOut);

	// Compute errors
	if (!DB.MPIrank)
		printf("  Computing errors\n");
	compute_errors();

	postproc.te = clock();

	free(fNameOut);
	free(string);

	memory_free();

	// End MPI and PETSC
	PetscFinalize();

	total.te = clock();

	printf("\n\n\nTotal time       : % .2f s\n\n",(total.te-total.ts)/(double) CLOCKS_PER_SEC);
	printf("  Preprocessing  : % .2f s\n",(preproc.te-preproc.ts)/(double) CLOCKS_PER_SEC);
	printf("  Solving        : % .2f s\n",(solving.te-solving.ts)/(double) CLOCKS_PER_SEC);
	printf("  Postprocessing : % .2f s\n",(postproc.te-postproc.ts)/(double) CLOCKS_PER_SEC);
	printf("\n\n\n");

	return 0;
}

#else // Run if -DTEST is passed as a compilation flag

#include <stdio.h>
#include <time.h>

#include "test.h"
#include "database.h"
#include "functions.h"

/*
 *	Purpose:
 *		Run test functions:
 *			1) Speed comparisons
 *			2) Correctness of implementation
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_DB DB;
struct S_TEST TestDB;

#define TOL 1e-15

int main(void)
{
	clock_t ts, te;

	TestDB.Ntest = 0;
	TestDB.Npass = 0;
	TestDB.Nwarnings = 0;


	printf("\n\nRunning Tests:\n\n\n");
	ts = clock();

	// Implementation tests

	test_imp_array_find_index();
	test_imp_array_norm();
	test_imp_array_sort();
	test_imp_array_swap();

	test_imp_math_factorial();
	test_imp_math_gamma();

	test_imp_matrix_diag();
	test_imp_matrix_identity();
	test_imp_matrix_inverse();
	test_imp_matrix_mm();
	test_imp_convert_to_CSR();

	test_imp_find_periodic_connections();

	test_imp_cubature_TP();
	test_imp_cubature_SI();
	test_imp_cubature_PYR();

	test_imp_basis_TP();
	test_imp_basis_SI();
	test_imp_basis_PYR();
	test_imp_grad_basis_TP();
	test_imp_grad_basis_SI();
	test_imp_grad_basis_PYR();

	test_imp_sum_factorization();
	test_imp_plotting();

	test_imp_fluxes_inviscid();
	test_imp_get_facet_ordering();


	te = clock();

	// Speed tests
//	test_speed_array_swap();
//	test_speed_mm_CTN();


	printf("\n\nRan %d test(s) in %.4f seconds.\n",TestDB.Ntest,(te-ts)/(float)CLOCKS_PER_SEC);

	unsigned int Nfail = TestDB.Ntest - TestDB.Npass;
	if (Nfail > 0) {
		printf("\n\n******** FAILED %d TEST(S) ********\n\n",Nfail);
	} else {
		printf("\nAll tests passed.\n\n");

		if (TestDB.Nwarnings)
			printf("Warnings (%d) were generated while running tests. "
			       "Scroll through test passing list and verify that all is OK.\n\n",TestDB.Nwarnings);
	}

	return 0;
}


#endif // End TEST