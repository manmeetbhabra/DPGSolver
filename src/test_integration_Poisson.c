// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_Poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACET.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "test_integration_linearization.h"
#include "compute_errors.h"
#include "array_norm.h"
#include "solver_Poisson.h"
#include "finalize_LHS.h"

/*
 *	Purpose:
 *		Test various aspects of the Poisson solver implementation:
 *			1) Linearization
 *			2) Optimal convergence orders
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_Poisson(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	TestDB.TestCase = malloc(STRLEN_MAX * sizeof *(TestDB.TestCase)); // free

	/*
	 *	Input:
	 *
	 *		ToBeModified
	 *
	 *	Expected Output:
	 *
	 *		ToBeModified
	 *
	 */

	unsigned int P, ML, PMin, PMax, MLMin, MLMax;

	PetscBool Symmetric;

	Mat A = NULL, A_cs = NULL, A_csc = NULL;
	Vec b = NULL, b_cs = NULL, b_csc = NULL,
	    x = NULL, x_cs = NULL, x_csc = NULL;

	strcpy(TestDB.TestCase,"Poisson");

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// **************************************************************************************************** //
	// 2D (Mixed TRI/QUAD mesh)
	TestDB.PGlobal = 3;
	TestDB.ML      = 0;
TestDB.PGlobal = 1;
TestDB.ML      = 0;

	strcpy(argvNew[1],"test/Test_Poisson_linearization_mixed2D");

	code_startup(nargc,argvNew,0,1);

	implicit_info_Poisson();

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);

	MatIsSymmetric(A,1e5*EPS,&Symmetric);
//	MatIsSymmetric(A_cs,1e5*EPS,&Symmetric);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf")     < 1e2*EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf") < 1e2*EPS &&
	    Symmetric)
		pass = 1, TestDB.Npass++;
	else
		printf("%e %e %d\n",PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf"),
		                    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf"),Symmetric);

	//     0         10        20        30        40        50
	printf("Linearization Poisson (2D - Mixed):              ");
	test_print(pass);
//	EXIT_MSG;

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();

/*
	// **************************************************************************************************** //
	// 3D (Mixed TET/PYR mesh)
	TestDB.PGlobal = 2;
	TestDB.ML      = 0;

//	strcpy(argvNew[1],"test/Test_Poisson_linearization_mixed3D_TP");
	strcpy(argvNew[1],"test/Test_Poisson_linearization_mixed3D_HW");

	code_startup(nargc,argvNew,0,1);

	implicit_info_Poisson();

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);

	MatIsSymmetric(A,1e5*EPS,&Symmetric);
//	MatIsSymmetric(A_cs,1e5*EPS,&Symmetric);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf")     < 1e2*EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf") < 1e2*EPS &&
	    Symmetric)
		pass = 1, TestDB.Npass++;
	else
		printf("%e %e %d\n",PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf"),
		                    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf"),Symmetric);

	//     0         10        20        30        40        50
	printf("Linearization Poisson (3D - Mixed TET/PYR):      ");
	test_print(pass);
//	EXIT_MSG;

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();
*/

	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
	strcpy(argvNew[1],"test/Test_Poisson_linearization_mixed3D_TP");
//	strcpy(argvNew[1],"test/Test_Poisson_linearization_mixed3D_HW");

/* Getting incorrect orders on TET mesh even with straight elements. Getting correct orders on structed Hex meshes.
 * Getting incorrect orders on unstructured Hex meshes with straight elements. All curved element trials are giving
 * suboptimal convergence. All results above using only Dirichlet BCs.
 *	-> Potentially have a bug with TET treatment and likely the parametrization using cube_to_sphere may introduce some
 *	problems. Also try on the unstructured hex mesh without using the cube_to_sphere projection and see if optimal
 *	orders are obtained.
 */
	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// Convergence orders
	PMin = 1;  PMax = 3;
	MLMin = 0; MLMax = 4;

	for (P = PMin; P <= PMax; P++) {
	for (ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		code_startup(nargc,argvNew,0,1);

		solver_Poisson();
		compute_errors_global();

		if (P == PMax && ML == MLMax)
			check_convergence_orders(MLMin,MLMax,PMin,PMax,&pass);

		code_cleanup();
	}}
	// test with all fluxes as they are implemented (ToBeModified)

	printf("Convergence Orders - Poisson (2D - TRI  ):       ");
	test_print(pass);





	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(TestDB.TestCase);
}
