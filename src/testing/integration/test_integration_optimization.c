
/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */

#include "petscsys.h"

#include "macros.h"

#include "test_base.h"
#include "test_integration.h"
#include "definitions_adaptation.h"
#include "test_integration_convergence_support.h"

#include "solve.h"
#include "definitions_visualization.h"
#include "visualization.h"
#include "compute_error.h"
#include "restart_writers.h"
#include "test_case.h"

#include "core.h"

#include "optimize.h"

// Static function declarations ************************************************************************************* //

/** \brief 	Run the given test case. Here, we will first run the flow solver 
 *	on the initial mesh. Then, we will run the optimization routine.
 */
static void run_test_case(
	int argc, ///< Standard.
	char** argv ///< Standard
	);


// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the solution convergence (\ref test_integration_run.c).
 *  \return 0 on success (when the solution converges).
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{

	// Ensure that this is a test case with only one command line argument
	assert_condition_message((argc == 3) || (argc == 4),"Invalid number of input arguments");

	// Initialize PETSc
	const char* petsc_options_name = set_petsc_options_name(argv[2]);
	PetscInitialize(&argc,&argv,petsc_options_name,PETSC_NULL);

	run_test_case(argc,argv);

	// Finalize PETSc
	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //

static void run_test_case(int argc, char** argv) {

	UNUSED(argc);

	// Read the case parameters from the .ctrl file. Set the parameters needed for the 
	// simulation object constructor and then construct the simulation object.
	const char* ctrl_name = argv[1];
	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	const char type_rc = 'r';

	bool ignore_static = false;

	int ml_prev = ml_ref[0]-1,
    	p_prev  = p_ref[0]-1,
    	ml 		= ml_ref[0],
    	p 		= p_ref[0];

	const int adapt_type = int_test_info->adapt_type;
	
	// To correct the ctrl file name as mesh level and p must be added
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,true,ctrl_name);

	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,ignore_static);
	solve_for_solution(sim);  // Solve flow over initial shape

	if (outputting_restart())
		output_restart(sim);

	optimize(sim);  // Perform the shape optimization

	// Output the results in graphical format
	output_visualization(sim,VIS_GEOM_EDGES);
	output_visualization(sim,VIS_GEOM_VOLUMES);
	output_visualization(sim,VIS_NORMALS);
	output_visualization(sim,VIS_SOLUTION);

	// Compute the error values
	output_error(sim);

	// Destroy Allocated Structures
	structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,type_rc,ignore_static);
	destructor_Integration_Test_Info(int_test_info);
}

