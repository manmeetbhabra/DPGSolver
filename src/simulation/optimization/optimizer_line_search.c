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

#include "optimizer_line_search.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "macros.h"

#include "simulation.h"

#include "matrix.h"
#include "matrix_constructors.h"
#include "matrix_math.h"
#include "multiarray.h"
#include "multiarray_constructors.h"
#include "multiarray_math.h"

#include "geometry.h"
#include "geometry_parametric.h"
#include "solve_implicit.h"

#include "math_functions.h"
#include "file_processing.h"

#include "optimization_case.h"
#include "adjoint.h"
#include "sensitivities.h"
#include "gradient.h"
#include "output_progress.h"


// Static function declarations ************************************************************************************* //


/** \brief Construct the data structure to hold the optimization line search information
 */
static struct Optimizer_Line_Search_Data* constructor_Optimizer_Line_Search_Data(
	struct Optimization_Case* optimization_case ///< Consult optimization_case.h
	);


/** \brief Destruct the data structure that holds the optimization line search information
 */
static void destructor_Optimizer_Line_Search_Data(
	struct Optimizer_Line_Search_Data* optimizer_line_search_data ///< Consult optimizer_line_search.h
	);


/** \brief 	Read the required data from the optimization.geo file and load it into the 
 * 	optimizer_line_search_data data structure.
 */
static void read_optimization_data(
	struct Optimizer_Line_Search_Data* optimizer_line_search_data ///< The data structure holding the optimizer line search data
	);


/** \brief Use the BFGS method to compute the search direction for the optimization. Then, use
 *	the search direction to compute the step length using the backtracking 
 *	algorithm.
 */
static void BFGS_minimizer(
	struct Optimization_Case *optimization_case, ///< Standard. Consult optimization_case.h
	struct Optimizer_Line_Search_Data *optimizer_line_search_data, ///< The optimizer_line_search_data data structure (stores BFGS related information)
	struct Gradient_Data *gradient_data, ///< Standard. Consult gradient.h (holds the gradient data)
	int design_iteration ///< The design iteration that is currently being processed
	);


/** \brief Perform the back tracking line search and ensure that the first Wolfe condition
 *	is satisfied. This function will exit once the optimal alpha (step length) is computed.
 *	Upon exiting the function, the geometry will have been updated to the desired
 *	shape for the given design iteration.
 * 
 * \return Alpha value used as the step length
 */
static double backtracking_line_search(
	struct Optimization_Case* optimization_case, ///< Standard. Consult optimization_case.h
	double* p_k, ///< The search direction for this iteration (kth iteration) as an array of doubles with num_design_pts_dofs values
	double* grad_f_k, ///< The gradient of the objective function for this iteration (kth iteration)
	double alpha, ///< The initial step length value
	double wolfe_condition_c, ///< The wolfe condition c value 
	double wolfe_condition_rho, ///< The wolfe condition rho value (how much to scale the search direction during the backtracking) 
	double wolfe_condition_alpha_min ///< The minimum step length to use
	);


/** \brief Traverse along the search direction p_k by the step length amount alpha. Modify the 
 *	control point positions.
 * 		
 * 		x_{k+1} = x_k + alpha * p_k
 *
 *	x_{k+1} = Modified design control points
 *	x_k = current design control points
 */
static void update_design_points(
	struct Optimization_Case* optimization_case, ///< Standard. Consult optimization_case.h
	double *p_k, ///< The search direction. Should have num_design_pts_dofs values
	double alpha ///< The step length
	);


// Interface functions ********************************************************************************************** //

void optimizer_line_search(struct Optimization_Case* optimization_case){

	// ================================
	//         Preprocessing
	// ================================
	
	struct Simulation *sim = optimization_case->sim;

	FILE *fp = constructor_optimization_progress_file(optimization_case);
	
	struct Optimizer_Line_Search_Data *optimizer_line_search_data = constructor_Optimizer_Line_Search_Data(optimization_case);

	struct Adjoint_Data *adjoint_data 		= constructor_Adjoint_Data(optimization_case);
	struct Sensitivity_Data *sensivity_data = constructor_Sensitivity_Data(optimization_case);
	struct Gradient_Data *gradient_data 	= constructor_Gradient_Data(adjoint_data, sensivity_data);

	// Keep track of the time taken for the optimization and number of design iterations
	clock_t start_t, iter_t;
	double t_elapse;
	start_t = clock();

	int design_iteration = 0;


	// ================================
	//      Optimization Routine
	// ================================


	while(true){

		// Solve the Adjoint equation
		setup_adjoint(adjoint_data, optimization_case->sim, optimization_case->sim_c,
			optimization_case->objective_function, optimization_case->objective_function_c);
		solve_adjoint(adjoint_data, optimization_case->sim);


		// Compute the Sensitivities
		compute_sensitivities(sensivity_data, optimization_case, optimization_case->objective_function,
			optimization_case->objective_function_c);


		// Compute the gradient using the sensitivities and the adjoint
		compute_gradient(gradient_data, adjoint_data, sensivity_data);


// TESTING
test_brute_force_gradient(optimization_case, optimization_case->objective_function, gradient_data->Gradient->data);
exit(0);


		// Keep track of the progress
		double L2_grad = norm_d(gradient_data->Gradient->ext_1, gradient_data->Gradient->data, "L2");
		double objective_func_value = optimization_case->objective_function(sim);
		iter_t = clock();
		t_elapse = (double)(iter_t - start_t)/(CLOCKS_PER_SEC);
		
		progress_file_add_information(fp, optimization_case, L2_grad, objective_func_value, 
			t_elapse, design_iteration, true);


		// Minimize the objective function
		if (strstr(optimization_case->optimizer_spec, "BFGS"))
			BFGS_minimizer(optimization_case, optimizer_line_search_data, gradient_data, design_iteration);
		else
			EXIT_UNSUPPORTED;


		// With the updated shape and flow, copy the data from the real to the complex structures
		// in order to perform the complex step at the next design iteration
		copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);
		destructor_Multiarray_c(optimization_case->geo_data.control_points_c);
		optimization_case->geo_data.control_points_c = 
			constructor_copy_Multiarray_c_Multiarray_d(optimization_case->geo_data.control_points);


		// Exit condition
		if (L2_grad 				< 	optimizer_line_search_data->exit_L2_norm_gradient 		|| 
			objective_func_value 	< 	optimizer_line_search_data->exit_objective_value 		||
			design_iteration 		>= 	optimizer_line_search_data->exit_max_design_iterations)
			break;

		design_iteration++;
	}

	// ================================
	//         Postprocessing
	// ================================

	destructor_optimization_progress_file(fp);

	// Destruct allocated data structures:
	destructor_Adjoint_Data(adjoint_data);
	destructor_Sensitivity_Data(sensivity_data);
	destructor_Gradient_Data(gradient_data);
	destructor_Optimizer_Line_Search_Data(optimizer_line_search_data);
}


// Static functions ************************************************************************************************* //


static struct Optimizer_Line_Search_Data* constructor_Optimizer_Line_Search_Data(struct Optimization_Case* optimization_case){


	struct Optimizer_Line_Search_Data* optimizer_line_search_data = calloc(1,sizeof *optimizer_line_search_data); // returned


	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;

	// Smoothing Matrix:
	// Create a laplacian smoothing matrix with epsilon = 0.6
	optimizer_line_search_data->Laplacian_Matrix = constructor_empty_Matrix_d('C', num_design_pts_dofs, num_design_pts_dofs);
	
	set_to_value_Matrix_d(optimizer_line_search_data->Laplacian_Matrix, 0.0);
	double eps_val = 0.6;

	for (int i = 0; i < num_design_pts_dofs; i++){

		get_col_Matrix_d(i, optimizer_line_search_data->Laplacian_Matrix)[i] = -1.0 * eps_val;
		if (i == 0){
			// First row
			get_col_Matrix_d(i+1, optimizer_line_search_data->Laplacian_Matrix)[i] = 0.5 * eps_val;
		} else if (i == num_design_pts_dofs-1){
			// Last row
			get_col_Matrix_d(i-1, optimizer_line_search_data->Laplacian_Matrix)[i] = 0.5 * eps_val;
		} else{
			// Middle row
			get_col_Matrix_d(i+1, optimizer_line_search_data->Laplacian_Matrix)[i] = 0.5 * eps_val;
			get_col_Matrix_d(i-1, optimizer_line_search_data->Laplacian_Matrix)[i] = 0.5 * eps_val;
		}
	}

	// BFGS Data Structures

	optimizer_line_search_data->B_k_inv = constructor_empty_Matrix_d('C', num_design_pts_dofs, num_design_pts_dofs);
	optimizer_line_search_data->s_k = constructor_empty_Matrix_d('C', num_design_pts_dofs, 1);
	optimizer_line_search_data->grad_f_k = constructor_empty_Matrix_d('C', num_design_pts_dofs, 1);

	// - Set the B_kMin1_inv matrix to be an identity matrix first
	set_to_value_Matrix_d(optimizer_line_search_data->B_k_inv, 0.0);
	for(int i = 0; i < num_design_pts_dofs; i++)
		get_col_Matrix_d(i, optimizer_line_search_data->B_k_inv)[i] = 1.0;


	// Read the optimization parameters
	read_optimization_data(optimizer_line_search_data);


	return optimizer_line_search_data;
}


static void destructor_Optimizer_Line_Search_Data(struct Optimizer_Line_Search_Data* optimizer_line_search_data){

	// Smoothing Matrix

	destructor_Matrix_d(optimizer_line_search_data->Laplacian_Matrix);

	// BFGS Data Structures

	destructor_Matrix_d(optimizer_line_search_data->B_k_inv);
	destructor_Matrix_d(optimizer_line_search_data->s_k);
	destructor_Matrix_d(optimizer_line_search_data->grad_f_k);

	free((void*)optimizer_line_search_data);

}


static void read_optimization_data(struct Optimizer_Line_Search_Data* optimizer_line_search_data){

	// Get the file pointer to the optimization file
	FILE* input_file = fopen_input('o',NULL,NULL); // closed
	char line[STRLEN_MAX];

	// Read in the information from the file
	while (fgets(line,sizeof(line),input_file)) {

		// Any line with a comment flag should be skipped
		if (strstr(line, "//"))
			continue;

		if (strstr(line, "exit_L2_norm_gradient")) 		read_skip_d_1(line, 1, &optimizer_line_search_data->exit_L2_norm_gradient, 1);
		if (strstr(line, "exit_objective_value")) 		read_skip_d_1(line, 1, &optimizer_line_search_data->exit_objective_value, 1);
		if (strstr(line, "exit_max_design_iterations")) read_skip_i_1(line, 1, &optimizer_line_search_data->exit_max_design_iterations, 1);

		if (strstr(line, "step_size_init")) 			read_skip_d_1(line, 1, &optimizer_line_search_data->step_size_init, 1);
		if (strstr(line, "wolfe_condition_c")) 			read_skip_d_1(line, 1, &optimizer_line_search_data->wolfe_condition_c, 1);
		if (strstr(line, "wolfe_condition_rho")) 		read_skip_d_1(line, 1, &optimizer_line_search_data->wolfe_condition_rho, 1);
		if (strstr(line, "wolfe_condition_alpha_min")) 	read_skip_d_1(line, 1, &optimizer_line_search_data->wolfe_condition_alpha_min, 1);
		if (strstr(line, "max_norm_search_vector")) 	read_skip_d_1(line, 1, &optimizer_line_search_data->max_norm_search_vector, 1);
	}

	fclose(input_file);
}


static void BFGS_minimizer(struct Optimization_Case *optimization_case, struct Optimizer_Line_Search_Data *optimizer_line_search_data,
	struct Gradient_Data *gradient_data, int design_iteration){

	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;

	struct Matrix_d *p_k;  // search direction
	p_k = constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free

	if(design_iteration > 0){
		
		// Compute the Hessian approximation using the information for the previous
		// design step (k) to find the information for this design step (k+1)

		struct Matrix_d *B_k_inv 	= optimizer_line_search_data->B_k_inv;
		struct Matrix_d *s_k 		= optimizer_line_search_data->s_k;
		struct Matrix_d *grad_f_k 	= optimizer_line_search_data->grad_f_k;

		double skT_dot_yk, ykT_dot_BKInv_dot_yk; 

		struct Matrix_d *y_k 					= constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free
		struct Matrix_d *BKInv_dot_yk 			= constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free
		struct Matrix_d *yk_dot_skT 			= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *BKInv_dot_yk_dot_skT 	= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *sk_dot_ykT 			= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *sk_dot_ykT_dot_BKInv 	= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *sk_dot_skT 			= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		
		// Compute y_k
		// y_k = grad_f_kPlus1 - grad_f_k
		//  - NOTE: gradient_data->Gradient holds the gradient for this design step (k+1)
		// 		so it is equivalent to grad_f_kPlus1
		for (int i = 0; i < num_design_pts_dofs; i++)
			y_k->data[i] = gradient_data->Gradient->data[i] - grad_f_k->data[i];


		// Compute skT_dot_yk
		skT_dot_yk = 0;	
		for(int i = 0; i < num_design_pts_dofs; i++){
			skT_dot_yk += s_k->data[i] * y_k->data[i];
		}


		// Compute ykT_dot_BKInv_dot_yk
		mm_d('N', 'N', 1.0, 0.0, 
			(const struct const_Matrix_d*const)B_k_inv, 
			(const struct const_Matrix_d*const)y_k, 
			BKInv_dot_yk);

		ykT_dot_BKInv_dot_yk = 0;
		for(int i = 0; i < num_design_pts_dofs; i++){
			ykT_dot_BKInv_dot_yk += y_k->data[i] * BKInv_dot_yk->data[i];
		}


		// Compute sk_dot_skT
		mm_d('N', 'T', 1.0, 0.0, 
			(const struct const_Matrix_d*const)s_k, 
			(const struct const_Matrix_d*const)s_k, 
			sk_dot_skT);


		// Compute BKInv_dot_yk_dot_skT
		mm_d('N', 'T', 1.0, 0.0, 
			(const struct const_Matrix_d*const)y_k, 
			(const struct const_Matrix_d*const)s_k, 
			yk_dot_skT);

		mm_d('N', 'N', 1.0, 0.0, 
			(const struct const_Matrix_d*const)B_k_inv, 
			(const struct const_Matrix_d*const)yk_dot_skT, 
			BKInv_dot_yk_dot_skT);


		// Compute sk_dot_ykT_dot_BKInv
		mm_d('N', 'T', 1.0, 0.0, 
			(const struct const_Matrix_d*const)s_k, 
			(const struct const_Matrix_d*const)y_k, 
			sk_dot_ykT);

		mm_d('N', 'N', 1.0, 0.0, 
			(const struct const_Matrix_d*const)sk_dot_ykT, 
			(const struct const_Matrix_d*const)B_k_inv, 
			sk_dot_ykT_dot_BKInv);


		// Compute B_inv_kPlus1. Place these values into B_k_inv since these new
		// values will form the inverse hessian for this design step coming up
		for (int i = 0; i < num_design_pts_dofs * num_design_pts_dofs; i++){

			B_k_inv->data[i] = B_k_inv->data[i]
				+ ((skT_dot_yk + ykT_dot_BKInv_dot_yk)/(skT_dot_yk*skT_dot_yk)) * sk_dot_skT->data[i] 
				- (1./skT_dot_yk) * (BKInv_dot_yk_dot_skT->data[i] + sk_dot_ykT_dot_BKInv->data[i]);
		}

		// Free allocated structures
		destructor_Matrix_d(y_k);
		destructor_Matrix_d(BKInv_dot_yk);
		destructor_Matrix_d(yk_dot_skT);
		destructor_Matrix_d(BKInv_dot_yk_dot_skT);
		destructor_Matrix_d(sk_dot_ykT);
		destructor_Matrix_d(sk_dot_ykT_dot_BKInv);
		destructor_Matrix_d(sk_dot_skT);

	}

	// NOTE: If on the first design iteration (design_iteration = 0) then B_k_inv is 
	// set to identity (in the constructor). Therefore, p_k will be just a simple
	// steepest descent.

	// Now that the Hessian for this design step has been found, store the data 
	// for this kth step so that we can use it to get the next steps approximation

	// Save grad_f_k so that B_kPlus1_inv can be found (next design iteration)
	for (int i = 0; i < num_design_pts_dofs; i++)
		optimizer_line_search_data->grad_f_k->data[i] = gradient_data->Gradient->data[i];


	// p_k = - B_k_inv * grad_f_k;
	mm_d('N', 'N', -1.0, 0.0,
		(const struct const_Matrix_d*const)optimizer_line_search_data->B_k_inv,
		(const struct const_Matrix_d*const)optimizer_line_search_data->grad_f_k,
		p_k);

	// TEMPORARY ADDITIONs:

	// Make the max norm (euclidean) of the P vector MAX_NORM_P
	double norm_P = 0.0;
	for (int i = 0; i < num_design_pts_dofs; i++){
		norm_P += p_k->data[i] * p_k->data[i];
	}
	norm_P = sqrt(norm_P);
	
	if (norm_P > optimizer_line_search_data->max_norm_search_vector){
		printf("Scale -> norm_P : %e \n", norm_P);
		for (int i = 0; i < num_design_pts_dofs; i++){
			p_k->data[i] = ((optimizer_line_search_data->max_norm_search_vector)/norm_P) * p_k->data[i];
		}	
	}


	// Perform the back tracking line search
	double alpha_k = backtracking_line_search(
		optimization_case, 
		p_k->data, 
		optimizer_line_search_data->grad_f_k->data, 
		optimizer_line_search_data->step_size_init, 
		optimizer_line_search_data->wolfe_condition_c,
		optimizer_line_search_data->wolfe_condition_rho, 
		optimizer_line_search_data->wolfe_condition_alpha_min);


	// Save s_k = alpha_k * p_k so that B_kPlus1_inv can be found (next design iteration)
	for (int i = 0; i < num_design_pts_dofs; i++)
		optimizer_line_search_data->s_k->data[i] = alpha_k * p_k->data[i];


	// Free allocated data structures:
	destructor_Matrix_d(p_k);

}


static double backtracking_line_search(struct Optimization_Case* optimization_case, double* p_k, double* grad_f_k,
	double alpha, double wolfe_condition_c, double wolfe_condition_rho, double wolfe_condition_alpha_min){

	printf("Second Wolfe Condition Check needs to be fixed here. In general just use NLPQLP to optimize until fized\n");
	exit(0);

	struct Simulation *sim = optimization_case->sim;

	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;

	// Compute the dot product of the search direction and gradient
	double p_k_dot_grad_f_k = 0; 
	for (int i = 0; i < num_design_pts_dofs; i++)
		p_k_dot_grad_f_k += p_k[i]*grad_f_k[i];


	// The value of the objective function (f) at x_k and x_k + alpha*p_k
	double f_x_k, f_x_k_plus_alpha_p_k;
	f_x_k = optimization_case->objective_function(sim);

	while (true){

		// Move along the search direction. Update the geometry and recompute the metric terms
		update_design_points(optimization_case, p_k, alpha);
		update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)optimization_case->geo_data.control_points);
		set_up_solver_geometry(sim);

		// Solve the flow on the updated geometry
		solve_implicit(sim);

		// Check Wolfe Condition
		f_x_k_plus_alpha_p_k = optimization_case->objective_function(sim);
		if (f_x_k_plus_alpha_p_k <= (f_x_k + wolfe_condition_c*alpha*p_k_dot_grad_f_k))
			break;

		// If the step length was below a constant minimum value stop the line search
		if (alpha <= wolfe_condition_alpha_min){
			printf("\n\n EXITING BACKTRACK -> WOLFE_CONDITION_ALPHA_MIN \n\n");
			break;
		}

		// Step length not valid (and not too small), so revert the geometry design points to their previous location
		update_design_points(optimization_case, p_k, -1.0*alpha);

		// Decrease step length by factor rho
		printf("\n\n (%e, %e) Alpha : %e -> Alpha : %e \n\n", f_x_k_plus_alpha_p_k, 
			f_x_k + wolfe_condition_c*alpha*p_k_dot_grad_f_k, alpha, alpha*wolfe_condition_rho);  // monitor progress
		alpha *= wolfe_condition_rho;
	}

	printf("\n\n Completed backtracking_line_search\n\n");  // monitor progress

	return alpha;
}


static void update_design_points(struct Optimization_Case* optimization_case, double *p_k, double alpha){

	struct Multiarray_d *control_points = optimization_case->geo_data.control_points;
	struct Multiarray_i *ctrl_pts_optimization = optimization_case->geo_data.control_points_optimization;
	int n_pts_optimization = (int)ctrl_pts_optimization->extents[0];

	int control_pt_index;
	int p_index = 0;

	for (int i = 0; i < n_pts_optimization; i++){
		// Loop over the design optimization points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_optimization)[i])
				continue;

			// Modify the control point position
			control_pt_index = get_col_Multiarray_i(0, ctrl_pts_optimization)[i];
			get_col_Multiarray_d(j-1, control_points)[control_pt_index] += alpha * p_k[p_index++];
		}
	}

	// Monitor the progress
	output_NURBS_patch_information(optimization_case);
}
