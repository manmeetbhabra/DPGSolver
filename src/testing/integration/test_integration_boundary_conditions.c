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

#include <assert.h>
#include <string.h>
#include "petscsys.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"
#include "definitions_test_integration.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support.h"
#include "test_support_multiarray.h"

#include "multiarray.h"

#include "face.h"
#include "face_solver.h"

#include "boundary.h"
#include "computational_elements.h"
#include "compute_face_rlhs.h"
#include "intrusive.h"
#include "numerical_flux.h"
#include "simulation.h"
#include "test_case.h"

#include "complex_multiarray.h"

#include "test_complex_face_solver.h"

#include "test_complex_compute_face_rhs.h"
#include "test_complex_numerical_flux.h"
#include "test_complex_math_functions.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref constructor_Numerical_Flux_Input_data_T, forcing the gradients to be computed using
 *         interpolation if required.
 */
void constructor_Numerical_Flux_Input_data_with_gradients
	(struct Numerical_Flux_Input*const num_flux_i, ///< See brief.
	 const struct Solver_Face*const s_face,        ///< See brief.
	 const struct Simulation*const sim             ///< See brief.
	);

/** \brief Constructor for the real \ref Numerical_Flux_Input_T computed using the complex step method.
 *  \return See brief. */
static struct Numerical_Flux_Input* constructor_Numerical_Flux_Input_cmplx_step
	(const struct Solver_Face*const s_face_r,   ///< Real \ref Solver_Face_T.
	 const struct Solver_Face_c*const s_face_c, ///< Complex \ref Solver_Face_T.
	 struct Simulation* sim[2]                  ///< Real and complex \ref Simulation containers.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the boundary conditions (\ref test_integration_boundary_conditions.c).
 *  \return 0 on success.
 *
 *  The linearizations are verified by comparing with the output when using the complex step method. Details of the
 *  complex step method can be found in Squire et al. \cite Squire1998 and Martins et al. \cite Martins2003.
 */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(nargc == 3,"Invalid number of input arguments");
	const char*const ctrl_name = argv[2];

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p          = int_test_info->p_ref[0],
	          ml         = int_test_info->ml[0],
	          adapt_type = int_test_info->adapt_type;
	assert(adapt_type == ADAPT_0);

	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	const char type_rc[2] = { 'r', 'c', };
	struct Simulation* sim[2] = { NULL, NULL, };
	int ind_rc = -1;

	// real
	ind_rc = 0;
	structor_simulation(&sim[ind_rc],'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc[ind_rc]); // destructed
	set_up_solver_geometry(sim[ind_rc]);
	set_initial_solution(sim[ind_rc]);

	((struct Test_Case*)sim[ind_rc]->test_case_rc->tc)->solver_method_curr = 'i';


	struct Numerical_Flux_Input* num_flux_i = constructor_Numerical_Flux_Input(sim[ind_rc]); // destructed

	// complex
	ind_rc = 1;
	structor_simulation(&sim[ind_rc],'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc[ind_rc]); // destructed
	convert_to_Test_Case_rc(sim[ind_rc],'c');
	set_up_solver_geometry(sim[ind_rc]);
	set_initial_solution_c(sim[ind_rc]);

	((struct Test_Case_c*)sim[ind_rc]->test_case_rc->tc)->solver_method_curr = 'e';

	bool pass_all = true;
	for (struct Intrusive_Link* curr_r = sim[0]->faces->first; curr_r; curr_r = curr_r->next) {
		struct Face* face_r          = (struct Face*) curr_r;
		struct Solver_Face* s_face_r = (struct Solver_Face*) curr_r;

		if (!face_r->boundary)
			continue;

		struct Face* face_c            = NULL;
		struct Solver_Face_c* s_face_c = NULL;
		for (struct Intrusive_Link* curr_c = sim[1]->faces->first; curr_c; curr_c = curr_c->next) {
			face_c = (struct Face*) curr_c;
			if (face_r->index == face_c->index) {
				s_face_c = (struct Solver_Face_c*) curr_c;
				break;
			}
		}

		constructor_Numerical_Flux_Input_data_with_gradients(num_flux_i,s_face_r,sim[0]); // destructed

		struct Numerical_Flux_Input* num_flux_i_cmplx_step =
			constructor_Numerical_Flux_Input_cmplx_step(s_face_r,s_face_c,sim); // destructed

		const bool* c_m = num_flux_i->bv_l.compute_member;

		assert(!c_m[4]); // Add support.
		assert(!c_m[5]); // Add support.
		bool pass        = false;
		const double tol = 1e1*EPS;
		const bool differences[] =
			{ c_m[0] ? diff_const_Multiarray_d(num_flux_i->bv_r.s,    num_flux_i_cmplx_step->bv_r.s,    tol) : false,
			  c_m[1] ? diff_const_Multiarray_d(num_flux_i->bv_r.ds_ds,num_flux_i_cmplx_step->bv_r.ds_ds,tol) : false,
			  c_m[2] ? diff_const_Multiarray_d(num_flux_i->bv_r.g    ,num_flux_i_cmplx_step->bv_r.g    ,tol) : false,
			  c_m[3] ? diff_const_Multiarray_d(num_flux_i->bv_r.dg_dg,num_flux_i_cmplx_step->bv_r.dg_dg,tol) : false,
			};
		const int len = COUNT_OF(differences);
		if (check_diff(len,differences,&pass)) {
			if (differences[0]) print_diff_const_Multiarray_d(num_flux_i->bv_r.s,    num_flux_i_cmplx_step->bv_r.s,    tol);
			if (differences[1]) print_diff_const_Multiarray_d(num_flux_i->bv_r.ds_ds,num_flux_i_cmplx_step->bv_r.ds_ds,tol);
			if (differences[2]) print_diff_const_Multiarray_d(num_flux_i->bv_r.g,    num_flux_i_cmplx_step->bv_r.g,    tol);
			if (differences[3]) print_diff_const_Multiarray_d(num_flux_i->bv_r.dg_dg,num_flux_i_cmplx_step->bv_r.dg_dg,tol);
		}
		char message[STRLEN_MAX];
		sprintf(message,"%s%d%s","boundary_condition_linearization (bc: ",face_r->bc,")");
		expect_condition(pass,message);
		if (!pass)
			pass_all = false;
#if 0
printf("%d %d %d %d %d %d\n",face_r->index,face_r->bc,differences[0],differences[1],differences[2],differences[3]);
print_const_Multiarray_d(num_flux_i->bv_r.ds_ds);
print_diff_const_Multiarray_d(num_flux_i->bv_r.g,    num_flux_i_cmplx_step->bv_r.g,    tol);
print_diff_const_Multiarray_d(num_flux_i->bv_r.dg_dg,num_flux_i_cmplx_step->bv_r.dg_dg,tol);
//print_const_Multiarray_d(num_flux_i->bv_r.dg_dg);
//print_const_Multiarray_d(num_flux_i_cmplx_step->bv_r.dg_dg);
#endif


		destructor_Numerical_Flux_Input_data(num_flux_i);
		destructor_Numerical_Flux_Input_data(num_flux_i_cmplx_step);
		destructor_Numerical_Flux_Input(num_flux_i_cmplx_step);
	}
	assert_condition(pass_all);

	// complex
	ind_rc = 1;
	convert_to_Test_Case_rc(sim[ind_rc],'r');
	structor_simulation(&sim[ind_rc],'d',adapt_type,p,ml,0,0,NULL,type_rc[ind_rc]);

	// real
	ind_rc = 0;
	destructor_Numerical_Flux_Input(num_flux_i);
	structor_simulation(&sim[ind_rc],'d',adapt_type,p,ml,0,0,NULL,type_rc[ind_rc]);

	destructor_Integration_Test_Info(int_test_info);
	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Set the members of the input real \ref mutable_Boundary_Value_T container to zero.
static void set_to_zero_Boundary_Value
	(struct mutable_Boundary_Value*const m_boundary_value, ///< \ref mutable_Boundary_Value_T.
	 const struct Numerical_Flux_Input*const num_flux_i    ///< \ref Numerical_Flux_Input_T.
	);

/// \brief Constructor for the 'l'eft data members of the complex \ref Numerical_Flux_Input_T container.
static void constructor_Numerical_Flux_Input_c_data_members_l
	(struct Numerical_Flux_Input_c*const num_flux_c_i, ///< Complex \ref Numerical_Flux_Input_T.
	 struct Numerical_Flux_Input*const num_flux_i      ///< Real    \ref Numerical_Flux_Input_T.
	);

/// \brief Destructor for the 'l'eft data members of the complex \ref Numerical_Flux_Input_T container.
static void destructor_Numerical_Flux_Input_c_data_members_l
	(struct Numerical_Flux_Input_c*const num_flux_c_i ///< Complex \ref Numerical_Flux_Input_T.
	);

/// \brief Version of \ref constructor_Numerical_Flux_Input_data_T constructing only the 'r'ight data members.
void constructor_Boundary_Value_c_data
	(struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 const struct Solver_Face_c* s_face,        ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Version of \ref destructor_Numerical_Flux_Input_data_T destructing only the 'r'ight data members.
void destructor_Boundary_Value_c_data
	(struct Numerical_Flux_Input_c* num_flux_i ///< See brief.
	);

void constructor_Numerical_Flux_Input_data_with_gradients
	(struct Numerical_Flux_Input*const num_flux_i, const struct Solver_Face*const s_face,
	 const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*)sim->test_case_rc->tc;

	constructor_Boundary_Value_Input_face_fptr
		constructor_Boundary_Value_Input_face_fcl = test_case->constructor_Boundary_Value_Input_face_fcl;
	if (test_case->has_2nd_order) {
		if (constructor_Boundary_Value_Input_face_fcl == constructor_Boundary_Value_Input_face_s_fcl_interp)
			constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_sg_fcl_interp;
	}

	constructor_Boundary_Value_Input_face_fcl(&num_flux_i->bv_l,s_face,sim);                // destructed
	s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

static struct Numerical_Flux_Input* constructor_Numerical_Flux_Input_cmplx_step
	(const struct Solver_Face*const s_face_r, const struct Solver_Face_c*const s_face_c, struct Simulation* sim[2])
{
	struct Numerical_Flux_Input* num_flux_i = constructor_Numerical_Flux_Input(sim[0]); // returned

	constructor_Numerical_Flux_Input_data_with_gradients(num_flux_i,s_face_r,sim[0]); // destructed

	struct mutable_Boundary_Value* m_boundary_value = (struct mutable_Boundary_Value*) &num_flux_i->bv_r;
	set_to_zero_Boundary_Value(m_boundary_value,num_flux_i);

	struct Numerical_Flux_Input_c* num_flux_i_c = constructor_Numerical_Flux_Input_c(sim[1]); // destructed
 	constructor_Numerical_Flux_Input_c_data_members_l(num_flux_i_c,num_flux_i); // destructed

	const bool* c_m = num_flux_i->bv_l.compute_member;

	struct Multiarray_c*const s = (struct Multiarray_c*) num_flux_i_c->bv_l.s;
	const ptrdiff_t n_n  = s->extents[0],
	                n_vr = s->extents[1];

	for (int vr = 0; vr < n_vr; ++vr) {
		double complex*const data_s = get_col_Multiarray_c(vr,s);
		add_to_c(data_s,CX_STEP*I,n_n);
		constructor_Boundary_Value_c_data(num_flux_i_c,s_face_c,sim[1]); // destructed

		// s[NVAR]
		assert(c_m[0] == true);
		if (vr == 0) {
			const double complex*const data_c = get_col_const_Multiarray_c(0,num_flux_i_c->bv_r.s);
			double*const data_r               = get_col_Multiarray_d(0,m_boundary_value->s);
			for (int n = 0; n < n_n*n_vr; ++n)
				data_r[n] = creal(data_c[n]);
		}

		// ds_ds[NVAR,NVAR]
		if (c_m[1]) {
			for (int vr_r = 0; vr_r < n_vr; ++vr_r) {
				const ptrdiff_t ind_c = vr_r,
				                ind_r = vr_r+n_vr*(vr);

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,num_flux_i_c->bv_r.s);
				double*const data_r               = get_col_Multiarray_d(ind_r,m_boundary_value->ds_ds);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = cimag(data_c[n])/CX_STEP;
			}
		}

		destructor_Boundary_Value_c_data(num_flux_i_c);
		add_to_c(data_s,-CX_STEP*I,n_n);
	}

	struct Multiarray_c*const g = (struct Multiarray_c*) num_flux_i_c->bv_l.g;

	if (c_m[2]) {
	for (int dx = 0; dx < DIM; ++dx) {
	for (int vr = 0; vr < n_vr; ++vr) {
		double complex*const data_g = get_col_Multiarray_c(vr+n_vr*(dx),g);
		add_to_c(data_g,CX_STEP*I,n_n);
		constructor_Boundary_Value_c_data(num_flux_i_c,s_face_c,sim[1]); // destructed

		// g[NVAR,DIM]
		if (c_m[2] && vr == 0) {
			const double complex*const data_c = get_col_const_Multiarray_c(0+n_vr*dx,num_flux_i_c->bv_r.g);
			double*const data_r               = get_col_Multiarray_d(0+n_vr*dx,m_boundary_value->g);
			for (int n = 0; n < n_n*n_vr; ++n)
				data_r[n] = creal(data_c[n]);
		}

		// dg_dg[NVAR,NVAR]
		if (c_m[3]) {
			for (int dx_r = 0; dx_r < DIM; ++dx_r) {
			for (int vr_r = 0; vr_r < n_vr; ++vr_r) {
				const ptrdiff_t ind_c = vr_r+n_vr*dx_r,
				                ind_r = vr_r+n_vr*(dx_r+DIM*(vr+n_vr*(dx)));

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,num_flux_i_c->bv_r.g);
				double*const data_r               = get_col_Multiarray_d(ind_r,m_boundary_value->dg_dg);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = cimag(data_c[n])/CX_STEP;
			}}
		}

		assert(!c_m[4]); // Add support.
		assert(!c_m[5]); // Add support.

		destructor_Boundary_Value_c_data(num_flux_i_c);
		add_to_c(data_g,-CX_STEP*I,n_n);
	}}}

 	destructor_Numerical_Flux_Input_c_data_members_l(num_flux_i_c);
	destructor_Numerical_Flux_Input_c(num_flux_i_c);

	return num_flux_i;
}

// Level 1 ********************************************************************************************************** //

static void set_to_zero_Boundary_Value
	(struct mutable_Boundary_Value*const m_boundary_value, const struct Numerical_Flux_Input*const num_flux_i)
{
	const bool*const c_m = num_flux_i->bv_l.compute_member;

	if (c_m[0])
		set_to_value_Multiarray_d(m_boundary_value->s,0.0);
	if (c_m[1])
		set_to_value_Multiarray_d(m_boundary_value->ds_ds,0.0);
	if (c_m[2]) {
		assert(m_boundary_value->g != NULL);
		set_to_value_Multiarray_d(m_boundary_value->g,0.0);
	}
	if (c_m[3]) {
		assert(m_boundary_value->dg_dg != NULL);
		set_to_value_Multiarray_d(m_boundary_value->dg_dg,0.0);
	}
	if (c_m[4])
		EXIT_ADD_SUPPORT;
	if (c_m[5])
		EXIT_ADD_SUPPORT;
}

static void constructor_Numerical_Flux_Input_c_data_members_l
	(struct Numerical_Flux_Input_c*const num_flux_c_i, struct Numerical_Flux_Input*const num_flux_i)
{
	num_flux_c_i->bv_l.normals = constructor_copy_const_Multiarray_d(num_flux_i->bv_l.normals);
	num_flux_c_i->bv_l.xyz     = constructor_copy_const_Multiarray_d(num_flux_i->bv_l.xyz);
	num_flux_c_i->bv_l.s       = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_l.s);
	if (num_flux_i->bv_l.g)
		num_flux_c_i->bv_l.g = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_l.g);
}

static void destructor_Numerical_Flux_Input_c_data_members_l (struct Numerical_Flux_Input_c*const num_flux_c_i)
{
	destructor_const_Multiarray_d(num_flux_c_i->bv_l.normals);
	destructor_const_Multiarray_d(num_flux_c_i->bv_l.xyz);
	destructor_const_Multiarray_c(num_flux_c_i->bv_l.s);
	destructor_conditional_const_Multiarray_c(num_flux_c_i->bv_l.g);
}

void constructor_Boundary_Value_c_data
	(struct Numerical_Flux_Input_c* num_flux_i, const struct Solver_Face_c* s_face, const struct Simulation* sim)
{
	s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Boundary_Value_c_data (struct Numerical_Flux_Input_c* num_flux_i)
{
	destructor_Boundary_Value_c(&num_flux_i->bv_r);
}