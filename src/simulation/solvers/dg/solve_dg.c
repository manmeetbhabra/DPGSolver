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

#include "solve_dg.h"

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_math.h"
#include "petscvec.h"

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

#include "face.h"
#include "face_solver.h"
#include "element.h"
#include "element_solver_dg.h"
#include "volume.h"
#include "volume_solver_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_face_rlhs.h"
#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "compute_face_rlhs_dg.h"
#include "compute_source_rlhs_dg.h"
#include "const_cast.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution_euler.h"
#include "solution_navier_stokes.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Call the common functions when computing the rlhs values for explicit and implicit dg schemes.
static void compute_rlhs_common_dg
	(const struct Simulation*const sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit*const ssi ///< \ref Solver_Storage_Implicit.
	);

/// \brief Scale the rhs terms by the the inverse mass matrices (for explicit schemes).
static void scale_rhs_by_m_inv
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute the maximum value of the rhs term for the first variable.
 *  \return See brief. */
static double compute_max_rhs
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Fill \ref Solver_Storage_Implicit::b with the negative of rhs values.
static void fill_petsc_Vec_b_dg
	(const struct Simulation*const sim,       ///< Defined for \ref compute_rlhs_dg.
	 struct Solver_Storage_Implicit*const ssi ///< Defined for \ref compute_rlhs_dg.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_dg_T.c"

double compute_rhs_dg (const struct Simulation* sim)
{
	compute_rlhs_common_dg(sim,NULL);
	scale_rhs_by_m_inv(sim);

	return compute_max_rhs(sim);
}

double compute_rlhs_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	compute_rlhs_common_dg(sim,ssi);
	fill_petsc_Vec_b_dg(sim,ssi);

	return compute_max_rhs(sim);
}

void set_petsc_Mat_row_col
	(struct Solver_Storage_Implicit*const ssi, const struct Solver_Volume* v_l, const int eq,
	 const struct Solver_Volume* v_r, const int vr)
{
	ssi->row = (int)(v_l->ind_dof+v_l->sol_coef->extents[0]*eq);
	ssi->col = (int)(v_r->ind_dof+v_r->sol_coef->extents[0]*vr);
}

void add_to_petsc_Mat (const struct Solver_Storage_Implicit*const ssi, const struct const_Matrix_d* lhs)
{
	const ptrdiff_t ext_0 = lhs->ext_0,
	                ext_1 = lhs->ext_1;

	PetscInt idxm[ext_0],
	         idxn[ext_1];

	for (int i = 0; i < ext_0; ++i)
		idxm[i] = ssi->row+i;

	for (int i = 0; i < ext_1; ++i)
		idxn[i] = ssi->col+i;

	const PetscScalar*const vv = lhs->data;
	MatSetValues(ssi->A,(PetscInt)ext_0,idxm,(PetscInt)ext_1,idxn,vv,ADD_VALUES);
}

void compute_flux_imbalances_dg (const struct Simulation*const sim)
{
	compute_flux_imbalances_faces_dg(sim);
	compute_flux_imbalances_source_dg(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref scale_rhs_by_m_inv used for non-collocated runs.
static void scale_rhs_by_m_inv_std
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Version of \ref scale_rhs_by_m_inv used for collocated runs.
static void scale_rhs_by_m_inv_col
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute and add diagonal block terms to the LHS matrix in case \ref LHS_CFL_RAMPING was selected.
 *
 *  The addition of these terms corresponds to the use of the backwards (implicit) Euler method with bounded time step,
 *  as opposed to the full Newton method which is the limit of the backwards Euler as the time step goes to infinity.
 *  The 'C'ourant'F'riedrichs'L'ewy number is used to determine the time step and is increased as the residual
 *  decreases.
 *
 *  \todo Add documentation for this procedure.
 */
static void compute_CFL_ramping
	(struct Solver_Storage_Implicit*const ssi, ///< \ref Solver_Storage_Implicit.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

static void compute_rlhs_common_dg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	zero_memory_volumes(sim->volumes);
	compute_grad_coef_dg(sim,sim->volumes,sim->faces);
	compute_volume_rlhs_dg(sim,ssi,sim->volumes);
	compute_face_rlhs_dg(sim,ssi,sim->faces);
	compute_source_rhs_dg(sim);

	compute_CFL_ramping(ssi,sim);
}

static void scale_rhs_by_m_inv (const struct Simulation*const sim)
{
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	assert((test_case->solver_proc == SOLVER_E) ||
	       (test_case->solver_proc == SOLVER_EI));

	if (!sim->collocated)
		scale_rhs_by_m_inv_std(sim);
	else
		scale_rhs_by_m_inv_col(sim);
}

static double compute_max_rhs (const struct Simulation*const sim)
{
	double max_rhs = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

		struct Multiarray_d* rhs = dg_s_vol->rhs;
		double max_rhs_curr = norm_d(rhs->extents[0],rhs->data,"Inf");
		if (max_rhs_curr > max_rhs)
			max_rhs = max_rhs_curr;
	}
	return max_rhs;
}

static void fill_petsc_Vec_b_dg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const int ind_dof        = (int)((struct Solver_Volume*)curr)->ind_dof;
		struct Multiarray_d* rhs = ((struct DG_Solver_Volume*)curr)->rhs;

		const int ni = (int)compute_size(rhs->order,rhs->extents);

		PetscInt    ix[ni];
		PetscScalar y[ni];

		for (int i = 0; i < ni; ++i) {
			ix[i] = ind_dof+i;
			y[i]  = -(rhs->data[i]);
		}

		VecSetValues(ssi->b,ni,ix,y,INSERT_VALUES);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Return the time step relating to the CFL ramping.
 *  \return See brief. */
static double compute_dt_cfl_constrained
	(const double max_rhs,                   ///< The maximum of the RHS terms (negated steady residual).
	 const struct Solver_Volume*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim       ///< \ref Simulation.
	);

static void scale_rhs_by_m_inv_std (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;
		mm_NN1C_overwrite_Multiarray_d(dg_s_vol->m_inv,&dg_s_vol->rhs);
	}
}

static void scale_rhs_by_m_inv_col (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol       = (struct Solver_Volume*) curr;
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

		const struct const_Vector_d jac_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
		scale_Multiarray_by_Vector_d('L',1.0,dg_s_vol->rhs,&jac_det_vc,true);
	}
}

static void compute_CFL_ramping (struct Solver_Storage_Implicit*const ssi, const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;

	if (test_case->lhs_terms != LHS_CFL_RAMPING)
		return;

	const double max_rhs = compute_max_rhs(sim);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		const double dt = compute_dt_cfl_constrained(max_rhs,s_vol,sim);
UNUSED(dt);
//		const double max_wave_speed = compute_max_wave_speed(s_vol,sim);
//		const double max_viscosity  = compute_max_viscosity(s_vol,sim);

//		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

EXIT_ADD_SUPPORT; UNUSED(ssi);
	}
}

// Level 2 ********************************************************************************************************** //

/** \brief Return the minimum measure of the length of the input volume.
 *  \return See brief. */
static double compute_min_length_measure
	(const struct Solver_Volume*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim       ///< \ref Simulation.
	);

/** \brief Return the ratio of the initial max_rhs divided by the current max_rhs.
 *  \return See brief. */
static double compute_max_rhs_ratio
	(const double max_rhs ///< The current maximum rhs value.
	);

static double compute_dt_cfl_constrained
	(const double max_rhs, const struct Solver_Volume*const s_vol, const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;

	struct Multiarray_d* s_coef_b = constructor_s_coef_bezier(s_vol,sim); // destructed

	const ptrdiff_t ext_0 = s_coef_b->extents[0];

	// Max wave speed.
	struct Multiarray_d* V_p_c = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,1}); // destructed
	compute_max_wavespeed(V_p_c,(struct const_Multiarray_d*)s_coef_b,'c');

	const double max_wave_speed = maximum_dd(V_p_c->data,ext_0);
	destructor_Multiarray_d(V_p_c);

	// Max viscosity.
	double max_viscosity = EPS;
	if (test_case->has_2nd_order) {
		struct Multiarray_d* mu = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,1}); // destructed
		compute_viscosity(mu,(struct const_Multiarray_d*)s_coef_b,'c',sim->input_path);

		max_viscosity = maximum_dd(mu->data,ext_0);
		destructor_Multiarray_d(mu);
	}
	destructor_Multiarray_d(s_coef_b);

	// Min length measure.
	double dx = compute_min_length_measure(s_vol,sim);

EXIT_ADD_SUPPORT;

	const double max_rhs_ratio = compute_max_rhs_ratio(max_rhs);
	const double cfl = test_case->cfl_initial * ( max_rhs_ratio < 1.0 ? 1.0 : sqrt(max_rhs_ratio) );

	return cfl*GSL_MIN(dx/max_wave_speed,( !test_case->has_2nd_order ? DBL_MAX : dx*dx/max_viscosity ));
}

// Level 3 ********************************************************************************************************** //

static double compute_min_length_measure (const struct Solver_Volume*const s_vol, const struct Simulation*const sim)
{
	/** In the 1D case, the length of the volume is being returned. In higher dimensional cases, the minimum of the
	 *  lengths of each of the edges is returned.
	 *
	 *  Note that the ideal length measure is the distance from the outgoing characteristic wave to the edge of the
	 *  neighbouring volume (relating to the physics of the CFL constraint). The approximation used here is chosen for
	 *  efficiency/convenience as it can be computed locally to each volume but that this would require that a smaller
	 *  CFL number than the maximum allowable for linear stability be used.
	 */
	double min_length_measure = DBL_MAX;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	const struct Volume*const vol    = (struct Volume*) s_vol;
	const struct const_Element* el   = (struct const_Element*) vol->element;
	const struct Solver_Element* s_e = (struct Solver_Element*) el;

	const bool curved = vol->curved;
	const int p = s_vol->p_ref;
	const struct const_Multiarray_d*const g_coef = s_vol->geom_coef;

	struct Multiarray_d*const v_vals = constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){2,DIM}); // destructed
	if (el->d == 1) {
		const struct Operator*const cv0_vg_vv =
			( curved ? get_Multiarray_Operator(s_e->cv0_vg_vv[curved],(ptrdiff_t[]){0,0,1,p})
			         : get_Multiarray_Operator(s_e->cv0_vg_vv[curved],(ptrdiff_t[]){0,0,1,1}) );

		mm_NN1_Operator_Multiarray_d(cv0_vg_vv,g_coef,v_vals,op_format,2,NULL,NULL);

		assert(v_vals->extents[0] == 2);
		assert(v_vals->extents[1] == 1);
		min_length_measure = fabs(v_vals->data[0]-v_vals->data[1]);
	} else {
		const int n_e = el->n_e;
		for (int le = 0; le < n_e; ++le) {
			const struct Operator*const cv0_vg_ev =
				( curved ? get_Multiarray_Operator(s_e->cv0_vg_ev[curved],(ptrdiff_t[]){le,0,0,1,p})
				         : get_Multiarray_Operator(s_e->cv0_vg_ev[curved],(ptrdiff_t[]){le,0,0,1,1}) );

			mm_NN1_Operator_Multiarray_d(cv0_vg_ev,g_coef,v_vals,op_format,2,NULL,NULL);

			assert(v_vals->extents[0] == 2);
			double*const data[2] = { get_row_Multiarray_d(0,v_vals), get_row_Multiarray_d(1,v_vals), };
			for (int d = 0; d < DIM; ++d)
				data[0][d] -= data[1][d];

			const double edge_len = norm_d(DIM,data[0],"L2");
			if (edge_len < min_length_measure)
				min_length_measure = edge_len;
		}
	}
	destructor_Multiarray_d(v_vals);
	return min_length_measure;
}

static double compute_max_rhs_ratio (const double max_rhs)
{
	static double max_rhs0 = 0.0;
	if (max_rhs0 == 0.0)
		max_rhs0 = max_rhs;
	return max_rhs0/max_rhs;
}
