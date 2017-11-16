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

#include "compute_face_rlhs.h"

#include <assert.h>

#include "face_solver.h"
#include "element_solver.h"
#include "volume_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

const struct Operator* get_operator__tw0_vt_fc (const int side_index, const struct Solver_Face* s_face)
{
	const struct Face* face           = (struct Face*) s_face;
	const struct Volume* vol          = face->neigh_info[side_index].volume;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;
	const struct Solver_Element* e    = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->tw0_vt_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

const struct Operator* get_operator__cv0_vs_fc (const int side_index, const struct Solver_Face* s_face)
{
	const struct Face* face           = (struct Face*) s_face;
	const struct Volume* vol          = face->neigh_info[side_index].volume;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;
	const struct Solver_Element* e    = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->cv0_vs_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

void permute_Matrix_d_fc
	(struct Matrix_d* data, const char perm_layout, const int side_index_dest, const struct Solver_Face* s_face)
{
	assert(perm_layout == 'R');
	assert(data->layout == 'R');

	const struct const_Vector_i* nc_fc = get_operator__nc_fc(side_index_dest,s_face);
	permute_Matrix_d_V(data,nc_fc);
}

const struct const_Vector_i* get_operator__nc_fc (const int side_index_dest, const struct Solver_Face* s_face)
{
	const struct Neigh_Info* neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];

	struct Volume* vol = neigh_info->volume;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int ind_ord = neigh_info->ind_ord,
	          ind_e   = get_face_element_index((struct Face*)s_face),
	          p_f     = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );

	return get_const_Multiarray_Vector_i(e->nc_fc[curved],(ptrdiff_t[]){ind_ord,ind_e,ind_e,0,0,p_f,p_f});
}

void constructor_Numerical_Flux_Input_data
	(struct Numerical_Flux_Input* num_flux_i, const struct Solver_Face* s_face, const struct Simulation* sim)
{
	struct Test_Case* test_case = sim->test_case;

	test_case->constructor_Boundary_Value_Input_face_fcl(&num_flux_i->bv_l,s_face,sim);     // destructed
	s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Numerical_Flux_Input_data (struct Numerical_Flux_Input* num_flux_i)
{
	destructor_Boundary_Value_Input(&num_flux_i->bv_l);
	destructor_Boundary_Value(&num_flux_i->bv_r);
}

struct Matrix_d* constructor_lhs_f_1
	(const int side_index[2], const struct Numerical_Flux* num_flux, const struct Solver_Face* s_face)
{
	const struct Operator* tw0_vt_fc    = get_operator__tw0_vt_fc(side_index[0],s_face);
	const struct Operator* cv0_vs_fc_op = get_operator__cv0_vs_fc(side_index[1],s_face);

	const struct const_Matrix_d* cv0_vs_fc = cv0_vs_fc_op->op_std;
	bool need_free_cv0 = false;
	if (side_index[0] != side_index[1]) {
		need_free_cv0 = true;
		cv0_vs_fc = constructor_copy_const_Matrix_d(cv0_vs_fc); // destructed
		permute_Matrix_d_fc((struct Matrix_d*)cv0_vs_fc,'R',side_index[0],s_face);
	}

	const struct const_Multiarray_d* dnnf_ds_Ma = num_flux->neigh_info[side_index[1]].dnnf_ds;

	const ptrdiff_t ext_0 = tw0_vt_fc->op_std->ext_0,
	                ext_1 = tw0_vt_fc->op_std->ext_1;
	const int n_eq = dnnf_ds_Ma->extents[1],
	          n_vr = dnnf_ds_Ma->extents[2];

	struct Matrix_d* tw0_nf = constructor_empty_Matrix_d('R',ext_0,ext_1);                         // destructed
	struct Matrix_d* lhs_l  = constructor_empty_Matrix_d('R',ext_0,cv0_vs_fc->ext_1);              // destructed
	struct Matrix_d* lhs    = constructor_empty_Matrix_d('R',n_eq*lhs_l->ext_0,n_vr*lhs_l->ext_1); // returned
	set_to_value_Matrix_d(tw0_nf,0.0);

	struct Vector_d dnnf_ds = { .ext_0 = dnnf_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		const ptrdiff_t ind =
			compute_index_sub_container(dnnf_ds_Ma->order,1,dnnf_ds_Ma->extents,(ptrdiff_t[]){eq,vr});
		dnnf_ds.data = (double*)&dnnf_ds_Ma->data[ind];
		mm_diag_d('R',1.0,0.0,tw0_vt_fc->op_std,(struct const_Vector_d*)&dnnf_ds,tw0_nf,false);

		mm_d('N','N',-1.0,0.0,(struct const_Matrix_d*)tw0_nf,cv0_vs_fc,lhs_l);

		set_block_Matrix_d(lhs,(struct const_Matrix_d*)lhs_l,eq*lhs_l->ext_0,vr*lhs_l->ext_1,'i');
//printf("var, eq: %d %d\n",vr,eq);
//print_Matrix_d(lhs_l);
	}}
	destructor_Matrix_d(tw0_nf);
	destructor_Matrix_d(lhs_l);

	if (need_free_cv0)
		destructor_const_Matrix_d(cv0_vs_fc);

	return lhs;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //