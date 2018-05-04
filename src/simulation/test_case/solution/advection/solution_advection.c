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

#include "solution_advection.h"

#include <math.h>

#include "macros.h"
#include "definitions_solution.h"

#include "multiarray.h"
#include "vector.h"

#include "boundary.h"
#include "compute_error.h"
#include "compute_error_advection.h"
#include "const_cast.h"
#include "file_processing.h"
#include "flux_advection.h"
#include "math_functions.h"
#include "numerical_flux_advection.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

#include "solution_advection_default.h"
#include "free_stream_advection/solution_free_stream_advection.h"
#include "peterson/solution_peterson.h"
#include "vortex_advection/solution_vortex_advection.h"
#include "hyperbolic_tan/solution_hyperbolic_tan.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_advection_T.c"

struct Sol_Data__Advection get_sol_data_advection ( )
{
	static bool need_input = true;

	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(&sol_data);
	}

	return sol_data;
}

void read_data_advection (struct Sol_Data__Advection*const sol_data)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed

	int advection_type = 0;

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_convert_i(line,"advection_type",&advection_type,&count_found);
		if (strstr(line,"u_scale"))
			read_skip_d_1(line,1,&sol_data->u_scale,1);
		if (strstr(line,"u_coef_polynomial4"))
			read_skip_d_1(line,1,sol_data->u_coef_polynomial4,5);
		if (strstr(line,"b_coef_polynomial_odd"))
			sol_data->b_coef_polynomial_odd = constructor_file_const_Vector_d(input_file,true); // leaked
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");

	switch (advection_type) {
		case ADVECTION_TYPE_CONST:       sol_data->compute_b_adv = compute_b_adv_constant;    break;
		case ADVECTION_TYPE_VORTEX:      sol_data->compute_b_adv = compute_b_adv_vortex;      break;
		case ADVECTION_TYPE_VORTEX_POLY: sol_data->compute_b_adv = compute_b_adv_vortex_poly; break;
		default:                         EXIT_ERROR("Unsupported: %d\n",advection_type);      break;
	}
}

const double* compute_b_adv_constant (const double*const xyz)
{
	UNUSED(xyz);
	static bool need_input = true;
	static double b_adv[DIM] = {0,};

	if (need_input) {
		need_input = false;

		const int count_to_find = 1;

		FILE* input_file = fopen_input('s',NULL,NULL); // closed

		int count_found = 0;
		char line[STRLEN_MAX];
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"b_adv")) {
				read_skip_d_1(line,1,b_adv,DIM);
				++count_found;
			}
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	return b_adv;
}

const double* compute_b_adv_vortex (const double*const xyz)
{
	static bool need_input = true;
	static double b_mag = 0;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;

		FILE* input_file = fopen_input('s',NULL,NULL); // closed

		int count_found = 0;
		char line[STRLEN_MAX];
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"b_magnitude")) {
				read_skip_d_1(line,1,&b_mag,1);
				++count_found;
			}
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	static double b_adv[DIM] = {0,};
	assert(DIM == 2);

	const double x = xyz[0],
	             y = xyz[1],
	             t = atan2(y,x);

	IF_DIM_GE_1( b_adv[0] =  b_mag*sin(t); ) // Note:  sin(atan2(y,x)) ==  y/sqrt(x^2+y^2).
	IF_DIM_GE_2( b_adv[1] = -b_mag*cos(t); ) //       -cos(atan2(y,x)) == -x/sqrt(x^2+y^2).

	return b_adv;
}

const double* compute_b_adv_vortex_poly (const double*const xyz)
{
	assert(DIM == 2);

	const struct Sol_Data__Advection sol_data = get_sol_data_advection();
	const ptrdiff_t n_coef  = sol_data.b_coef_polynomial_odd->ext_0;
	const double*const coef = sol_data.b_coef_polynomial_odd->data;
	assert(n_coef > 0);
	assert(coef[0] != 0.0);

	if (n_coef > 1)
		EXIT_ERROR("The advection velocity is no longer divergence free in this case.\n");
	// Based on currenlty only being usable for n_coef = 1, there is nothing added by using this function over
	// compute_b_adv_vortex... Possibly investigate higher-order divergence free, rotating field in future.

/// \todo Make the exponent flexible here but remove the flexibility in the polynomial order.
	const double xyz_norm = pow(xyz[0]*xyz[0]+xyz[1]*xyz[1],0.0),
	             x        = xyz[0],
	             y        = xyz[1];
UNUSED(x);

	static double b_adv[DIM] = {0,};
	for (int d = 0; d < DIM; ++d)
		b_adv[d] = 0.0;

	for (int i = 0; i < n_coef; ++i) {
/// \todo Delete non-vortex and place in separate functions if used in future. Also remove non-divergence free options.
// vortex
		IF_DIM_GE_1( b_adv[0] += coef[i]*pow(y/xyz_norm,2*i+1) );
		IF_DIM_GE_2( b_adv[1] -= coef[i]*pow(x/xyz_norm,2*i+1) );
// doublet
//b_adv[0] =  (y*y-x*x)/pow(x*x+y*y,2); b_adv[1] = -(2*x*y)/pow(x*x+y*y,2);
// source
//b_adv[0] = x/(x*x+y*y); b_adv[1] = y/(x*x+y*y);
// other
//b_adv[0] = -x*x; b_adv[1] = 2*x*y;
//b_adv[0] = 2*x*y; //b_adv[1] = -y*y;
//b_adv[0] = x*(3*y*y-x*x); b_adv[1] = y*(3*x*x-y*y);
	}

	return b_adv;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
