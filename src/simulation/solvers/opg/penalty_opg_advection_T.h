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
 *  \brief Provides functions used for the penalty term for both right and left-hand side (rlhs) terms of the OPG scheme.
 */

#include "def_templates_penalty_opg_advection.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver.h"

struct Flux_Ref_T;
struct Numerical_Flux_T;
struct Solver_Face_T;
struct Solver_Storage_Implicit;

/** \brief Version of \ref compute_rlhs_opg_f_fptr_T for rhs and lhs terms which adds penalty terms relating to the use
 * of an upwind numerical flux for the linear advection equation. */
void constructor_rlhs_f_test_penalty_advection_upwind_T
	(const struct Flux_Ref_T*const flux_r,         ///< See brief.
	 const struct Numerical_Flux_T*const num_flux, ///< See brief.
	 struct Solver_Face_T*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi      ///< See brief.
	 );

#include "undef_templates_penalty_opg_advection.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver.h"
