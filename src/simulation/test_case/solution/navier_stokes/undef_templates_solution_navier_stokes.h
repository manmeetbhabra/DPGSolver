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
 *  \brief Undefine macro definitions for c-style templated containers/functions relating to solution for the
 *         navier-stokes equations.
 */

#include "taylor_couette/undef_templates_solution_taylor_couette.h"

#undef compute_mu_fptr_T

#undef set_function_pointers_solution_navier_stokes_T
#undef convert_variables_gradients_T
#undef get_compute_mu_fptr_T
#undef set_viscosity_type_T
#undef compute_mu_constant_T
#undef compute_mu_sutherland_T

#undef set_function_pointers_num_flux
