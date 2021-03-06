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
 *  \brief Provides templated functions relating to the linear Diffusion solutions.
 */

#include "def_templates_solution_diffusion.h"
#include "def_templates_test_case.h"

struct Test_Case_T;
struct Simulation;

/// \brief Set the solution function pointer members of an Diffusion \ref Test_Case_T.
void set_function_pointers_solution_diffusion_T
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

#include "undef_templates_solution_diffusion.h"
#include "undef_templates_test_case.h"
