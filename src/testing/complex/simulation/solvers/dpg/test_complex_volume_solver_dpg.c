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

#include "test_complex_volume_solver_dpg.h"

#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "volume_solver_dpg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //