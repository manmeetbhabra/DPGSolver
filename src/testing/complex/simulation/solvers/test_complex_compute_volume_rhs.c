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

#include "test_complex_compute_volume_rhs.h"

#include <assert.h>
#include <stdlib.h>

#include "macros.h"

#include "test_complex_compute_volume_rhs_dg.h"
#include "test_complex_flux.h"
#include "test_support_math_functions.h"

#include "complex_multiarray.h"
#include "multiarray.h"

#include "volume_solver.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "compute_volume_rlhs_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
