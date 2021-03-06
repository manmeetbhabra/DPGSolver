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
 *         navier-stokes equations (test case: taylor-couette).
 */

#undef set_sol_taylor_couette_T
#undef set_grad_taylor_couette_T
#undef constructor_const_sol_taylor_couette_T
#undef constructor_const_grad_taylor_couette_T

#undef constructor_sol_taylor_couette
#undef constructor_grad_taylor_couette
#undef Sol_Data__tc
#undef get_sol_data
#undef read_data_taylor_couette
#undef set_data_taylor_couette
