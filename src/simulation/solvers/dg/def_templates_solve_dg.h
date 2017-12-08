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
 *  \brief Provides the macro definitions used for c-style templating related to the dg solver functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define update_ind_dof_dg_T     update_ind_dof_dg
#define permute_Multiarray_T_fc permute_Multiarray_d_fc
#define constructor_nnz_dg_T    constructor_nnz_dg
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define update_ind_dof_dg_T     update_ind_dof_dg_c
#define permute_Multiarray_T_fc permute_Multiarray_c_fc
#define constructor_nnz_dg_T    constructor_nnz_dg_c
///\}

#endif
