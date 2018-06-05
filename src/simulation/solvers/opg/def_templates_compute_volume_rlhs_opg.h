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
 *  \brief Provides the macro definitions used for c-style templating related to the rlhs computing functions for the opg
 *         volumes.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define S_Params_T S_Params
///\}

///\{ \name Function names
#define compute_volume_rlhs_opg_T compute_volume_rlhs_opg
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define S_Params_T S_Params_c
///\}

///\{ \name Function names
#define compute_volume_rlhs_opg_T compute_volume_rlhs_opg_c
///\}

#endif
