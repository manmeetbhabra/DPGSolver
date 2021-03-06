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
 *  \brief Provides the macro definitions used for c-style templating related to the matrix print functions.
 */

#if defined TYPE_RC
#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define check_Matrix_extents_zero_T check_Matrix_extents_zero_d

#define print_Matrix_T_tol           print_Matrix_d_tol
#define print_const_Matrix_T_tol     print_const_Matrix_d_tol
#define print_Matrix_T               print_Matrix_d
#define print_const_Matrix_T         print_const_Matrix_d
#define print_to_file_Matrix_T       print_to_file_Matrix_d
#define print_to_file_const_Matrix_T print_to_file_const_Matrix_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define check_Matrix_extents_zero_T check_Matrix_extents_zero_c

#define print_Matrix_T_tol           print_Matrix_c_tol
#define print_const_Matrix_T_tol     print_const_Matrix_c_tol
#define print_Matrix_T               print_Matrix_c
#define print_const_Matrix_T         print_const_Matrix_c
#define print_to_file_Matrix_T       print_to_file_Matrix_c
#define print_to_file_const_Matrix_T print_to_file_const_Matrix_c
///\}

#endif

///\{ \name Real Data types/Function names
#define print_Matrix_R           print_Matrix_d
#define print_const_Matrix_R     print_const_Matrix_d
///\}


#elif defined TYPE_I
#if TYPE_I == TYPE_II

///\{ \name Function names
#define check_Matrix_extents_zero_T check_Matrix_extents_zero_i

#define print_Matrix_T_tol           print_Matrix_i_tol
#define print_const_Matrix_T_tol     print_const_Matrix_i_tol
#define print_Matrix_T               print_Matrix_i
#define print_const_Matrix_T         print_const_Matrix_i
#define print_to_file_Matrix_T       print_to_file_Matrix_i
#define print_to_file_const_Matrix_T print_to_file_const_Matrix_i
///\}

#endif

#endif
