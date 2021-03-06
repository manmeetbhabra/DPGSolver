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
 *  \brief Provides the macro definitions used for c-style templating related to exact normals.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function pointers
#define correct_for_exact_normal_fptr_T correct_for_exact_normal_fptr_d
///\}

///\{ \name Function names
#define set_correct_for_exact_normal_fptr_T      set_correct_for_exact_normal_fptr_d
#define correct_for_exact_normal_cylinder_T      correct_for_exact_normal_cylinder_d
#define correct_for_exact_normal_gaussian_bump_T correct_for_exact_normal_gaussian_bump_d
///\}

///\{ \name Static names
#define correct_normals_fc correct_normals_fc
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers
#define correct_for_exact_normal_fptr_T correct_for_exact_normal_fptr_c
///\}

///\{ \name Function names
#define set_correct_for_exact_normal_fptr_T      set_correct_for_exact_normal_fptr_c
#define correct_for_exact_normal_cylinder_T      correct_for_exact_normal_cylinder_c
#define correct_for_exact_normal_gaussian_bump_T correct_for_exact_normal_gaussian_bump_c
///\}

///\{ \name Static names
#define correct_normals_fc correct_normals_fc_c
///\}

#endif
