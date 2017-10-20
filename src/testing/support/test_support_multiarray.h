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

#ifndef DPG__test_support_multiarray_h__INCLUDED
#define DPG__test_support_multiarray_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the containers defined in \ref multiarray.h.
 */

#include <stdio.h>
#include <stdbool.h>

struct Multiarray_Matrix_d;
struct const_Multiarray_d;
struct const_Multiarray_Matrix_d;

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a \ref Multiarray_d\* from data in the input file of the given name.
 *  \return Standard. */
struct Multiarray_d* constructor_file_name_Multiarray_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief `const` version of \ref constructor_file_name_Multiarray_d.
 *  \return Standard. */
const struct const_Multiarray_d* constructor_file_name_const_Multiarray_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Multiarray_Vector_i\* from data in the input file of the given name.
 *  \return Standard. */
struct Multiarray_Vector_i* constructor_file_name_Multiarray_Vector_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief `const` version of \ref constructor_file_name_Multiarray_Vector_i.
 *  \return Standard. */
const struct const_Multiarray_Vector_i* constructor_file_name_const_Multiarray_Vector_i
	(const char*const var_name,      ///< Defined for \ref constructor_file_name_Multiarray_Vector_i.
	 const char*const file_name_full ///< Defined for \ref constructor_file_name_Multiarray_Vector_i.
	);

/** \brief Constructor for a \ref const_Multiarray_Matrix_d\* from data in the input file of the given name.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_file_name_const_Multiarray_Matrix_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Multiarray_d\* from the current line in the input file.
 *  \return Standard. */
struct Multiarray_d* constructor_file_Multiarray_d
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief Constructor for a \ref Multiarray_Vector_i\* as read from a file.
 *  \return Standard. */
struct Multiarray_Vector_i* constructor_file_Multiarray_Vector_i
	(FILE* data_file ///< The file containing the data.
	);

// Difference functions ********************************************************************************************* //

/** \brief Check the difference between entries in the input \ref Multiarray_Vector_i\*s.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a, ///< Input 0.
	 const struct Multiarray_Vector_i*const b  ///< Input 1.
	);

/** \brief `const` version of \ref diff_Multiarray_Vector_i.
 *  \return See brief. */
bool diff_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const a, ///< Input 0.
	 const struct const_Multiarray_Vector_i*const b  ///< Input 1.
	);

/** \brief Check the relative difference between entries in the input \ref Multiarray_d\*s up to the input tolerance.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Multiarray_d
	(const struct Multiarray_d*const a, ///< Input 0.
	 const struct Multiarray_d*const b, ///< Input 1.
	 const double tol                   ///< The tolerance.
	);

/** \brief `const` version of \ref diff_Multiarray_d.
 *  \return See brief. */
bool diff_const_Multiarray_d
	(const struct const_Multiarray_d*const a, ///< Defined for \ref diff_Multiarray_d.
	 const struct const_Multiarray_d*const b, ///< Defined for \ref diff_Multiarray_d.
	 const double tol                         ///< Defined for \ref diff_Multiarray_d.
	);

/** \brief Check the relative difference between entries in the input \ref Multiarray_Matrix_d\*s up to the input
 *         tolerance.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Multiarray_Matrix_d
	(const struct Multiarray_Matrix_d*const a, ///< Input 0.
	 const struct Multiarray_Matrix_d*const b, ///< Input 1.
	 const double tol                          ///< The tolerance.
	);

/** \brief `const` version of \ref diff_Multiarray_Matrix_d.
 *  \return See brief. */
bool diff_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const a, ///< Defined for \ref diff_Multiarray_Matrix_d.
	 const struct const_Multiarray_Matrix_d*const b, ///< Defined for \ref diff_Multiarray_Matrix_d.
	 const double tol                                ///< Defined for \ref diff_Multiarray_Matrix_d.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print the difference of the input \ref Multiarray_Vector_i\*s.
void print_diff_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a, ///< Input 0.
	 const struct Multiarray_Vector_i*const b  ///< Input 1.
	);

/// \brief `const` version of \ref print_diff_Multiarray_Vector_i.
void print_diff_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const a, ///< Defined for \ref print_diff_Multiarray_Vector_i.
	 const struct const_Multiarray_Vector_i*const b  ///< Defined for \ref print_diff_Multiarray_Vector_i.
	);

/// \brief Print the relative difference of the input \ref Multiarray_d\*s, outputting 0 if less than the tolerance.
void print_diff_Multiarray_d
	(const struct Multiarray_d*const a, ///< Input 0.
	 const struct Multiarray_d*const b, ///< Input 1.
	 const double tol                   ///< The tolerance.
	);

/// \brief `const` version of \ref print_diff_Multiarray_d.
void print_diff_const_Multiarray_d
	(const struct const_Multiarray_d*const a, ///< Input 0.
	 const struct const_Multiarray_d*const b, ///< Input 1.
	 const double tol                         ///< The tolerance.
	);

/** \brief Print the relative difference of the input \ref Multiarray_Matrix_d\*s, outputting 0 if less than the
 *         tolerance. */
void print_diff_Multiarray_Matrix_d
	(const struct Multiarray_Matrix_d*const a, ///< Input 0.
	 const struct Multiarray_Matrix_d*const b, ///< Input 1.
	 const double tol                          ///< The tolerance.
	);

/// \brief `const` version of \ref print_diff_Multiarray_Matrix_d.
void print_diff_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const a, ///< Defined for \ref print_diff_Multiarray_Matrix_d.
	 const struct const_Multiarray_Matrix_d*const b, ///< Defined for \ref print_diff_Multiarray_Matrix_d.
	 const double tol                                ///< Defined for \ref print_diff_Multiarray_Matrix_d.
	);

#endif // DPG__test_support_multiarray_h__INCLUDED
