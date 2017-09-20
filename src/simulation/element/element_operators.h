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

#ifndef DPG__element_operators_h__INCLUDED
#define DPG__element_operators_h__INCLUDED
/** \file
 *  \brief Provides the functions relating to element operators.
 */

#include <stddef.h>
#include <stdbool.h>

struct Simulation;
struct const_Element;

/// Container for operator range related information.
struct Operator_Info {
	const struct const_Element* element; ///< \ref Element.

	const int range_d, ///< Range of dimensions (For differentiation operators).
	          range_f, ///< Range of faces.
	          range_h, ///< Range of h-refinement related operators.
	          range_p; ///< Range of orders.

	const int op_spec_sc

	const int cub_type; ///< The type of cubature.

	const int p_ref[2]; ///< Reference polynomial orders from \ref Simulation.

	struct Vector_i* extents_cub; ///< The extents of the associated \ref Multiarray_Cubature\*.
	struct Vector_i* extents_op;  ///< The extents of the associated \ref Multiarray_Matrix_d\* of operators.

	struct Matrix_i* values_op ///< The values of d, f, h, p_in, and p_out for each operator.
};

/// Container for a Multiarray of \ref Cubature\* data.
struct Multiarray_Cubature {
	int order;          ///< Defined in \ref Multiarray_d.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_d.

	bool owns_data;                    ///< Defined in \ref Multiarray_d.
	const struct const_Cubature** data; ///< Defined in \ref Multiarray_d.
};

/// `const` version of \ref Multiarray_Cubature.
struct const_Multiarray_Cubature {
	const int order;               ///< Defined in \ref Multiarray_d.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_d.

	const bool owns_data;                         ///< Defined in \ref Multiarray_d.
	const struct const_Cubature*const*const data; ///< Defined in \ref Multiarray_d.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for the \ref Operator_Info\* having the given inputs.
 *  \return Standard. */
struct Operator_Info* constructor_Operator_Info
	(const int cmp_rng,                  /**< Holds the compound range value. Supported values can be found in
	                                      *   \ref definitions_element_operators.h */
	 const int*const cub_type_info,      ///< Holds the Operator_Info:: \todo update comment.
	 const int p_ref[2],                 ///< Defined in \ref Operator_Info.
	 const struct const_Element* element ///< Defined in \ref Operator_Info.
	);

/// \brief Destructor for a \ref Operator_Info\*.
void destructor_Operator_Info
	(struct Operator_Info* op_ranges ///< Standard.
	);

/** \brief Constructor for a \ref const_Multiarray_Cubature\* holding the cubature nodes (and weights if applicable) for
 *         the range of supported hp adaptive operators.
 *  \return Standard.
 *
 *  The required nodes are determined from the last `n_hp` values in the `ext_v1_V` vector. Cubature nodes on
 *  element sub-regions are computed by:
 *  1. Computing the barycentric coordinates of the nodes on the complete reference element;
 *  2. Multiplying with the appropriate vertices of the element sub-regions.
 */
const struct const_Multiarray_Cubature* constructor_const_Multiarray_Cubature
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Operator_Info* op_info  ///< \ref Operator_Info.
	);

/// \brief Destructor for a \ref const_Multiarray_Cubature\* container.
void destructor_const_Multiarray_Cubature
	(const struct const_Multiarray_Cubature*const a ///< Standard.
	);

/** \brief Constructor for a \ref const_Multiarray_Matrix_d\* of operators.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_operators_Multiarray_Matrix_d_V
	(const int op_type,                   ///< Operator type. See options in \ref definitions_element_operators.h.
	 const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

#endif // DPG__element_operators_h__INCLUDED
