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

#include "multiarray_operator.h"

#include <assert.h>

#include "macros.h"
#include "definitions_tol.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "operator.h"

// Static function declarations ************************************************************************************* //

/// `mutable` version of \ref Multiarray_Operator.
struct mutable_Multiarray_Operator {
	int order;          ///< Defined in \ref Multiarray_Operator.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_Operator.

	bool owns_data;                 ///< Defined in \ref Multiarray_Operator.
	struct mutable_Operator** data; ///< Defined in \ref Multiarray_Operator.
};

/** \brief Move constructor for a \ref Multiarray_Operator\* with the input extents having been previously
 *         dynamically allocated.
 *  \return See brief. */
const struct Multiarray_Operator* constructor_move_Multiarray_Operator_dyn_extents
	(const int order,            ///< Standard.
	 ptrdiff_t*const extents,    ///< Standard.
	 const bool owns_data,       ///< Standard.
	 struct Operator**const data ///< Standard.
	);

/// \brief `mutable` version of \ref destructor_Multiarray_Operator.
static void destructor_mutable_Multiarray_Operator
	(struct mutable_Multiarray_Operator* a ///< Standard.
	);

/// \brief `mutable` version of \ref set_MO_from_MO.
void set_mutable_MO_from_MO
	(struct mutable_Multiarray_Operator* dest, ///< The destination.
	 struct mutable_Multiarray_Operator* src,  ///< The source.
	 const int order_o,                        ///< The order of the output (destination).
	 const ptrdiff_t*const sub_indices         ///< The sub-indices specifying which part of the source to extract.
	);

// Constructor functions ******************************************************************************************** //
// Default constructors ********************************************************************************************* //

const struct Multiarray_Operator* constructor_default_Multiarray_Operator ()
{
	return constructor_move_Multiarray_Operator_dyn_extents(0,NULL,true,NULL);
}

// Empty constructors *********************************************************************************************** //

const struct Multiarray_Operator* constructor_empty_Multiarray_Operator_V
	(const struct const_Vector_i*const extents_i_V)
{
	const ptrdiff_t order = extents_i_V->ext_0;

	ptrdiff_t extents_i[order];
	for (ptrdiff_t i = 0; i < order; ++i)
		extents_i[i] = extents_i_V->data[i];

	return constructor_empty_Multiarray_Operator(order,extents_i);
}

const struct Multiarray_Operator* constructor_empty_Multiarray_Operator
	(const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep

	const ptrdiff_t size = compute_size(order,extents);
	struct Operator** data = calloc(size , sizeof *data); // keep
	for (int i = 0; i < size; ++i)
		data[i] = calloc(1,sizeof *data[i]); // keep

	return constructor_move_Multiarray_Operator_dyn_extents(order,extents,true,data);
}

// Destructors ****************************************************************************************************** //

void destructor_Multiarray_Operator (const struct Multiarray_Operator* a)
{
	destructor_mutable_Multiarray_Operator((struct mutable_Multiarray_Operator*)a);
}

void destructor_Multiarray2_Operator (const struct Multiarray_Operator* a[2])
{
	destructor_Multiarray_Operator(a[0]);
	destructor_Multiarray_Operator(a[1]);
}

// Setter functions ************************************************************************************************* //

void set_MO_from_MO
	(const struct Multiarray_Operator* dest, const struct Multiarray_Operator* src, const int order_o,
	 const ptrdiff_t*const sub_indices)
{
	set_mutable_MO_from_MO(
		(struct mutable_Multiarray_Operator*)dest,(struct mutable_Multiarray_Operator*)src,order_o,sub_indices);
}

const struct Operator* get_Multiarray_Operator
	(const struct Multiarray_Operator* src, const ptrdiff_t*const sub_indices)
{
	assert(src != NULL);
	return src->data[compute_index_sub_container(src->order,0,src->extents,sub_indices)];
}

// Printing functions *********************************************************************************************** //

void print_Multiarray_Operator (const struct Multiarray_Operator*const a)
{
	print_Multiarray_Operator_tol(a,EPS);
}

void print_Multiarray_Operator_tol (const struct Multiarray_Operator*const a, const double tol)
{
	print_Multiarray_extents(a->order,a->extents);

	const ptrdiff_t size = compute_size(a->order,a->extents);

	for (ptrdiff_t i = 0; i < size; i++) {
		printf("\nIndex (MO) % 3td:\n",i);

		print_Operator_tol(a->data[i],tol);
	}
	printf("\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

const struct Multiarray_Operator* constructor_move_Multiarray_Operator_dyn_extents
	(const int order, ptrdiff_t*const extents, const bool owns_data, struct Operator**const data)
{
	struct mutable_Multiarray_Operator* dest = calloc(1,sizeof *dest); // returned

	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = (struct mutable_Operator**) data;

	return (struct Multiarray_Operator*) dest;
}

static void destructor_mutable_Multiarray_Operator (struct mutable_Multiarray_Operator* a)
{
	// Multiarray_Operators may be NULL for tensor-product operator sub-operator components of tensor-product
	// elements.
	if (a == NULL)
		return;

	assert(a->data != NULL);

	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (int i = 0; i < size; ++i) {
		struct mutable_Operator* op = a->data[i];
		if (a->owns_data)
			destructor_mutable_Operator(op);
		else
			free(op);
	}
	free(a->data);
	free(a->extents);
	free(a);
}

void set_mutable_MO_from_MO
	(struct mutable_Multiarray_Operator* dest, struct mutable_Multiarray_Operator* src, const int order_o,
	 const ptrdiff_t*const sub_indices)
{
	assert(src != NULL);
	assert(order_o == 1);

	dest->owns_data = false;
	dest->order     = order_o;
	dest->extents   = src->extents;
	dest->data      = &src->data[compute_index_sub_container(src->order,dest->order,src->extents,sub_indices)];
}