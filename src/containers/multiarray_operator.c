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

// Static function declarations ************************************************************************************* //

/// `mutable` version of \ref Operator
struct mutable_Operator {
	struct Matrix_d* op_std;            ///< The standard dense matrix operator.
	struct Multiarray_Matrix_d* ops_tp; ///< The multiarray of tensor-product sub-operators.
	struct Matrix_CSR_d* op_csr;        ///< The sparse matrix operator in CSR format.
};

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
		const struct Operator* op = a->data[i];

		printf("%-35s","\tdense operator:");
		if (op->op_std) {
			printf("\n\n");
/// \todo add '_tol' to names of printing functions with tolerance.
			print_const_Matrix_d(op->op_std,tol);
		} else {
			printf("*** NULL ***\n");
		}

		printf("%-35s","\ttensor-product sub-operators:");
		if (op->ops_tp) {
			printf("\n{\n\n");
			print_const_Multiarray_Matrix_d(op->ops_tp,tol);
			printf("}\n");
		} else {
			printf("*** NULL ***\n");
		}

		printf("%-35s","\tsparse (CSR) operator:");
		if (op->op_csr) {
			printf("\n\n");
			EXIT_ADD_SUPPORT;
		} else {
			printf("*** NULL ***\n");
		}
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
	assert(a != NULL);
	assert(a->data != NULL);

	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (int i = 0; i < size; ++i) {
		struct mutable_Operator* op = a->data[i];
		if (a->owns_data) {
			if (op->op_std)
				destructor_Matrix_d(op->op_std);
			if (op->ops_tp) {
				assert(op->ops_tp->owns_data == false);
				destructor_Multiarray_Matrix_d(op->ops_tp);
			}
			if (op->op_csr)
				EXIT_ADD_SUPPORT;
//				destructor_Matrix_CSR_d(op->op_csr);
		}
		free(op);
	}
	free(a->data);
	free(a->extents);
	free(a);
}