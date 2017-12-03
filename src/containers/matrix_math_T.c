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

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

Type compute_norm_Matrix_T_row
	(const ptrdiff_t row, const struct Matrix_T*const a, const char*const norm_type)
{
	const Type*const data = get_row_Matrix_T(row,a);
	const ptrdiff_t i_max   = a->ext_1;

	Type norm = 0.0;
	if (strstr(norm_type,"L2")) {
		for (ptrdiff_t i = 0; i < i_max; ++i)
			norm += data[i]*data[i];
		return sqrt(norm);
	}
	EXIT_UNSUPPORTED;
}

void transpose_Matrix_T (struct Matrix_T* a, const bool mem_only)
{
	if (a->layout == 'R')
		mkl_dimatcopy(a->layout,'T',(size_t)a->ext_0,(size_t)a->ext_1,1.0,a->data,(size_t)a->ext_1,(size_t)a->ext_0);
	else if (a->layout == 'C')
		mkl_dimatcopy(a->layout,'T',(size_t)a->ext_0,(size_t)a->ext_1,1.0,a->data,(size_t)a->ext_0,(size_t)a->ext_1);

	if (mem_only) {
		swap_layout(&a->layout);
	} else {
		ptrdiff_t tmp = a->ext_0;
		a->ext_0 = a->ext_1;
		a->ext_1 = tmp;
	}
}

void invert_sub_block_Matrix_T (struct Matrix_T* a, const ptrdiff_t row0, const ptrdiff_t col0, const ptrdiff_t ext)
{
	struct Matrix_T* a_sub = constructor_sub_block_Matrix_T(row0,col0,ext,ext,a); // destructed

	struct Matrix_T* a_sub_inv = constructor_inverse_Matrix_T(a_sub); // destructed
	destructor_Matrix_T(a_sub);

	set_block_Matrix_T(a,(struct const_Matrix_T*)a_sub_inv,row0,col0,'i');
	destructor_Matrix_T(a_sub_inv);
}

void scale_Matrix_T (struct Matrix_T* a, const Type val)
{
	ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] *= val;
}

void permute_Matrix_T (struct Matrix_T* a, const ptrdiff_t* p)
{
	assert((a->layout == 'R') || a->layout == 'C');
	assert(p != NULL);

	// There was a problem with the gsl_permute_matrix function call when setting tda outside of the gsl_matrix
	// definition. Not sure what the problem was but this is the reason for the redundant code.
	/// \todo Attempt to do this again and update comments.

	// It seems like gsl_permute_matrix naturally permutes the columns of a row-major stored matrix from the right.
	// Thus, transposition is required for permutation of row-major from the left/column-major from the right.
	if (a->layout == 'R') {
		const int n_p = (int)a->ext_0;

		size_t perm_st[n_p];
		for (int i = 0; i < n_p; ++i)
			perm_st[i] = (size_t)p[i];

		const gsl_permutation perm = { .size = (size_t)n_p, .data = perm_st, };

		transpose_Matrix_T(a,false);
		gsl_matrix A =
			{ .size1 = (size_t)a->ext_0,
			  .size2 = (size_t)a->ext_1,
			  .tda   = (size_t)a->ext_1,
			  .data  = a->data,
			  .block = NULL,
			  .owner = 0, };

		gsl_permute_matrix(&perm,&A);
		transpose_Matrix_T(a,false);
	} else {
		const int n_p = (int)a->ext_1;

		size_t perm_st[n_p];
		for (int i = 0; i < n_p; ++i)
			perm_st[i] = (size_t)p[i];

		const gsl_permutation perm = { .size = (size_t)n_p, .data = perm_st, };

		transpose_Matrix_T(a,true);
		gsl_matrix A =
			{ .size1 = (size_t)a->ext_0,
			  .size2 = (size_t)a->ext_1,
			  .tda   = (size_t)a->ext_1, /// \todo Ensure that this is correct.
			  .data  = a->data,
			  .block = NULL,
			  .owner = 0, };

		gsl_permute_matrix(&perm,&A);
		transpose_Matrix_T(a,true);
	}
}

void permute_Matrix_T_V (struct Matrix_T* a, const struct const_Vector_i* p_V)
{
	const ptrdiff_t ext_0 = p_V->ext_0;
	assert((a->layout == 'R' && a->ext_0 == ext_0) || (a->layout == 'C' && a->ext_1 == ext_0));

	ptrdiff_t p[ext_0];
	for (int i = 0; i < ext_0; ++i)
		p[i] = p_V->data[i];

	permute_Matrix_T(a,p);
}

void mm_T
	(const char trans_a_i, const char trans_b_i, const Real alpha, const Real beta,
	 const struct const_Matrix_T*const a, const struct const_Matrix_T*const b, struct Matrix_T*const c)
{
	const CBLAS_LAYOUT    layout = ( c->layout == 'R' ? CBRM : CBCM );
	const CBLAS_TRANSPOSE transa = ( (c->layout == a->layout) == (trans_a_i == 'N') ? CBNT : CBT ),
	                      transb = ( (c->layout == b->layout) == (trans_b_i == 'N') ? CBNT : CBT );
	const MKL_INT m   = (MKL_INT) c->ext_0,
	              n   = (MKL_INT) c->ext_1,
	              k   = (MKL_INT) ( trans_a_i == 'N' ? a->ext_1 : a->ext_0 ),
	              lda = (MKL_INT) ( a->layout == 'R' ? a->ext_1 : a->ext_0 ),
	              ldb = (MKL_INT) ( b->layout == 'R' ? b->ext_1 : b->ext_0 ),
	              ldc = (MKL_INT) ( c->layout == 'R' ? c->ext_1 : c->ext_0 );

	assert(m == ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ));
	assert(n == ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 ));
	assert(k == ( trans_b_i == 'N' ? b->ext_0 : b->ext_1 ));

	cblas_dgemm(layout,transa,transb,m,n,k,alpha,a->data,lda,b->data,ldb,beta,c->data,ldc);
}

void mv_T
	(const char trans_a_i, const Real alpha, const Real beta,
	 const struct const_Matrix_T*const a, const struct const_Vector_d*const b, struct Vector_d*const c)
{
	const CBLAS_LAYOUT    layout = ( a->layout == 'R' ? CBRM : CBCM );
	const CBLAS_TRANSPOSE transa = ( trans_a_i == 'N' ? CBNT : CBT );

	/// \note Unlike the \ref mm_d function, m and n here represent the dimensions of A and not op(A).
	const MKL_INT m   = (MKL_INT) a->ext_0,
	              n   = (MKL_INT) a->ext_1,
	              lda = (MKL_INT) ( a->layout == 'R' ? a->ext_1 : a->ext_0 ),
	              ldb = 1,
	              ldc = 1;

	assert(m > 0);
	assert(n > 0);
	assert(m == ( transa == CBNT ? c->ext_0 : b->ext_0 ));
	assert(n == ( transa == CBNT ? b->ext_0 : c->ext_0 ));

	cblas_dgemv(layout,transa,m,n,alpha,a->data,lda,b->data,ldb,beta,c->data,ldc);
}

void scale_Matrix_by_Vector_T
	(const char side, const Real alpha, struct Matrix_T*const a, const struct const_Vector_d*const b,
	 const bool invert_diag)
{
	if (invert_diag)
		invert_Vector_d((struct Vector_d*)b);

	if (alpha != 1.0)
		scale_Matrix_T(a,alpha);

	const ptrdiff_t n_row = a->ext_0,
	                n_col = a->ext_1;

	bool transpose_a = false;
	if (side == 'L') {
		assert(b->ext_0 == a->ext_0);

		if (a->layout == 'C') {
			transpose_a = true;
			transpose_Matrix_T(a,true);
		}

		for (ptrdiff_t row = 0; row < n_row; ++row) {
			const Type val = b->data[row];
			Type* data_row = get_row_Matrix_T(row,a);
			for (ptrdiff_t col = 0; col < n_col; ++col)
				*data_row++ *= val;
		}
	} else if (side == 'R') {
		assert(b->ext_0 == a->ext_1);

		if (a->layout == 'R') {
			transpose_a = true;
			transpose_Matrix_T(a,true);
		}

		for (ptrdiff_t col = 0; col < n_col; ++col) {
			const Type val = b->data[col];
			Type* data_col = get_col_Matrix_T(col,a);
			for (ptrdiff_t row = 0; row < n_row; ++row)
				*data_col++ *= val;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
	if (transpose_a)
		transpose_Matrix_T(a,true);
	if (invert_diag)
		invert_Vector_d((struct Vector_d*)b);
}

void mm_diag_T
	(const char side, const Real alpha, const Real beta, const struct const_Matrix_T*const a,
	 const struct const_Vector_d*const b, struct Matrix_T* c, const bool invert_diag)
{
	assert(a->ext_0 == c->ext_0);
	assert(a->ext_1 == c->ext_1);
	assert(a->layout == c->layout);

	if (invert_diag)
		invert_Vector_d((struct Vector_d*)b);

	if (beta != 1.0)
		scale_Matrix_T(c,beta);

	const ptrdiff_t n_row = a->ext_0,
	                n_col = a->ext_1;

	if (side == 'L') {
		assert(b->ext_0 == a->ext_0);

		if (a->layout == 'R') {
			for (ptrdiff_t row = 0; row < n_row; ++row) {
				const Type val = b->data[row];
				const Type* data_a = get_row_const_Matrix_T(row,a);
				Type* data_c       = get_row_Matrix_T(row,c);
				for (ptrdiff_t col = 0; col < n_col; ++col)
					data_c[col] += alpha*data_a[col]*val;
			}
		} else if (a->layout == 'C') {
			for (ptrdiff_t col = 0; col < n_col; ++col) {
				const Type* data_a = get_col_const_Matrix_T(col,a);
				Type* data_c       = get_col_Matrix_T(col,c);
				for (ptrdiff_t row = 0; row < n_row; ++row) {
					const Type val = b->data[row];
					data_c[col] += alpha*data_a[col]*val;
				}
			}
		}
	} else if (side == 'R') {
		assert(b->ext_0 == a->ext_1);

		if (a->layout == 'R') {
			for (ptrdiff_t row = 0; row < n_row; ++row) {
				const Type* data_a = get_row_const_Matrix_T(row,a);
				Type* data_c       = get_row_Matrix_T(row,c);
				for (ptrdiff_t col = 0; col < n_col; ++col) {
					const Type val = b->data[col];
					data_c[col] += alpha*data_a[col]*val;
				}
			}
		} else if (a->layout == 'C') {
			for (ptrdiff_t col = 0; col < n_col; ++col) {
				const Type val = b->data[col];
				const Type* data_a = get_col_const_Matrix_T(col,a);
				Type* data_c       = get_col_Matrix_T(col,c);
				for (ptrdiff_t row = 0; row < n_row; ++row)
					data_c[row] += alpha*data_a[row]*val;
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}

	if (invert_diag)
		invert_Vector_d((struct Vector_d*)b);
}

void reinterpret_const_Matrix_T (const struct const_Matrix_T* a, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	assert(ext_0*ext_1 == ((a->ext_0)*(a->ext_1)));
	const_cast_ptrdiff(&a->ext_0,ext_0);
	const_cast_ptrdiff(&a->ext_1,ext_1);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //