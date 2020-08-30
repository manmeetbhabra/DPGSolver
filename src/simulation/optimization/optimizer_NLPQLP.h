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

#ifndef DPG__optimizer_NLPQLP_h__INCLUDED
#define DPG__optimizer_NLPQLP_h__INCLUDED

#include "definitions_alloc.h"
struct Optimization_Case;


/** \brief Container for information relating to the NLPQLP optimizer. For more details on 
 * each variable, consult the NLPQLP documentation.
 */
struct Optimizer_NLPQLP_Data {

	// Optimization Parameters
	int NP, ///< Number of processors
		N, ///< Number of design variable dofs
		NMAX, ///< Must be greater than N
		M, ///< Total number of constraints
		ME, ///< Number of equality constraints
		IFAIL, ///< Optimization progress flag. Consult NLPQLP Documentation
		MODE, ///< Consult NLPQLP Documentation
		MNN2; ///< Equal to M+N+N+2

	int IOUT, ///< Consult NLPQLP Documentation
		MAXIT, ///< Consult NLPQLP Documentation
		MAXFUN, ///< Consult NLPQLP Documentation
		MAXNM, ///< Consult NLPQLP Documentation
		LQL, ///< Consult NLPQLP Documentation
		IPRINT; ///< Consult NLPQLP Documentation

	double 	ACC, ///< Consult NLPQLP Documentation
			ACCQP, ///< Consult NLPQLP Documentation
			STPMIN, ///< Consult NLPQLP Documentation
			RHO; ///< Consult NLPQLP Documentation

	// Design Points and their limits
	struct Multiarray_d *X, ///< Design point values
						*XL, ///< Design point lower limit
						*XU; ///< Design point upper limits

	// Objective and Constraint function values
	double F; ///< Objective function value
	struct Multiarray_d *G, ///< Constraint function values
						*dF, ///< Objective function gradient
						*dG; ///< Constraint function gradient

	// NLPQLP Structures
	struct Multiarray_d *C,	///< Consult NLPQLP Documentation
						*U, ///< Consult NLPQLP Documentation
						*D,	///< Consult NLPQLP Documentation
						*WA; ///< Consult NLPQLP Documentation

	struct Multiarray_i *KWA, ///< Consult NLPQLP Documentation
						*ACTIVE; ///< Consult NLPQLP Documentation

	int LWA, ///< Consult NLPQLP Documentation
		LKWA, ///< Consult NLPQLP Documentation
		LACTIV; ///< Consult NLPQLP Documentation
};


/**	\brief Use a line search method (gradient descent or BFGS) to find the optimal shape
 */
void optimizer_NLPQLP(
	struct Optimization_Case* optimization_case ///< Standard. Consult optimization_case.h
	);


#endif // DPG__optimizer_NLPQLP_h__INCLUDED