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

#ifndef DPG__compute_error_euler_h__INCLUDED
#define DPG__compute_error_euler_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for error computation and output relating to the Euler variables.
 */

struct Simulation;
struct Error_CE_Data;

///\{ \name Parameters related to removal of lift/drag functionals for error computation.
#define REMOVE_CD_CL__NONE 0   ///< Remove neither.
#define REMOVE_CD_CL__D    101 ///< Remove drag.
///\}

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of all supported Euler variables.
 *  \return See brief.
 *
 *  This function should be used when the exact solution is available.
 */
struct Error_CE* constructor_Error_CE_euler_all
	(const struct Simulation* sim ///< See brief.
	);

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of all supported Euler variables and each of
 *         their residuals.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_euler_all_p_rhs
	(const struct Simulation* sim ///< See brief.
	);

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of the entropy.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_euler_entropy
	(const struct Simulation* sim ///< Defined for \ref constructor_Error_CE_fptr.
	);

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of the entropy and each of the residuals.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_euler_entropy_p_rhs
	(const struct Simulation* sim ///< Defined for \ref constructor_Error_CE_fptr.
	);

/// \brief Add the computed and exact specified variable to the solution data Multiarrays.
void add_euler_variable_Error_CE_Data
	(const char var_type,               ///< The variable type to add. Options: 's' (entropy), 't'emperature.
	 struct Error_CE_Data*const e_ce_d, ///< \ref Error_CE_Data.
	 const struct Simulation*const sim  ///< \ref Simulation.
	);

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of pressure drag and lift coefficients.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_functionals__cd_cl
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Version of \ref constructor_Error_CE_functionals__cd_cl with drag coefficient removed.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_functionals__cl
	(const struct Simulation*const sim ///< See brief.
	);

#endif // DPG__compute_error_euler_h__INCLUDED
