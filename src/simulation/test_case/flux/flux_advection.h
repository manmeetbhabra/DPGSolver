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

#ifndef DPG__flux_advection_h__INCLUDED
#define DPG__flux_advection_h__INCLUDED
/** \file
 *  \brief Provides functions relating to linear advection fluxes.
 */

struct Flux_Input;
struct mutable_Flux;

/// \brief Compute the fluxes for the linear advection equation.
void compute_Flux_advection
	(const struct Flux_Input* flux_i, ///< \ref Flux_Input.
	 struct mutable_Flux* flux        ///< \ref Flux.
	);

/// \brief Compute the fluxes (and optionally the flux Jacobians) for the linear advection equation.
void compute_Flux_advection_jacobian
	(const struct Flux_Input* flux_i, ///< \ref Flux_Input.
	 struct mutable_Flux* flux        ///< \ref Flux.
	);

#endif // DPG__flux_advection_h__INCLUDED