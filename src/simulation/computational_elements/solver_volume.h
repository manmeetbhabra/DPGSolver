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

#ifndef DPG__solver_volume_h__INCLUDED
#define DPG__solver_volume_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Solver_Volume container and associated functions.
 */

#include "volume.h"

/// \brief Container for data relating to the solver volumes.
struct Solver_Volume {
	struct Volume volume; ///< The base \ref Volume.

	/// The reference order of the volume. Need not be equal to the order of the solution in the volume.
	const int p_ref;

	/** The geometry coefficients of the volume in the \ref Simulation::basis_geom. For each of the supported
	 *  \ref Simulation::domain_type options, geom_coef represents:
	 *	- DOM_STRAIGHT: the projection of xyz_ve into the geometry basis of order 1.
	 *	- DOM_CURVED:
	 *		- straight volumes: [See DOM_STRAIGHT];
	 *		- boundary volumes: the coefficients of the *blended* face geometry of order k_g.
	 *	- DOM_PARAMETRIC: the coefficients of the mapped geometry of order k_g.
	 */
	const struct const_Multiarray_d*const geom_coef;

	/// The coefficients of the solution in the \ref Simulation::basis_sol.
	struct Multiarray_d* sol_coef;

	/// The coefficients of the solution gradient in the \ref Simulation::basis_sol.
	struct Multiarray_d* grad_coef;

	/** The metric terms (cofactors of the geometry Jacobian) used for transformation of integrals between physical
	 *  and computational space stored at the (v)olume (m)etric nodes. */
	const struct const_Multiarray_d*const metrics_vm;

	/** The metric terms (cofactors of the geometry Jacobian) used for transformation of integrals between physical
	 *  and computational space stored at the (v)olume (c)ubature nodes. */
	const struct const_Multiarray_d*const metrics_vc;

	/// The determinant of the geometry mapping Jacobian evaluated at the volume cubature nodes.
	const struct const_Multiarray_d*const jacobian_det_vc;
};

/// \brief Constructor for a derived \ref Solver_Volume.
void constructor_derived_Solver_Volume
	(struct Volume* volume_ptr,   ///< Pointer to the volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Solver_Volume.
void destructor_derived_Solver_Volume
	(struct Volume* volume_ptr ///< Pointer to the volume.
	);

#endif // DPG__solver_volume_h__INCLUDED
