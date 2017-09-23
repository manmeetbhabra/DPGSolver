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

#ifndef DPG__mesh_connectivity_h__INCLUDED
#define DPG__mesh_connectivity_h__INCLUDED
/**	\file
 *	\brief Provides the interface to mesh connectivity containers and functions.
 */

struct const_Intrusive_List;
struct Mesh_Data;

/// \brief The number of 'M'aster and 'S'lave entities (always two but defined to avoid use of the magic number).
#define N_MS 2

/// \brief Container for local connectivity related information.
struct Conn_info {
	// Available from mesh_data:
	const int d;                                    ///< The dimension.
	const struct const_Vector_i*const elem_per_dim; ///< The number of elements of each dimension.

	struct const_Vector_i* volume_types; ///< Pointer to the first volume entry in \ref Mesh_Data::elem_types.
	struct Vector_i* v_n_lf;             ///< The number of local faces for each volume.

	// Computed here:
	struct Multiarray_Vector_i* f_ve;     ///< Global face to vertex correspondence.
	struct Vector_i*            ind_f_ve; ///< Indices of \ref f_ve after sorting.
};

/// \brief Holds data relating to the mesh connectivity.
struct Mesh_Connectivity {
	const struct const_Multiarray_Vector_i*const v_to_v; ///< Volume to volume connectivity.

	/** Volume to local face connectivity.
	 *	Redundant self-reference entries are replaced with the number corresponding to the appropriate boundary
	 *	condition. */
	const struct const_Multiarray_Vector_i*const v_to_lf;
};

/** \brief Constructor for the \ref Mesh_Connectivity.
 *	\return Standard. */
struct Mesh_Connectivity* constructor_Mesh_Connectivity
	(const struct Mesh_Data*const mesh_data,     ///< \ref Mesh_Data.
	 const struct const_Intrusive_List* elements ///< The base \ref Element list.
	);

/// \brief Destructor for \ref Mesh_Connectivity.
void destructor_Mesh_Connectivity
	(struct Mesh_Connectivity* mesh_conn ///< Standard.
	);

/// \brief Set the face node numbers and sort them.
void set_f_node_nums
	(struct Vector_i**const f_node_nums,         ///< The face node numbers.
	 const struct const_Vector_i*const node_nums ///< The entry of \ref Mesh_Data::node_nums for the current face.
	);

#endif // DPG__mesh_connectivity_h__INCLUDED