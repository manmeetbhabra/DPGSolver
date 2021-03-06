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

#ifndef DPG__test_support_computational_elements_h__INCLUDED
#define DPG__test_support_computational_elements_h__INCLUDED
/** \file
 *  \brief Provides supporting functions relating to the computational element processing.
 */

#include <stddef.h>

struct Intrusive_Link;

/** \brief Constructor for a \ref Intrusive_Link\* to be inserted in a list copied from the base.
 *  \return See brief. */
struct Intrusive_Link* constructor_copied_Intrusive_Link
	(struct Intrusive_Link* base, ///< Pointer to the base link.
	 const size_t sizeof_base,    ///< Value of std::sizeof(base).
	 const size_t sizeof_derived  ///< Value of std::sizeof(derived).
	);

#endif // DPG__test_support_computational_elements_h__INCLUDED
