// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include <stdlib.h>
#include <stdio.h>

#include "test_info.h"
#include "Macros.h"

#include "test_integration_euler.h"
#include "test_integration_linearization.h"
//#include "test_support.h"

void test_integration_euler (struct Test_Info*const test_info)
{
	if (test_info->t_int.equivalence_real_complex) {
		;
	} else {
//		test_print_warning("Euler equivalence real/complex testing currently disabled");
	}

	if (test_info->t_int.equivalence_algorithms) {
		;
	} else {
//		test_print_warning("Euler equivalence algorithms testing currently disabled");
	}

	if (test_info->t_int.linearization) {
#if 1
	#if 1
		#if 0
test_integration_linearization("test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
		#else
test_integration_linearization("test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
		#endif
	#else
		test_integration_linearization("test/Euler/Test_Euler_PeriodicVortex_Stationary_QUAD");
	#endif
#endif
	} else {
//		test_print_warning("Euler linearization testing currently disabled");
	}

	if (test_info->t_int.conv_order) {
		;
	} else {
//		test_print_warning("Euler convergence order testing currently disabled");
	}
}
