// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_unit.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "test_base.h"
#include "test_unit_containers.h"
#include "test_unit_cubature.h"
#include "test_unit_bases.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_unit (struct Test_Info*const test_info)
{
	printf("\n\nRunning Unit Tests:\n");
	printf("-------------------------------------------------------------------------------------------------\n\n");

	test_unit_containers(test_info);
	test_unit_cubature(test_info);
	test_print_warning(test_info,"Need to add tests for plotting nodes");
	test_unit_bases(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
