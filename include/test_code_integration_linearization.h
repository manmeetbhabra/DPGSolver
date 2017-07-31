// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_linearization_h__INCLUDED
#define DPG__test_code_integration_linearization_h__INCLUDED

#include <stdbool.h>

#include "petscmat.h"

struct S_linearization {
	bool         PrintEnabled, PrintTimings,
	             CheckFullLinearization, CheckWeakGradients,
	             update_argv, omit_root;
	char         **argvNew, *PrintName;
	unsigned int Nref, PGlobal, ML, PG_add, IntOrder_mult, IntOrder_add, check_level;
	int          nargc;

	Mat A, A_cs, A_csc;
	Vec b, b_cs, b_csc, x, x_cs, x_csc;
};

extern void test_linearization (struct S_linearization *const data, char const *const TestName);

#endif // DPG__test_code_integration_linearization_h__INCLUDED
