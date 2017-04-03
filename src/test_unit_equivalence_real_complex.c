// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_equivalence_real_complex.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "test_code_fluxes_inviscid.h"
#include "test_support.h"
#include "array_norm.h"
#include "initialize_test_case.h"
#include "fluxes_inviscid.h"
#include "fluxes_inviscid_c.h"
#include "boundary_conditions.h"
#include "boundary_conditions_c.h"
#include "variable_functions.h"
#include "variable_functions_c.h"

/*
 *	Purpose:
 *		Test equivalence of output between real and complex versions of the same functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static unsigned int compare_flux_inviscid (const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                           const unsigned int Neq, double *Wr);
static unsigned int compare_flux_Num      (const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                           const unsigned int Neq, double *Wr, double *nL, const char *nFType);
static unsigned int compare_boundary      (const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                           const unsigned int Neq, double *WLr, double *nL, double *XYZ,
                                           const char *TestCase);
static unsigned int compare_variables     (const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                           const unsigned int Neq, double *Wr);

void test_unit_equivalence_real_complex(void)
{
	unsigned int pass;

	char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free

	/*
	 *	Input:
	 *		Real functions input.
	 *
	 *	Expected Output:
	 *		No difference between outputs.
	 *
	 */

	printf("\nWarning: Ensure that tests for real/complex equivalence for all relevant functions are implemented.\n\n");
	TestDB.Nwarnings++;

	unsigned int Nn, Nel, d, Neq;
	double       *W, *nL, *XYZ;

	DB.TestCase     = calloc(STRLEN_MAX , sizeof *(DB.TestCase));     // free
	DB.PDE          = calloc(STRLEN_MAX , sizeof *(DB.PDE));          // free
	DB.PDESpecifier = calloc(STRLEN_MAX , sizeof *(DB.PDESpecifier)); // free
	DB.Geometry     = calloc(STRLEN_MAX , sizeof *(DB.Geometry));     // free
	DB.MeshType     = calloc(STRLEN_MAX , sizeof *(DB.MeshType));     // free

	strcpy(DB.TestCase,"SupersonicVortex");
	strcpy(DB.PDE,"Euler");
	strcpy(DB.PDESpecifier,"Internal");
	strcpy(DB.Geometry,"n-Cylinder_HollowSection");
	strcpy(DB.MeshType,"Curved");

	initialize_test_case_parameters();

	for (d = 1; d <= 3; d++) {
		Neq  = d+2;

		W   = initialize_W(&Nn,&Nel,d); // free
		nL  = initialize_n(Nn,Nel,d);   // free
		XYZ = initialize_XYZ(Nn,Nel,d); // free

		// flux_inviscid
		pass = compare_flux_inviscid(Nn,Nel,d,Neq,W);
		if (d == 1)
			sprintf(PrintName,"equivalence_flux_inviscid     (d = %d):",d);
		else
			sprintf(PrintName,"            flux_inviscid     (d = %d):",d);
		test_print2(pass,PrintName);

		// flux_LF
		pass = compare_flux_Num(Nn,Nel,d,Neq,W,nL,"LF");
		sprintf(PrintName,"            flux_LF                  :");
		test_print2(pass,PrintName);

		// flux_Roe
		pass = compare_flux_Num(Nn,Nel,d,Neq,W,nL,"Roe");
		sprintf(PrintName,"            flux_Roe                 :");
		test_print2(pass,PrintName);

		// boundary_SlipWall
		pass = compare_boundary(Nn,Nel,d,Neq,W,nL,XYZ,"SlipWall");
		sprintf(PrintName,"            boundary_SlipWall        :");
		test_print2(pass,PrintName);

		// boundary_Riemann
		if (d != 1) {
			pass = compare_boundary(Nn,Nel,d,Neq,W,nL,XYZ,"Riemann");
			sprintf(PrintName,"            boundary_Riemann         :");
			test_print2(pass,PrintName);
		}

		// convert_variables
		pass = compare_variables(Nn,Nel,d,Neq,W);
		sprintf(PrintName,"            convert_variables        :");
		test_print2(pass,PrintName);

		free(W);
		free(nL);
		free(XYZ);
	}

	free(DB.TestCase);
	free(DB.PDE);
	free(DB.PDESpecifier);
	free(DB.Geometry);
	free(DB.MeshType);
	free(DB.SolverType);

	free(PrintName);
}

static unsigned int compare_flux_inviscid(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                          const unsigned int Neq, double *Wr)
{
	unsigned int pass = 0;

	unsigned int   i, iMax, NnTotal, Nvar;
	double         *Fr, *Fctr;
	double complex *Wc, *Fc;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	Wc   = malloc(NnTotal*Nvar  * sizeof *Wc);   // free
	Fr   = malloc(NnTotal*d*Neq * sizeof *Fr);   // free
	Fc   = malloc(NnTotal*d*Neq * sizeof *Fc);   // free
	Fctr = malloc(NnTotal*d*Neq * sizeof *Fctr); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		Wc[i] = Wr[i];

	flux_inviscid(Nn,Nel,Wr,Fr,d,Neq);
	flux_inviscid_c(Nn,Nel,Wc,Fc,d,Neq);

	for (i = 0, iMax = NnTotal*d*Neq; i < iMax; i++)
		Fctr[i] = creal(Fc[i]);

	if (array_norm_diff_d(NnTotal*d*Neq,Fr,Fctr,"Inf") < EPS)
		pass = 1;

	free(Wc);
	free(Fr);
	free(Fc);
	free(Fctr);

	return pass;
}

static unsigned int compare_flux_Num(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                     const unsigned int Neq, double *Wr, double *nL, const char *nFType)
{
	unsigned int pass = 0;

	unsigned int   i, iMax, n, var, Nvar, IndWLR;
	double         *W_ptr, *WLr, *WRr, *nFr, *nFctr;
	double complex *WLc, *WRc, *nFc;

	if (Nel != 2)
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	Nvar = Neq;

	WLr   = malloc(Nn*Nvar * sizeof *WLr);   // free
	WRr   = malloc(Nn*Nvar * sizeof *WRr);   // free
	WLc   = malloc(Nn*Nvar * sizeof *WLc);   // free
	WRc   = malloc(Nn*Nvar * sizeof *WRc);   // free
	nFr   = malloc(Nn*Nvar * sizeof *nFr);   // free
	nFc   = malloc(Nn*Nvar * sizeof *nFc);   // free
	nFctr = malloc(Nn*Nvar * sizeof *nFctr); // free

	W_ptr = Wr;
	for (var = 0; var < Nvar; var++) {
		IndWLR = var*Nn;
		for (n = 0; n < Nn; n++)
			WLr[IndWLR+n] = *W_ptr++;
		for (n = 0; n < Nn; n++)
			WRr[IndWLR+n] = *W_ptr++;
	}

	for (i = 0, iMax = Nn*Nvar; i < iMax; i++) {
		WLc[i] = WLr[i];
		WRc[i] = WRr[i];
	}

	if (strstr(nFType,"LF")) {
		flux_LF(Nn,1,WLr,WRr,nFr,nL,d,Neq);
		flux_LF_c(Nn,1,WLc,WRc,nFc,nL,d,Neq);
	} else if (strstr(nFType,"Roe")) {
		flux_Roe(Nn,1,WLr,WRr,nFr,nL,d,Neq);
		flux_Roe_c(Nn,1,WLc,WRc,nFc,nL,d,Neq);
	} else {
		printf("Error: Unsupported nFType.\n"), EXIT_MSG;
	}

	for (i = 0, iMax = Nn*Neq; i < iMax; i++)
		nFctr[i] = creal(nFc[i]);

	if (array_norm_diff_d(Nn*Neq,nFr,nFctr,"Inf") < EPS)
		pass = 1;

	free(WLr);
	free(WRr);
	free(WLc);
	free(WRc);
	free(nFr);
	free(nFc);
	free(nFctr);

	return pass;
}

static unsigned int compare_boundary(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                     const unsigned int Neq, double *WLr, double *nL, double *XYZ, const char *TestCase)
{
	unsigned int pass = 0;

	unsigned int   i, iMax, NnTotal, Nvar;
	double         *WBr, *WBctr;
	double complex *WLc, *WBc;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	WBr   = malloc(NnTotal*Nvar * sizeof *WBr);   // free
	WLc   = malloc(NnTotal*Nvar * sizeof *WLc);   // free
	WBc   = malloc(NnTotal*Nvar * sizeof *WBc);   // free
	WBctr = malloc(NnTotal*Nvar * sizeof *WBctr); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WLc[i] = WLr[i];

	if (strstr(TestCase,"SlipWall")) {
		boundary_SlipWall(Nn,Nel,WLr,WBr,nL,d);
		boundary_SlipWall_c(Nn,Nel,WLc,WBc,nL,d);
	} else if (strstr(TestCase,"Riemann")) {
		boundary_Riemann(Nn,Nel,XYZ,WLr,NULL,WBr,nL,d);
		boundary_Riemann_c(Nn,Nel,XYZ,WLc,NULL,WBc,nL,d);
	}

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WBctr[i] = creal(WBc[i]);

	if (array_norm_diff_d(NnTotal*Nvar,WBr,WBctr,"Inf") < 1e1*EPS)
		pass = 1;

	free(WBr);
	free(WLc);
	free(WBc);
	free(WBctr);

	return pass;
}

static unsigned int compare_variables(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                      const unsigned int Neq, double *Wr)
{
	unsigned int pass = 0;

	unsigned int   i, iMax, NnTotal, Nvar;
	double         *Ur, *Wctr, *Uctr;
	double complex *Wc, *Uc;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	Ur   = malloc(NnTotal*Nvar * sizeof *Ur);   // free
	Wc   = malloc(NnTotal*Nvar * sizeof *Wc);   // free
	Uc   = malloc(NnTotal*Nvar * sizeof *Uc);   // free
	Wctr = malloc(NnTotal*Nvar * sizeof *Wctr); // free
	Uctr = malloc(NnTotal*Nvar * sizeof *Uctr); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		Wc[i] = Wr[i];

	convert_variables(Wr,Ur,d,d,Nn,Nel,'c','p');
	convert_variables(Ur,Wr,d,d,Nn,Nel,'p','c');

	convert_variables_c(Wc,Uc,d,d,Nn,Nel,'c','p');
	convert_variables_c(Uc,Wc,d,d,Nn,Nel,'p','c');

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++) {
		Wctr[i] = creal(Wc[i]);
		Uctr[i] = creal(Uc[i]);
	}

	if (array_norm_diff_d(NnTotal*Nvar,Wr,Wctr,"Inf") < EPS &&
	    array_norm_diff_d(NnTotal*Nvar,Ur,Uctr,"Inf") < EPS)
			pass = 1;

	free(Ur);
	free(Wc);
	free(Uc);
	free(Wctr);
	free(Uctr);

	return pass;
}
