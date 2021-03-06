// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_equivalence_real_complex.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "test_code_fluxes.h"
#include "test_code_boundary_conditions.h"
#include "test_support.h"
#include "initialize_test_case.h"
#include "fluxes_structs.h"
#include "fluxes_inviscid.h"
#include "fluxes_inviscid_c.h"
#include "fluxes_viscous.h"
#include "fluxes_viscous_c.h"
#include "boundary_conditions.h"
#include "boundary_conditions_c.h"
#include "variable_functions.h"
#include "variable_functions_c.h"

#include "array_norm.h"
#include "array_free.h"
#include "array_print.h"

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

static unsigned int compare_flux_inviscid(unsigned int const Nn, unsigned int const Nel, double const *const Wr);
static unsigned int compare_flux_viscous  (unsigned int const Nn, unsigned int const Nel, unsigned int const d,
                                           double const *const Wr, double const *const *const Qr);
static unsigned int compare_flux_Num      (unsigned int const Nn, unsigned int const Nel, double const *const nL,
                                           double const *const XYZ, double const *const Wr);
static unsigned int compare_boundary      (const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                           double *const WLr, double *const *const QLr, double *nL, double *XYZ,
                                           const char *BType);
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

	test_print_warning("Ensure that tests for real/complex equivalence for all relevant functions are implemented");

	unsigned int Nn, Nel, d, Neq;
	double       *W, *nL, *XYZ;

	set_memory_test_jacobians('A');

	unsigned int NBTypes = 0;
	char         **BType;
	set_BTypes(&NBTypes,&BType); // free

	unsigned int NFiTypes = 0;
	char         **FiType;
	set_FiTypes(&NFiTypes,&FiType); // free

	unsigned int NFNumTypes = 0;
	char         **FNumType;
	set_FNumTypes(&NFNumTypes,&FNumType); // free

	strcpy(DB.MeshType,"ToBeCurved"); // Meshes are not used for this test (and thus need not exist)

	unsigned int const dMin = 2, dMax = 3;
	for (d = dMin; d <= dMax; d++) {
		W   = initialize_W(&Nn,&Nel,d); // free
		double **Q = initialize_Q(Nn,Nel,d); // free
		nL  = initialize_n(Nn,Nel,d);   // free
		XYZ = initialize_XYZ(Nn,Nel,d); // free

		// flux_inviscid
		for (size_t i = 0; i < NFiTypes; i++) {
			set_parameters_test_flux_inviscid(FiType[i],d);
			initialize_test_case_parameters();

			pass = compare_flux_inviscid(Nn,Nel,W);
			if (i == 0) {
				if (d == dMin)
					sprintf(PrintName,"equivalence_flux_%s             (d = %d):",FiType[i],d);
				else
					sprintf(PrintName,"            flux_%s             (d = %d):",FiType[i],d);
			} else {
				sprintf(PrintName,"                 %s                    :",FiType[i]);
			}
			test_print2(pass,PrintName);

			free(DB.SolverType); // Initialized in "initialize_test_case_parameters"
		}

		// flux_num_inviscid
		for (size_t i = 0; i < NFNumTypes; i++) {
			set_parameters_test_flux_Num(FNumType[i],d);
			initialize_test_case_parameters();

			pass = compare_flux_Num(Nn,Nel,nL,XYZ,W);
			if (i == 0)
				sprintf(PrintName,"                 Num_%s                  :",FNumType[i]);
			else
				sprintf(PrintName,"                     %s                  :",FNumType[i]);
			test_print2(pass,PrintName);

			free(DB.SolverType); // Initialized in "initialize_test_case_parameters"
		}

		Neq  = d+2;

		// flux_viscous
		DB.PDE_index = PDE_NAVIERSTOKES;
		pass = compare_flux_viscous(Nn,Nel,d,W,(double const *const *const) Q);
		sprintf(PrintName,"                 viscous                      :");
		test_print2(pass,PrintName);

		for (size_t i = 0; i < NBTypes; i++) {
			set_parameters_test_boundary_conditions(BType[i],d);
			initialize_test_case_parameters();

			reset_entered_test_boundary_conditions(BType[i]);

			if (strstr(BType[i],"BackPressure"))
				update_values_BackPressure(Nn,Nel,W,nL,d);
			pass = compare_boundary(Nn,Nel,d,W,Q,nL,XYZ,BType[i]);
			if (i == 0)
				sprintf(PrintName,"            boundary_%s        :",BType[i]);
			else
				sprintf(PrintName,"                     %s        :",BType[i]);
			test_print2(pass,PrintName);

			free(DB.SolverType); // Initialized in "initialize_test_case_parameters"
		}

		// convert_variables
		pass = compare_variables(Nn,Nel,d,Neq,W);
		sprintf(PrintName,"            convert_variables                 :");
		test_print2(pass,PrintName);

		free(W);
		array_free2_d(d,Q);
		free(nL);
		free(XYZ);
	}

	set_memory_test_jacobians('F');

	array_free2_c(NBTypes,BType);
	array_free2_c(NFiTypes,FiType);
	array_free2_c(NFNumTypes,FNumType);

	free(PrintName);
}

static unsigned int compare_flux_inviscid(unsigned int const Nn, unsigned int const Nel, double const *const Wr)
{
	unsigned int pass = 0;

	unsigned int const d       = DB.d,
	                   Nvar    = DB.Nvar,
	                   Neq     = DB.Neq,
	                   NnTotal = Nn*Nel;

	double complex *const Wc = malloc(NnTotal*Nvar  * sizeof *Wc), // free
	               *const Fc = malloc(NnTotal*d*Neq * sizeof *Fc); // free
	double *const Fr   = malloc(NnTotal*d*Neq * sizeof *Fr),   // free
	       *const Fctr = malloc(NnTotal*d*Neq * sizeof *Fctr); // free

	for (size_t i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		Wc[i] = Wr[i];

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free

	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = Nn;
	FLUXDATA->Nel = Nel;
	FLUXDATA->W   = Wr;
	FLUXDATA->F   = Fr;
	FLUXDATA->W_c = Wc;
	FLUXDATA->F_c = Fc;

	flux_inviscid(FLUXDATA);
	flux_inviscid_c(FLUXDATA);
	free(FLUXDATA);

	for (size_t i = 0, iMax = NnTotal*d*Neq; i < iMax; i++)
		Fctr[i] = creal(Fc[i]);

	if (array_norm_diff_d(NnTotal*d*Neq,Fr,Fctr,"Inf") < EPS)
		pass = 1;

	free(Wc);
	free(Fc);
	free(Fr);
	free(Fctr);

	return pass;
}

static unsigned int compare_flux_viscous(unsigned int const Nn, unsigned int const Nel, unsigned int const d,
                                         double const *const Wr, double const *const *const Qr)
{
	DB.d        = d;
	DB.Pr       = 0.72;
	DB.mu       = 1e-0;
	DB.Const_mu = 1;

	unsigned int pass = 0;

	unsigned int const Neq     = d+2,
	                   Nvar    = d+2,
	                   NnTotal = Nn*Nel;

	double complex *const Wc = malloc(NnTotal*Nvar * sizeof *Wc); // free
	for (size_t i = 0; i < NnTotal*Nvar; i++)
		Wc[i] = Wr[i];

	double complex **const Qc = malloc(d * sizeof *Qc); // free
	for (size_t dim = 0; dim < d; dim++) {
		Qc[dim] = malloc(NnTotal*Nvar * sizeof *Qc[dim]); // free
		for (size_t i = 0; i < NnTotal*Nvar; i++)
			Qc[dim][i] = Qr[dim][i];
	}

	double         *const Fr   = malloc(NnTotal*d*Neq * sizeof *Fr),   // free
	               *const Fctr = malloc(NnTotal*d*Neq * sizeof *Fctr); // free
	double complex *const Fc   = malloc(NnTotal*d*Neq * sizeof *Fc);   // free

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = Nn;
	FLUXDATA->Nel = Nel;
	FLUXDATA->W   = Wr;
	FLUXDATA->Q   = Qr;
	FLUXDATA->F   = Fr;
	FLUXDATA->W_c = Wc;
	FLUXDATA->Q_c = (double complex const *const *const) Qc;
	FLUXDATA->F_c = Fc;

	flux_viscous(FLUXDATA);
	flux_viscous_c(FLUXDATA);
	free(FLUXDATA);

	for (size_t i = 0; i < NnTotal*d*Neq; i++)
		Fctr[i] = creal(Fc[i]);

	if (array_norm_diff_d(NnTotal*d*Neq,Fr,Fctr,"Inf") < EPS)
		pass = 1;

	free(Wc);
	array_free2_cmplx(d,Qc);
	free(Fr);
	free(Fctr);
	free(Fc);

	return pass;
}

static unsigned int compare_flux_Num(unsigned int const Nn, unsigned int const Nel, double const *const nL,
                                     double const *const XYZ, double const *const Wr)
{
	if (Nel != 2)
		EXIT_UNSUPPORTED;

	unsigned int pass = 0;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	double *const WLr   = malloc(Nn*Nvar * sizeof *WLr),   // free
	       *const WRr   = malloc(Nn*Nvar * sizeof *WRr),   // free
	       *const nFr   = malloc(Nn*Nvar * sizeof *nFr),   // free
	       *const nFctr = malloc(Nn*Nvar * sizeof *nFctr); // free
	double complex *const WLc = malloc(Nn*Nvar * sizeof *WLc), // free
	               *const WRc = malloc(Nn*Nvar * sizeof *WRc), // free
	               *const nFc = malloc(Nn*Nvar * sizeof *nFc); // free

	double const *W_ptr = Wr;
	for (size_t var = 0; var < Nvar; var++) {
		size_t const IndWLR = var*Nn;
		for (size_t n = 0; n < Nn; n++)
			WLr[IndWLR+n] = *W_ptr++;
		for (size_t n = 0; n < Nn; n++)
			WRr[IndWLR+n] = *W_ptr++;
	}

	for (size_t i = 0, iMax = Nn*Nvar; i < iMax; i++) {
		WLc[i] = WLr[i];
		WRc[i] = WRr[i];
	}

	struct S_NUMERICALFLUX *const NUMFLUXDATA = malloc(sizeof *NUMFLUXDATA); // free

	NUMFLUXDATA->NumFluxInviscid_index = DB.InviscidFluxType;
	NUMFLUXDATA->d   = d;
	NUMFLUXDATA->Nn  = Nn;
	NUMFLUXDATA->Nel = 1;
	NUMFLUXDATA->nL  = nL;
	NUMFLUXDATA->XYZ = XYZ;

	NUMFLUXDATA->WL         = WLr;
	NUMFLUXDATA->WR         = WRr;
	NUMFLUXDATA->nFluxNum   = nFr;
	NUMFLUXDATA->WL_c       = WLc;
	NUMFLUXDATA->WR_c       = WRc;
	NUMFLUXDATA->nFluxNum_c = nFc;

	flux_num_inviscid(NUMFLUXDATA);
	flux_num_inviscid_c(NUMFLUXDATA);
	free(NUMFLUXDATA);

	for (size_t i = 0, iMax = Nn*Neq; i < iMax; i++)
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
                                     double *const WLr, double *const *const QLr, double *nL, double *XYZ,
                                     const char *BType)
{
	unsigned int pass = 0;

	unsigned int const NnTotal = Nn*Nel,
	                   Nvar    = DB.Nvar;

	double *const WBr   = malloc(NnTotal*Nvar * sizeof *WBr),   // free
	       *const WBctr = malloc(NnTotal*Nvar * sizeof *WBctr); // free

	double complex *const WLc = malloc(NnTotal*Nvar * sizeof *WLc),   // free
	               *const WBc = malloc(NnTotal*Nvar * sizeof *WBc);   // free

	for (size_t i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WLc[i] = WLr[i];

	struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

	set_BC_from_BType(BCdata,BType);

	BCdata->d   = d;
	BCdata->Nn  = Nn;
	BCdata->Nel = Nel;

	BCdata->XYZ  = XYZ;
	BCdata->nL   = nL;
	BCdata->WL   = WLr;
	BCdata->WB   = WBr;
	BCdata->WL_c = WLc;
	BCdata->WB_c = WBc;

	double **QBr = NULL, **QBctr = NULL;
	double complex **QLc = NULL, **QBc = NULL;
	if (BCdata->ComputeQ) {
		QBr   = malloc(d * sizeof *QBr);   // free
		QBctr = malloc(d * sizeof *QBctr); // free
		QLc   = malloc(d * sizeof *QLc);   // free
		QBc   = malloc(d * sizeof *QBc);   // free
		for (size_t dim = 0; dim < d; dim++) {
			QBr[dim]   = malloc(NnTotal*Nvar * sizeof *QBr[dim]);   // free
			QBctr[dim] = malloc(NnTotal*Nvar * sizeof *QBctr[dim]); // free
			QLc[dim]   = malloc(NnTotal*Nvar * sizeof *QLc[dim]);   // free
			QBc[dim]   = malloc(NnTotal*Nvar * sizeof *QBc[dim]);   // free

			for (size_t i = 0; i < NnTotal*Nvar; i++)
				QLc[dim][i] = QLr[dim][i];
		}

		BCdata->QL   = (double const *const *const) QLr;
		BCdata->QB   = QBr;
		BCdata->QL_c = (double complex const *const *const) QLc;
		BCdata->QB_c = QBc;
	}

	correct_XYZ_for_exact_normal(BCdata,BType);
	compute_boundary_values(BCdata);
	compute_boundary_values_c(BCdata);

	for (size_t i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WBctr[i] = creal(WBc[i]);

	if (!BCdata->ComputeQ) {
		if (array_norm_diff_d(NnTotal*Nvar,WBr,WBctr,"Inf") < 1e1*EPS) {
			pass = 1;
		} else {
			printf("% .3e\n",array_norm_diff_d(NnTotal*Nvar,WBr,WBctr,"Inf"));
			array_print_d(NnTotal,Nvar,WBr,'C');
			array_print_d(NnTotal,Nvar,WBctr,'C');
		}
	} else {
		double diff = array_norm_diff_d(NnTotal*Nvar,WBr,WBctr,"Inf");
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NnTotal*Nvar; i++)
				QBctr[dim][i] = creal(QBc[dim][i]);
			diff += array_norm_diff_d(NnTotal*Nvar,QBr[dim],QBctr[dim],"Inf");
		}

		if (diff < 1e1*EPS) {
			pass = 1;
		} else {
			for (size_t i = 0; i < NnTotal*Nvar; i++)
				WBr[i] -= WBctr[i];
			array_print_d(NnTotal,Nvar,WBr,'C');

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0; i < NnTotal*Nvar; i++)
					QBr[dim][i] -= QBctr[dim][i];
				printf("%zu\n",dim);
				array_print_d(NnTotal,Nvar,QBr[dim],'C');
			}
		}

		array_free2_d(d,QBr);
		array_free2_d(d,QBctr);
		array_free2_cmplx(d,QLc);
		array_free2_cmplx(d,QBc);
	}
	free(BCdata);

	// Ensure that all settings were entered for the boundary conditions (if applicable)
	bool CheckedAll = 0;
	check_entered_test_boundary_conditions(&CheckedAll,BType);

	if (!CheckedAll) {
		pass = 0;
		printf("Did not check all boundary condition settings for boundary %s\n",BType);
	}

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
