// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "implicit_VOLUME_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_OpCSR.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "jacobian_fluxes_inviscid.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the LHS term.
 *
 *	Comments:
 *		Certain multiplications can be avoided when computing either Fr from F or when computing LHS terms based on the
 *		sparsity of the flux Jacobian. Test performance improvement if these terms are neglected (ToBeModified).
 *		As the LHS terms for each (eq,var) combination are stored as row-major arrays, mm_d is used to compute these
 *		terms for the uncollocated scheme. Potentially investigate whether a custom mm implementation would be faster
 *		than the BLAS call (ToBeModified).
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnI, NvnS, NvnS_SF, NvnI_SF;
	double       *ChiS_vI, **D_Weak, *I_Weak;

	struct S_OpCSR **D_Weak_sp;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass);
static void compute_VOLUME_LHS_EFE(void);

void implicit_VOLUME_info(void)
{
	// Initialize DB Parameters
	unsigned int EFE        = DB.EFE,
	             Vectorized = DB.Vectorized;

	if (EFE) {
		switch (Vectorized) {
		case 0:
			compute_VOLUME_LHS_EFE();
			break;
		default:
			printf("Error: Unsupported Vectorized (%d).\n",Vectorized), EXIT_MSG;
			break;
		}
	} else {
		;
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int ***SF_BE = DB.SF_BE;

	// Standard datatypes
	unsigned int P, type, curved, Eclass;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	// silence
	ELEMENT_OPS = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;
	Eclass = VOLUME->Eclass;

	ELEMENT = get_ELEMENT_type(type);
	if ((Eclass == C_TP && SF_BE[P][0][0]) || (Eclass == C_WEDGE && SF_BE[P][1][0]))
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else
		ELEMENT_OPS = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[P];
	OPS->NvnS_SF = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->NvnI    = ELEMENT->NvnIs[P];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIs[P];

		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P][P][0];
		OPS->D_Weak  = ELEMENT_OPS->Ds_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT_OPS->Is_Weak_VV[P][P][0];

		OPS->D_Weak_sp = ELEMENT->Ds_Weak_VV_sp[P][P][0];
	} else {
		OPS->NvnI    = ELEMENT->NvnIc[P];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIc[P];

		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P][P][0];
		OPS->D_Weak  = ELEMENT_OPS->Dc_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT_OPS->Ic_Weak_VV[P][P][0];

		OPS->D_Weak_sp = ELEMENT->Dc_Weak_VV_sp[P][P][0];
	}
}

static void compute_VOLUME_LHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
				 Nvar       = DB.Nvar,
				 Neq        = DB.Neq,
	             ***SF_BE   = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, eq, var, dim1, dim2, P, iMax,
	             Indeqvar, InddFrdW, InddFdW, IndC, IndD, IndLHS, Eclass,
	             NvnI, NvnS, NvnI_SF[2], NvnS_SF[2], NIn[3], NOut[3], Diag[3];
	double       *W_vI, *dFdW_vI, *dFrdW_vI, *C_vI, *LHS, *DdFrdW, **D, *Ddim, *OP[3];

//	struct S_OpCSR **D_sp;

	struct S_OPERATORS *OPS[2];
	struct S_VOLUME    *VOLUME;

	// silence
	NvnI_SF[1] = 0;
	NvnS_SF[1] = 0;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(Form,"Weak")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			P = VOLUME->P;

			// Obtain operators
			init_ops(OPS[0],VOLUME,0);
			if (VOLUME->type == WEDGE)
				init_ops(OPS[1],VOLUME,1);

			Eclass = VOLUME->Eclass;

			// Obtain W_vI
			NvnI = OPS[0]->NvnI;
			if (Collocated) {
				W_vI = VOLUME->What;
			} else {
				W_vI = malloc(NvnI*Nvar * sizeof *W_vI); // free

				if (Eclass == C_TP && SF_BE[P][0][0]) {
					for (i = 0; i < 1; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,0,0,NULL,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;

					sf_apply_d(VOLUME->What,W_vI,NIn,NOut,Nvar,OP,Diag,d);
				} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
					for (i = 0; i < 2; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,
					                  NvnS_SF[1],NvnI_SF[1],OPS[1]->ChiS_vI,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;
					Diag[1] = 2;

					sf_apply_d(VOLUME->What,W_vI,NIn,NOut,Nvar,OP,Diag,d);
				} else {
					mm_CTN_d(NvnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_vI,VOLUME->What,W_vI);
				}
			}
/*
if (VOLUME->indexg == 0) {
printf("VOLUME %d\n",VOLUME->indexg);
array_print_d(OPS[0]->NvnS,Nvar,VOLUME->What,'C');
array_print_d(NvnI,Nvar,W_vI,'C');
}
*/

			// Compute Flux in reference space
			dFdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFdW_vI); // free
			jacobian_flux_inviscid(NvnI,1,W_vI,dFdW_vI,d,Neq);
/*
for (eq = 0; eq < Neq; eq++) {
for (var = 0; var < Nvar; var++) {
	printf("%d %d\n",eq,var);
	array_print_d(NvnI,d,&dFdW_vI[(eq*Nvar+var)*NvnI*d],'C');
}}
EXIT_MSG;
*/

			if (!Collocated)
				free(W_vI);

			C_vI = VOLUME->C_vI;
//array_print_d(NvnI,d,C_vI,'C');
//array_print_d(NvnI,d,&C_vI[NvnI*d],'C');
//EXIT_MSG;

			dFrdW_vI = calloc(NvnI*d*Nvar*Neq , sizeof *dFrdW_vI); // free
			for (eq = 0; eq < Neq; eq++) {
			for (var = 0; var < Nvar; var++) {
				Indeqvar = (eq*Nvar+var)*d;
				for (dim1 = 0; dim1 < d; dim1++) {
					InddFrdW = (Indeqvar+dim1)*NvnI;
					for (dim2 = 0; dim2 < d; dim2++) {
						InddFdW = (Indeqvar+dim2)*NvnI;
						IndC    = (dim1*d+dim2)*NvnI;
						for (i = 0; i < NvnI; i++)
							dFrdW_vI[InddFrdW+i] += dFdW_vI[InddFdW+i]*C_vI[IndC+i];
					}
				}
			}}
			free(dFdW_vI);
/*
for (eq = 0; eq < Neq; eq++) {
for (var = 0; var < Nvar; var++) {
	printf("%d %d\n",eq,var);
	array_print_d(NvnI,d,&dFrdW_vI[(eq*Nvar+var)*NvnI*d],'C');
}}
EXIT_MSG;
*/

			// Compute LHS terms
			NvnS = OPS[0]->NvnS;

			if (VOLUME->LHS)
				free(VOLUME->LHS);
			LHS = calloc(NvnS*NvnS*Neq*Nvar , sizeof *LHS); // keep (requires external free)
			VOLUME->LHS = LHS;

			if (Collocated && (Eclass == C_TP || Eclass == C_WEDGE)) {
// Note: Collocated here implies that ChiS_vI == Identity
				printf("Error: Modifications required.\n"), EXIT_MSG;

/*
				D_sp = OPS[0]->D_Weak_sp;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_CSR_d(NvnS,1,NvnI,D_sp[dim1],&Fr_vI[(eq*d+dim1)*NvnI],DFr);

					IndRHS = eq*NvnS;
					for (i = 0; i < NvnS; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
*/
			} else {
				D = OPS[0]->D_Weak;

				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					Indeqvar = (eq*Nvar+var)*d;

					DdFrdW = calloc(NvnS*NvnI , sizeof *DdFrdW); // free
					for (dim1 = 0; dim1 < d; dim1++) {
						Ddim = D[dim1];
//array_print_d(OPS[0]->NvnS,NvnI,Ddim,'R');

						InddFrdW = (Indeqvar+dim1)*NvnI;
						for (i = 0; i < NvnS; i++) {
							IndD = i*NvnI;
							for (j = 0; j < NvnI; j++)
								DdFrdW[IndD+j] += Ddim[IndD+j]*dFrdW_vI[InddFrdW+j];
						}
					}
//EXIT_MSG;

					IndLHS = (eq*Nvar+var)*NvnS*NvnS;
					if (Collocated) {
						for (i = 0, iMax = NvnS*NvnS; i < iMax; i++)
							LHS[IndLHS+i] = DdFrdW[i];
					} else {
						mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,DdFrdW,OPS[0]->ChiS_vI,&LHS[IndLHS]);
					}
					free(DdFrdW);
				}}
			}
			free(dFrdW_vI);
/*
if (VOLUME->indexg == 0) {
printf("VOLUME %d\n",VOLUME->indexg);
for (eq = 0; eq < Neq; eq++) {
for (var = 0; var < Nvar; var++) {
	printf("%d %d\n",eq,var);
	array_print_d(NvnS,NvnS,&LHS[(eq*Nvar+var)*NvnS*NvnS],'R');
}}
}
EXIT_MSG;
*/
		}
	} else if (strstr(Form,"Strong")) {
		printf("Exiting: Implement the strong form in compute_VOLUME_RHS_EFE.\n"), exit(1);
	}

	for (i = 0; i < 2; i++)
		free(OPS[i]);
}