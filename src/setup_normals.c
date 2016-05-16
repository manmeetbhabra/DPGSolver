#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up normals at integration nodes on element facets.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
 */

struct S_OPERATORS {
	unsigned int NvnC, NvnI, NfnI;
	double       *I_vC_vI, **I_vC_fI, *nr;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass);

void setup_normals(struct S_FACET *FACET)
{
	// Initialize DB Parameters
	unsigned int d        = DB.d,
	             NfrefMax = DB.NfrefMax;

	// Standard datatypes
	unsigned int i, fn, fnMax, curved, dim, dim1, dim2,
	             VfIn,
	             NvnC0, NvnI0, NfnI0, NvnC1, NvnI1, NfnI1, NnI;
	double       nSum, nSum2,
	             *C_fI, *C_vC, *nrIn, *n;

	struct S_OPERATORS *OPS[2];
	struct S_VOLUME    *VIn, *VOut;

	// silence
//	fnMax = 0;
	NnI = 0;
	C_fI = NULL;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	VIn  = FACET->VIn;
	VOut = FACET->VOut;

	VfIn = FACET->VfIn;
	curved = FACET->curved;

	init_ops(OPS[0],VIn,FACET,0);
	if (VIn->type == WEDGE)
		init_ops(OPS[1],VIn,FACET,1);

	NvnC0 = OPS[0]->NvnC;
	NvnI0 = OPS[0]->NvnI;
	NfnI0 = OPS[0]->NfnI;

	C_vC = VIn->C_vC;
	if (VIn->Eclass == C_TP) {
		printf("Add in support for C_TP in setup_normals.\n");
		exit(1);
	} else if (VIn->Eclass == C_SI || VIn->Eclass == C_PYR) {

		C_fI = malloc(NvnI0*d*d * sizeof *C_fI); // free

		mm_CTN_d(NfnI0,d*d,NvnC0,OPS[0]->I_vC_fI[VfIn],C_vC,C_fI);

		NnI = NfnI0;
	} else if (VIn->Eclass == C_WEDGE) {
		NvnC1 = OPS[1]->NvnC;
		NvnI1 = OPS[1]->NvnI;
		NfnI1 = OPS[1]->NfnI;

		printf("Add in support for C_WEDGE in setup_normals.\n");
		exit(1);
	}
	nrIn = &OPS[0]->nr[(VfIn/NfrefMax)*d];

	// Store a single normal on straight FACETs
	if (!curved) fnMax = 1;
	else         fnMax = NnI;

	n = calloc(fnMax*d , sizeof *n); // keep
	for (fn = 0; fn < fnMax; fn++) {
		for (dim1 = 0; dim1 < d; dim1++) {
		for (dim2 = 0; dim2 < d; dim2++) {
			n[fn*d+dim1] += nrIn[dim2]*C_fI[NnI*(dim1+d*dim2)+fn];
		}}
	}

	for (fn = 0; fn < fnMax; fn++) {
		nSum2 = 0;
		for (dim = 0; dim < d; dim++)
			nSum2 += pow(n[fn*d+dim],2.0);

		nSum = sqrt(nSum2); // == detJF_vI
		for (dim = 0; dim < d; dim++)
			n[fn*d+dim] /= nSum;
	}

	FACET->n = n;

//printf("%d %d\n",FACET->indexg,FACET->curved);
//array_print_d(fnMax,d,n,'R');

	free(C_fI);
	for (i = 0; i < 2; i++)
		free(OPS[i]);
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	// silence
	ELEMENT_OPS = NULL;

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;
	FtypeInt = FACET->typeInt;

	ELEMENT = get_ELEMENT_type(Vtype);
	if (Vtype == LINE || Vtype == QUAD || Vtype == HEX || Vtype == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else if (Vtype == TRI || Vtype == TET || Vtype == PYR)
		ELEMENT_OPS = ELEMENT;

	if (!Vcurved) {
		// Straight VOLUME
		OPS->NvnC = ELEMENT_OPS->NvnCs[PV];
		if (FtypeInt == 's') {
			// Straight FACET Integration
			OPS->NvnI = ELEMENT_OPS->NvnIs[PV];
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCs_vIs[PV];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCs_fIs[PV][PF];
		} else {
			// Curved FACET Integration
			OPS->NvnI = ELEMENT_OPS->NvnIc[PV];
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCs_vIc[PV];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCs_fIc[PV][PF];
		}
	} else {
		// Curved VOLUME
		OPS->NvnC = ELEMENT_OPS->NvnCc[PV];
		if (FtypeInt == 's') {
			// Straight FACET Integration
			OPS->NvnI = ELEMENT_OPS->NvnIs[PV];
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCc_vIs[PV];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCc_fIs[PV][PF];
		} else {
			// Curved FACET Integration
			OPS->NvnI = ELEMENT_OPS->NvnIc[PV];
			OPS->NfnI = ELEMENT_OPS->NfnIc[PV][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCc_vIc[PV];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCc_fIc[PV][PF];
		}
	}
	OPS->nr = ELEMENT->nr;
}
