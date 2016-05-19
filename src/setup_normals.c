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
	unsigned int fn, fnMax, curved, dim, dim1, dim2,
	             VfIn, Eclass, IndFType, fIn,
	             NvnC0, NvnI0, NfnI0, NnI;
	double       nSum, nSum2,
	             *C_fI, *C_vC, *nrIn, *n;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VIn, *VOut;

	OPS = malloc(sizeof *OPS); // free

	VIn  = FACET->VIn;
	VOut = FACET->VOut;

	VfIn = FACET->VfIn;
	curved = FACET->curved;

	fIn = VfIn/NfrefMax;

	Eclass = get_Eclass(VIn->type);
	IndFType = get_IndFType(Eclass,fIn);

	init_ops(OPS,VIn,FACET,IndFType);

	NvnC0 = OPS->NvnC;
	NvnI0 = OPS->NvnI;
	NfnI0 = OPS->NfnI;

	C_vC = VIn->C_vC;
	C_fI = malloc(NvnI0*d*d * sizeof *C_fI); // free

	mm_CTN_d(NfnI0,d*d,NvnC0,OPS->I_vC_fI[VfIn],C_vC,C_fI);

	NnI = NfnI0;
	nrIn = &(OPS->nr[fIn*d]);

	// Store a single normal on straight FACETs
	if (!curved) fnMax = 1;
	else         fnMax = NnI;

	n = calloc(fnMax*d , sizeof *n); // keep
	for (fn = 0; fn < fnMax; fn++) {
		for (dim1 = 0; dim1 < d; dim1++) {
		for (dim2 = 0; dim2 < d; dim2++) {
			n[fn*d+dim1] += nrIn[dim2]*C_fI[NnI*(d*dim1+dim2)+fn];
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

//printf("%d %d %d\n",FACET->indexg,VfIn,IndFType);
//array_print_d(fnMax,d,n,'R');

	free(C_fI);
	free(OPS);
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;
	FtypeInt = FACET->typeInt;

	ELEMENT     = get_ELEMENT_type(Vtype);
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