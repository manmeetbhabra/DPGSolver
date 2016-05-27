// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up VOLUME and FACET structures.
 *
 *	Comments:
 *		Will need two different setup_structures functions: One using the initial global arrays and one updating
 *		elements individually after hp refinement. (ToBeDeleted)
 *		Probably better to change the name of XYZc to XYZ_vC to avoid confusion that this is XYZc(urved). (ToBeDeleted)
 *
 *	Notation:
 *		XYZ_(1)(2) : Physical node locations (XYZ) of (1) nodes of (2) type
 *		             (1) : (v)olume, (f)acet
 *		             (2) : (C)orner
 *		IndOrd(In/Out)(Out/In) : (Ind)ex of (Ord)ering relating projection from (In)ner VOLUME to the FACET to the
 *		                         projection from the (Out)er VOLUME to the FACET
 *
 *	References:
*/

static void compute_distance_matrix(const unsigned int Nn, const unsigned int BC, const unsigned int d, double *XYZIn,
                                    double *XYZOut, double *DXYZ)
{
	// Standard datatypes
	unsigned int i, j, k, kMax,
	             tmp_ui, IndDXYZ, IndComp[2];
	double       tmp_d;

	if (d == 1)
		return;

	for (i = 0; i < Nn; i++) {
	for (j = 0; j < Nn; j++) {
		tmp_ui = BC % BC_STEP_SC;
		if (tmp_ui > BC_PERIODIC_MIN) { // special case for periodic
			if (d == 3) {
				if      ((tmp_ui - BC_PERIODIC_MIN)/2 == 1) IndComp[0] = 1, IndComp[1] = 2; // Periodic in X
				else if ((tmp_ui - BC_PERIODIC_MIN)/2 == 2) IndComp[0] = 0, IndComp[1] = 2; // Periodic in Y
				else if ((tmp_ui - BC_PERIODIC_MIN)/2 == 3) IndComp[0] = 0, IndComp[1] = 1; // Periodic in Z
			} else if (d == 2) {
				if      ((tmp_ui - BC_PERIODIC_MIN)/2 == 1) IndComp[0] = 1; // Periodic in X
				else if ((tmp_ui - BC_PERIODIC_MIN)/2 == 2) IndComp[0] = 0; // Periodic in Y
			}
			IndDXYZ = i*Nn+j;
			for (k = 0, kMax = d-1; k < kMax; k++) {
				tmp_d = fabs(XYZIn[i*d+IndComp[k]]-XYZOut[j*d+IndComp[k]]);
				if (tmp_d > DXYZ[IndDXYZ])
					DXYZ[IndDXYZ] = tmp_d;
			}
		} else {
			DXYZ[i*Nn+j] = array_norm_diff_d(d,&XYZIn[i*d],&XYZOut[j*d],"Inf");
		}
	}}
}

static void get_ordering_index(const unsigned int Nn, const unsigned int d, double *DXYZ, unsigned int *IndOrdInOut,
                               unsigned int *IndOrdOutIn)
{
	/*
	 *	Purpose:
	 *		Return the ordering index corresponding to the match between two FACETs based on DXYZ.
	 *
	 *	Comments:
	 *		In 1D, IndOrdInOut == IndOrdOutIn == 0.
	 *		In 2D, IndOrdInOut == IndOrdOutIn.
	 *
	 *	Notation:
	 *		IndZerosP : (Ind)ices of (Zeros) which are (P)ossible.
	 */

	if (d == 1) {
		*IndOrdInOut = 0;
		*IndOrdOutIn = 0;
		return;
	} else {
		unsigned int i, j,
					 IndZerosInOut[Nn], IndZerosOutIn[Nn];
		double       zero[1] = {0.0};

		// Find indices of zeros in DXYZ
		for (i = 0; i < Nn; i++) {
		for (j = 0; j < Nn; j++) {
//printf("%d %d\n",i,j);
			if (array_norm_diff_d(1,&DXYZ[i*Nn+j],zero,"Inf") < EPS) {
				IndZerosInOut[i] = j;
				break;
			}
		}}
//array_print_ui(1,Nn,IndZerosInOut,'R');

		if (d == 3) {
			// Transpose DXYZ and find Out->In Ordering as well
			mkl_dimatcopy('R','T',Nn,Nn,1.0,DXYZ,Nn,Nn);

			// Find indices of zeros in DXYZ'
			for (i = 0; i < Nn; i++) {
			for (j = 0; j < Nn; j++) {
				if (array_norm_diff_d(1,&DXYZ[i*Nn+j],zero,"Inf") < EPS) {
					IndZerosOutIn[i] = j;
					break;
				}
			}}


			if (Nn == 4) { // QUAD FACET
				unsigned int IndZerosP[32] = { 0, 1, 2, 3,
				                               2, 0, 3, 1,
				                               3, 2, 1, 0,
				                               1, 3, 0, 2,
				                               0, 2, 1, 3,
				                               1, 0, 3, 2,
				                               3, 1, 2, 0,
				                               2, 3, 0, 1};

				for (i = 0; i < 8; i++) {
					if (array_norm_diff_ui(Nn,IndZerosInOut,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdInOut = i;
					if (array_norm_diff_ui(Nn,IndZerosOutIn,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdOutIn = i;
				}
			} else if (Nn == 3) { // TRI FACET
				unsigned int IndZerosP[18] = { 0, 1, 2,
				                               2, 0, 1,
				                               1, 2, 0,
				                               0, 2, 1,
				                               1, 0, 2,
				                               2, 1, 0};

				for (i = 0; i < 6; i++) {
					if (array_norm_diff_ui(Nn,IndZerosInOut,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdInOut = i;
					if (array_norm_diff_ui(Nn,IndZerosOutIn,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdOutIn = i;
				}
			}
		} else if (d == 2) {
			unsigned int IndZerosP[4] = { 0, 1,
			                              1, 0};
			for (i = 0; i < 2; i++) {
				if (array_norm_diff_ui(Nn,IndZerosInOut,&IndZerosP[i*Nn],"Inf") < EPS) {
					*IndOrdInOut = i;
					*IndOrdOutIn = i;
					return;
				}
			}
		}
	}
}

void setup_structures(void)
{
	// Initialize DB Parameters
	unsigned int d        = DB.d,
	             NveMax   = DB.NveMax,
	             NfveMax  = DB.NfveMax,
	             NfrefMax = DB.NfrefMax,
	             AC       = DB.AC,
	             P        = DB.P,
	             NP       = DB.NP,
	             NfMax    = DB.NfMax,
	             NV       = DB.NV,
	             NGF      = DB.NGF,
	             NGFC     = DB.NGFC,
	             NVC      = DB.NVC,
	             *NE      = DB.NE,
	             *EToVe   = DB.EToVe,
	             *EType   = DB.EType,
	             *EToPrt  = DB.EToPrt,
//	             *VToV    = DB.VToV,
	             *VToGF   = DB.VToGF,
	             *VToBC   = DB.VToBC,
	             *VC      = DB.VC,
	             *GFC     = DB.GFC;
	int          MPIrank  = DB.MPIrank;
	double       *VeXYZ   = DB.VeXYZ;

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int i, iMax, f, v, dim, ve, gf, curved,
	             IndE, IndVC, IndVgrp, IndGFC, IndVIn, IndVeF, Indf, Indsf, IndOrdInOut, IndOrdOutIn,
	             Vs, vlocal, NVlocal, NECgrp, NVgrp, Vf,
	             Nve,*Nfve, *Nfref, Nfn,
				 indexg, NvnGs,
	             uMPIrank,
	             *GFToV, *GF_Nv;
	double       *XYZ_vC, *VeF, *XYZIn_fC, *XYZOut_fC, *DXYZ;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **Vgrp, **Vgrp_tmp;
	struct S_FACET   **FACET, **FoundFACET;

	// silence
	NECgrp = 0;
	IndOrdInOut = 0;
	IndOrdOutIn = 0;

	uMPIrank = MPIrank;

	if      (d == 1) NECgrp = 1;
	else if (d == 2) NECgrp = 2;
	else if (d == 3) NECgrp = 4;

	FACET      = calloc(2   , sizeof *FACET); // free
	FoundFACET = calloc(NGF , sizeof *FoundFACET); // free

	NVgrp = NECgrp*NP*2;
	Vgrp     = malloc(NVgrp * sizeof *Vgrp);     // keep
	Vgrp_tmp = malloc(NVgrp * sizeof *Vgrp_tmp); // free
	for (i = 0, iMax = NVgrp; iMax--; i++) {
		Vgrp[i]     = NULL;
		Vgrp_tmp[i] = NULL;
	}

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	VOLUME = New_VOLUME();
	DB.VOLUME = VOLUME;
	DB.FACET  = NULL;

	// Note: only initialize volumes on the current processor.
	for (v = 0, NVlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank)
			NVlocal++;
	}

	for (v = 0, IndE = Vs, IndGFC = IndVC = 0, vlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank) {
			// General
			VOLUME->indexl = vlocal;
			VOLUME->indexg = v;
			VOLUME->P      = P;
			VOLUME->type   = EType[IndE];
			VOLUME->Eclass = get_Eclass(VOLUME->type);

			if (AC || v == VC[IndVC]) {
				VOLUME->curved = 1;
				IndVC++;
			} else {
				VOLUME->curved = 0;
			}

			// Connectivity
			/*	Element connectivity numbering will be based on numbering supporting h-refinement; for this reason it may
			 *	seem overly complex while storing the connectivity information on the initially conforming mesh.
			 *	Supported h-refinements include:
			 *		TET: Isotropic refinement into 4 TET, 2 PYR (PYR orientation tbd)
			 *		HEX: Isotropic and anisotropic refinement of opposite facets into 2 horizontal, 2 vertical, 4 equal
			 *		     QUADs.
			 *		WEDGE: Isotropic refinement of opposing TRI facets, isotropic and anisotropic refinement of all QUAD
			 *		       facets.
			 *		PYR: Isotropic refinement into 6 PYR, 4 TET.
			 *	Max number of faces (including all h-refined variations:
			 *		TET: 4*(1+4) = 20
			 *		HEX: 6*(1+2*2+4) = 54
			 *		WEDGE: 2*(1+4)+3*(1+2*2+4) = 37
			 *		PYR: 4*(1+4)+1*(1+2*2+4) = 25
			 *
			 *			MAX FACES: 6*9 = 54
			 */

			// FACETs adjacent to VOLUMEs on the current processor.
			// Note that GFC and VToGF cycle through the global FACET indices in order
			for (f = 0; f < NfMax; f++) {
				gf = VToGF[v*NfMax+f];
				if (FoundFACET[gf] == NULL) {

					if (DB.FACET != NULL) {
						FACET[0]->next = New_FACET();
						FACET[0]       = FACET[0]->next;
					} else {
						DB.FACET = New_FACET();
						FACET[0] = DB.FACET;
					}

//					FACET[0]->indexl = gflocal;
					FACET[0]->indexg = gf;
					FACET[0]->P      = VOLUME->P;

					FACET[0]->VIn   = VOLUME; IndVIn = VOLUME->indexg;
					FACET[0]->VfIn  = NfrefMax*f;

					FACET[0]->BC    = VToBC[IndVIn*NfMax+f];

					// Overwritten if a second VOLUME is found adjacent to this FACET
					FACET[0]->VOut  = VOLUME;
					FACET[0]->VfOut = NfrefMax*f;

					if (!VOLUME->curved) {
						FACET[0]->typeInt = 's';
					} else {
						FACET[0]->typeInt = 'c';
						if (AC || (IndGFC < NGFC && gf == GFC[IndGFC])) {
							FACET[0]->curved = 1;
							IndGFC++;
						}
					}

					FoundFACET[gf] = FACET[0];
				} else {
					FACET[1] = FoundFACET[gf];

					FACET[1]->P = max(FACET[1]->P,VOLUME->P);
					FACET[1]->VOut  = VOLUME;
					FACET[1]->VfOut = NfrefMax*f;
					if (VOLUME->curved) {
						FACET[1]->typeInt = 'c';
						if (AC || (IndGFC < NGFC && gf == GFC[IndGFC])) {
							FACET[1]->curved = 1;
							IndGFC++;
						}
					}
				}
//				// Indexing from connectivity discussion above noting that the mesh is conforming at the start (ToBeDeleted)
//				VOLUME->GF[f*NfrefMax] = FoundFACET[gf];
			}


			// Geometry
			ELEMENT = get_ELEMENT_type(VOLUME->type);

			if (VOLUME->Eclass == C_TP)
				NvnGs = pow(ELEMENT->ELEMENTclass[0]->NvnGs[0],d);
			else if (VOLUME->Eclass == C_WEDGE)
				NvnGs = (ELEMENT->ELEMENTclass[0]->NvnGs[0])*(ELEMENT->ELEMENTclass[1]->NvnGs[0]);
			else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR)
				NvnGs = ELEMENT->NvnGs[0];
			else
				printf("Error: Unsupported element type setup_struct (NvnGs).\n"), exit(1);

			XYZ_vC = malloc(NvnGs*d * sizeof *XYZ_vC); // keep
			VOLUME->XYZ_vC = XYZ_vC;

			// XYZ_vC may be interpreted as [X Y Z] where each of X, Y, Z are column vectors (ToBeDeleted)
			// Add this comment to notation section.
			indexg = VOLUME->indexg;
			for (ve = 0; ve < NvnGs; ve++) {
			for (dim = 0; dim < d; dim++) {
				XYZ_vC[dim*NvnGs+ve] = VeXYZ[EToVe[(Vs+indexg)*8+ve]*d+dim];
			}}
//printf("%d\n",VOLUME->indexg);
//array_print_d(NvnGs,d,XYZ_vC,'C');

			// MPI
			IndVgrp = ((VOLUME->Eclass)*NP*2)+(VOLUME->P*2)+(VOLUME->curved);
			if (Vgrp[IndVgrp] == NULL)
				Vgrp[IndVgrp] = VOLUME;
			else
				Vgrp_tmp[IndVgrp]->grpnext = VOLUME;

			Vgrp_tmp[IndVgrp] = VOLUME;

			if (vlocal != NVlocal-1) {
				VOLUME->next = New_VOLUME();
				VOLUME = VOLUME->next;
			}

			vlocal++;
		} else {
			// Ensure that appropriate global indices are incremented if necessary
			if (v == VC[IndVC])
				IndVC++;
			for (f = 0; f < NfMax; f++) {
				gf = VToGF[v*NfMax+f];
				if (IndGFC < NGFC && gf == GFC[IndGFC])
					IndGFC++;
			}
		}

		IndE++;
	}
	free(Vgrp_tmp);
	free(FoundFACET);

	if (!AC && IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);

	for (FACET[0] = DB.FACET; FACET[0] != NULL; FACET[0] = FACET[0]->next) {
// May potentially have a problem for PYR-HEX interface due to PYR nodes being ordered for symmetry while QUAD nodes are
// order for TP extension. (ToBeDeleted)

		// Obtain XYZIn_fC/XYZOut_fC
		VOLUME  = FACET[0]->VIn;
printf("%d\n",VOLUME->indexg);
		Vf      = FACET[0]->VfIn;
		Indf    = Vf / NfrefMax; // face index (ToBeDeleted: Move to notation)
		Indsf   = Vf % NfrefMax; // sub face index (ToBeDeleted: Move to notation)

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nve   = ELEMENT->Nve;
		Nfve  = ELEMENT->Nfve;
		Nfref = ELEMENT->Nfref;
		VeF   = ELEMENT->VeF;

		NvnGs = ELEMENT->NvnGs[0];
		XYZ_vC = VOLUME->XYZ_vC;

		for (i = IndVeF = 0; i < Indsf; i++)
			IndVeF += Nfref[i]*Nfve[i];
		IndVeF *= Nve;
		IndVeF += Indf*(NveMax*NfveMax*NfrefMax);

//printf("%d %d %d\n",Indf,Indsf,IndVeF);
//array_print_d(Nfve[0],Nve,&VeF[IndVeF],'R');
		XYZIn_fC = malloc(Nfve[Indf]*d * sizeof *XYZIn_fC); // free
		mm_CTN_d(Nfve[Indf],d,Nve,&VeF[IndVeF],XYZ_vC,XYZIn_fC);

		VOLUME  = FACET[0]->VOut;
printf("%d\n",VOLUME->indexg);
		Vf      = FACET[0]->VfOut;
		Indf    = Vf / NfrefMax; // face index (ToBeDeleted: Move to notation)
		Indsf   = Vf % NfrefMax; // sub face index (ToBeDeleted: Move to notation)

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nve   = ELEMENT->Nve;
		Nfve  = ELEMENT->Nfve;
		Nfref = ELEMENT->Nfref;
		VeF   = ELEMENT->VeF;

		NvnGs = ELEMENT->NvnGs[0];
		XYZ_vC = VOLUME->XYZ_vC;

		for (i = IndVeF = 0; i < Indsf; i++)
			IndVeF += Nfve[i];
		IndVeF *= Nve;
		IndVeF += Indf*(NveMax*NfveMax*NfrefMax);

//printf("%d %d %d\n",Indf,Indsf,IndVeF);
//array_print_d(Nfve[0],Nve,&VeF[IndVeF],'R');
		XYZOut_fC = malloc(Nfve[Indf]*d * sizeof *XYZOut_fC); // free
		mm_CTN_d(Nfve[Indf],d,Nve,&VeF[IndVeF],XYZ_vC,XYZOut_fC);

array_print_d(NvnGs,d,XYZ_vC,'C');
array_print_d(Nfve[Indf],d,XYZIn_fC,'C');

array_print_d(NvnGs,d,FACET[0]->VOut->XYZ_vC,'C');
array_print_d(Nfve[Indf],d,XYZOut_fC,'C');

		// Compute distance matrix
		Nfn  = Nfve[Indf];
		DXYZ = calloc(Nfn*Nfn , sizeof *DXYZ); // free
		compute_distance_matrix(Nfn,FACET[0]->BC,d,XYZIn_fC,XYZOut_fC,DXYZ);

array_print_d(Nfn,Nfn,DXYZ,'R');

		// Obtain the index of corresponding ordering between FACETs
		get_ordering_index(Nfn,d,DXYZ,&IndOrdInOut,&IndOrdOutIn);

printf("%d %d\n",IndOrdInOut,IndOrdOutIn);

/* 1) Here: Find number of FACET nodes for integration (or whatever other operator is needed)
 * 2) Here: Find node reordering
 * 3) In explicit_FACET_info: Check using VOLUME->What projected to the FACET. Modify What to be different at each
 *    point so that you can check this as this may not be the case based on the initial solution.
 */

// 1) Given PF  => NfnI[PF][IndClass]
// 2) THINK

		free(XYZIn_fC);
		free(XYZOut_fC);

		free(DXYZ);

exit(1);
	}
	free(FACET);
/*
	// Determine GFToV array
	GFToV = malloc(NGF*2 * sizeof *GFToV); // tbd
	GF_Nv = calloc(NGF   , sizeof *GF_Nv); // free

	for (gf = 0; gf < NGF; gf++)
		GFToV[gf*2+1] = NV;

	for (v = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		gf = VToGF[v*NfMax+f];

		GFToV[gf*2+GF_Nv[gf]] = v;
		GF_Nv[gf]++;
	}}
	free(GF_Nv);

array_print_ui(NGF,2,GFToV,'R');
exit(1);
*/

/*
for (FACET[0] = DB.FACET; FACET[0] != NULL; FACET[0] = FACET[0]->next) {
	printf("%d %d %d %c %d %d %d %d\n",
	       FACET[0]->indexg,FACET[0]->curved,FACET[0]->typeInt,FACET[0]->VIn->indexg,FACET[0]->VOut->indexg,
		   FACET[0]->VfIn,FACET[0]->VfOut,FACET[0]->BC);
}
exit(1);
*/

/*
for (i = 0, iMax = NVgrp; iMax--; i++) {
	for (VOLUME = Vgrp[i]; VOLUME != NULL; VOLUME = VOLUME->grpnext) {
		printf("%d %d %d %d\n",i,VOLUME->Eclass,VOLUME->P,VOLUME->curved);
	}
	printf("\t\t%p\n",Vgrp[i]);
}
*/
//exit(1);


	// Assign/Overwrite DB parameters
	DB.NV     = NVlocal;
	DB.NECgrp = NECgrp;

	DB.Vgrp = Vgrp;

//VOLUME = DB.VOLUME; while(VOLUME != NULL) printf("%d %d\n",VOLUME->type,VOLUME->curved), VOLUME = VOLUME->next;
}
