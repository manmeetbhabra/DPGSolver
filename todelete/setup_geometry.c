#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up geometry.
 *
 *	Comments:
 *		Loops are placed outside of the functions listed below (ToBeModified) as they are reused to update individual
 *		VOLUMEs after refinement:
 *			setup_ToBeCurved
 *			setup_geom_factors
 *		Potentially clean up this routine in the future by adding in an initialize operators routine and putting the set
 *		up for XYZ_S in a sub function. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
*/

void setup_geometry(void)
{
	// Initialize DB Parameters
	char         *MeshType = DB.MeshType,
	             *TestCase = DB.TestCase;
	unsigned int ExactGeom = DB.ExactGeom,
	             d         = DB.d,
	             *NE       = DB.NE,

	             Testing   = DB.Testing;

	int          PrintTesting = 0;

	// Standard datatypes
	unsigned int i, dim, P, vn,
	             Vs, NvnG, NvnGs, NvnGc, u1,
	             NIn, NOut, NIn_SF[3], NOut_SF[3], NCols, Diag[3], NOut_Total;
	double       *XYZ_vC, *XYZ_S, *XYZ,
	             *I_vGs_vGc, *Input_SF, *OP_SF[3];

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;
	struct S_FACET   *FACET;

	// silence
	NvnGs = 0; NvnGc = 0;
	XYZ_S = NULL;
	I_vGs_vGc = NULL;

	u1 = 1;
	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	// Modify vertex locations if exact geometry is known
	if (ExactGeom) {
		if(!DB.MPIrank) printf("    Modify vertex nodes if exact geometry is known\n");
		printf("Did not yet verify the implementation.\n");
		vertices_to_exact_geom();
	}

	// Set up global XYZ_(S)tart VOLUME coordinates at (S)tart (i.e. before curving)

	/* For the vectorized version of the code, set up a linked list looping through each type of element and each
	 * polynomial order; do this for both VOLUMEs and FACETs (eventually). Then this list can be used to sequentially
	 * performed vectorized operations on each type of element at once.
	 */

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		P      = VOLUME->P;
		XYZ_vC = VOLUME->XYZ_vC;

		ELEMENT = get_ELEMENT_type(VOLUME->type);
		if (!VOLUME->curved) {
			// If not curved, the P1 geometry representation suffices to fully specify the element geometry.
			if (VOLUME->Eclass == C_TP)
				NvnG = pow(ELEMENT->ELEMENTclass[0]->NvnGs[0],d);
			else if (VOLUME->Eclass == C_WEDGE)
				NvnG = (ELEMENT->ELEMENTclass[0]->NvnGs[0])*(ELEMENT->ELEMENTclass[1]->NvnGs[0]);
			else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR)
				NvnG = ELEMENT->NvnGs[0];
			else
				printf("Error: Unsupported element type setup_geom (NvnG).\n"), exit(1);

			VOLUME->NvnG = NvnG;

			XYZ_S = malloc(NvnG*d * sizeof *XYZ_S); // keep
			XYZ  = malloc(NvnG*d * sizeof *XYZ);  // keep
			VOLUME->XYZ_S = XYZ_S;

			for (dim = 0; dim < d; dim++) {
			for (vn = 0; vn < NvnG; vn++) {
				XYZ_S[dim*NvnG+vn] = XYZ_vC[dim*NvnG+vn];
			}}
		} else {
			if (VOLUME->Eclass == C_TP) {
				NvnGs = ELEMENT->ELEMENTclass[0]->NvnGs[0];
				NvnGc = ELEMENT->ELEMENTclass[0]->NvnGc[P];

				I_vGs_vGc = ELEMENT->ELEMENTclass[0]->I_vGs_vGc[P];

				Input_SF = XYZ_vC; // note multi column input

				NIn  = NvnGs;
				NOut = NvnGc;
				NOut_Total = 1;
				for (dim = 0; dim < 3; dim++) {
					if (dim < d) {
						NIn_SF[dim]  = NIn;
						NOut_SF[dim] = NOut;

						OP_SF[dim] = I_vGs_vGc;
						Diag[dim]  = 0;
					} else {
						NIn_SF[dim]  = 1;
						NOut_SF[dim] = 1;

						OP_SF[dim] = NULL;
						Diag[dim]  = 2;
					}
					NOut_Total *= NOut_SF[dim];
				}

				NCols    = d*1; // d coordinates * 1 element

				VOLUME->NvnG = NOut_Total;

				XYZ_S = malloc(NOut_Total*NCols * sizeof *XYZ_S); // keep
				XYZ  = malloc(NOut_Total*NCols * sizeof *XYZ);  // keep
				sf_apply_d(Input_SF,XYZ_S,NIn_SF,NOut_SF,NCols,OP_SF,Diag,d);
			} else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR) {
				NvnGs = ELEMENT->NvnGs[0];
				NvnGc = ELEMENT->NvnGc[P];
				I_vGs_vGc = ELEMENT->I_vGs_vGc[P];

				NCols = d*1; // d coordinates * 1 element

				VOLUME->NvnG = NvnGc;

				XYZ_S = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep
				XYZ  = malloc(NvnGc*NCols * sizeof *XYZ);  // keep

				mm_d(CblasColMajor,CblasTrans,CblasNoTrans,NvnGc,NCols,NvnGs,1.0,I_vGs_vGc,XYZ_vC,XYZ_S);
			} else if (VOLUME->Eclass == C_WEDGE) {
				Input_SF = XYZ_vC;

				NOut_Total = 1;
				for (dim = 0; dim < 3; dim++) {
					if (dim == 0 || dim == 2) {
						NIn_SF[dim]  = ELEMENT->ELEMENTclass[min(dim,u1)]->NvnGs[0];
						NOut_SF[dim] = ELEMENT->ELEMENTclass[min(dim,u1)]->NvnGc[P];

						OP_SF[dim] = ELEMENT->ELEMENTclass[min(dim,u1)]->I_vGs_vGc[P];
						Diag[dim]  = 0;
					} else {
						NIn_SF[dim]  = 1;
						NOut_SF[dim] = 1;

						OP_SF[dim] = NULL;
						Diag[dim]  = 2;
					}
					NOut_Total  *= NOut_SF[dim];
				}

				NCols = d*1; // d coordinates * 1 element

				VOLUME->NvnG = NOut_Total;

				XYZ_S = malloc(NOut_Total*NCols * sizeof *XYZ_S); // keep
				XYZ  = malloc(NOut_Total*NCols * sizeof *XYZ);  // keep
				sf_apply_d(Input_SF,XYZ_S,NIn_SF,NOut_SF,NCols,OP_SF,Diag,d);
			}
		}
		VOLUME->XYZ_S = XYZ_S;

//array_print_d(VOLUME->NvnG,d,VOLUME->XYZ_S,'C');
//exit(1);
	}

	if (Testing) {
		// Output straight coordinates to paraview
		output_to_paraview("ZTest_Geom_straight");
	}

	// Set up curved geometry nodes
	if (strstr(MeshType,"ToBeCurved") != NULL) {
		printf("    Set geometry of VOLUME nodes in ToBeCurved Mesh\n");
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
			setup_ToBeCurved(VOLUME);

/*
VOLUME = DB.VOLUME;
array_print_d(VOLUME->NvnG,d,VOLUME->XYZ_S,'C');
array_print_d(VOLUME->NvnG,d,VOLUME->XYZ,'C');
exit(1);
*/

	} else {
		printf("Add in support for MeshType != ToBeCurved");
		exit(1);
	}

	if (Testing) {
		// Output curved coordinates to paraview
		output_to_paraview("ZTest_Geom_curved");
	}

	printf("    Set up geometric factors\n");
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		setup_geom_factors(VOLUME);

	printf("    Set up normals\n");
	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next)
		setup_normals(FACET);
//exit(1);

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		free(VOLUME->C_vC);

	if (Testing) {
		// Output normals to paraview
		output_to_paraview("ZTest_Normals");
	}


	// Find node index ordering on each FACET
	// This will be done in setup_structures using XYZ_vC.

	/* WRITE A TEST ROUTINE TO MAKE SURE THAT THIS IS WORKING? (ToBeDeleted)
	 * Note: Ideally, this information can be given without relying on matching of physical points as there will be a
	 *       potential for non-conforming elements.
	 * For the "sum-factorization" on simplex elements, the nodes must be ordered according to their symmetries, which
	 * is unrelated to the vertex positions? Perhaps define the nodes only in one section (1/3) of the reference
	 * triangle with the associated multiplicity.
	 * Note: After having read through Hesthaven(2000) on the sum-factorization on triangles, it seemed intuitive to
	 *       attempt a similar extension for the TP case, using the symmetry about 0. This succeeded, resulting in an
	 *       asymptotic complexity reduction of 2. This does require re-ordering of both the nodes and basis functions
	 *       however => Do this. Also, while I thought I might be onto a new result, a very similar demonstration seems
	 *       to have been made in Solomonoff(1992)-A_Fast_Algorithm_for_Spectral_Differentiation, although he does not
	 *       seem to have considered the interpolation.
	 */


}