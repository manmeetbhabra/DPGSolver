// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_operators.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_OpCSR.h"

#include "bases.h"
#include "cubature.h"

#include "setup_operators_support.h"
#include "select_functions.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "sum_factorization.h"
#include "array_free.h"
#include "plotting_element_info.h"
#include "adaptation.h"
#include "memory_destructors.h"
#include "array_norm.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Set up operators to be used throughout the code.
 *
 *	Comments:
 *
 *		In future versions of the code, it would be beneficial to cast operators const after initial computation here
 *		such that they not be modified inadvertently in other functions.
 *
 *		Different operators are set up depending on which higher-dimensional element dependencies are present as well as
 *		whether adaptivity is desired:
 *			LINE  : Standard operators + additional SF operators if QUADs are present (standard assembly)
 *			TRI   : Standard operators + additional SF operators if WEDGEs are present (standard assembly)
 *			QUAD  : Standard operators + sparse versions of heavily used operators (TP asssembly)
 *			TET   : Standard operators (standard assembly)
 *			HEX   : Standard operators + sparse versions of heavily used operators (TP asssembly)
 *			WEDGE : Standard operators + sparse versions of heavily used operators (TP asssembly)
 *			PYR   : Standard operators (standard assembly)
 *
 *			This is done so that the only operators computed are those which are used. The additional SF operators are
 *			needed to interpolate from VOLUME to VOLUME nodes in components of the TP operator assembly for higher
 *			dimensional differentiation and FACE operators.
 *			It must be noted, however, that the data structure is such that the most general case can be handled (3D
 *			mixed meshes using hp adaptivity and sum factorized operators) which results in notably excessive operator
 *			indexing in the much simpler cases (such as in 1D). The indexing convention is as follows:
 *				OP_v*_(v/f)*[1][2][3] : (OP)erator from order [1] to order [2] with (v)olume/(f)ace refinement index
 *				                        [3] where potential (OP)erators and '*'s are defined in the notation below.
 *			If sum factorization has not yet broken-even for a certain order, the associated operators are freed here.
 *			(ToBeModified: DON'T FORGET TO DO THIS. Remove the sum factorized output_to_paraview function as well. This
 *			can be done when sum factorization is included in the solver).
 *
 *		When collocation is enabled (such that mass matrices are diagonal), the inverse weights are prepended to the
 *		most commonly used solver operators ([I/D][s/c]_[Weak/Strong]_[VV/FV]). This is especially beneficial as it
 *		results in the I[s/c]_Weak_VV being identity (instead of diagonal), reducing the required sum factorized
 *		operations when applying the derivative operators; the impact of this saving is significant for low order, high
 *		dimension simulations. Despite this advantage, complications arise from this treatment, with those currently
 *		known listed below:
 *			1) Using the same operator for construction of the Poisson system of equations leads to an asymmetric global
 *			   matrix, which must be subsequently corrected through multiplication by the diagonal weights to recover
 *			   the symmetry.
 *
 *		Standard (i.e. non sum-factorized) operators are used for all functions which are not performance critical in
 *		order to improve code readability.
 *		Ensure that operators for hp refinement are only stored when refinement is enabled (ToBeDeleted).
 *
 *		For the PYR element rst_vGs != E_rst_vV because of the TP structures of the TP nodes and the rotational symmetry
 *		ordering of the PYR nodes in layers of 't'. This means that special consideration must be made if attempting to
 *		use the rotational symmetry of rst_vGs PYR nodes.
 *
 *		The most intuitive method for obtaining projections between different h-refined spaces is to use the
 *		(Q)uadrature (M)irror (F)ilter coefficient matrices from the multiresolution framework. As the code discards the
 *		detail coefficients when projecting from fine to coarser spaces, the wavelet basis functions are not needed.
 *		Furthermore, the information embedded in the detail coefficients is simply an elegant representation of lack of
 *		smoothness.
 *		When imagining the use of a modal basis in the discussion of Kopriva(1996) in the context of h-adaptation
 *		related L2 projections, it is easily seen that the QMF coefficient matrices perform exactly such projections.
 *		See Alpert(2002) for their definitions.
 *
 *		Consider removing the "-ve" sign included in the I()_Weak_FV operator to avoid confusion. (ToBeModified)
 *
 *	Notation:
 *
 *	General:
 *
 *		Eclass : (E)lement (class)
 *		         Options: TP, SI, PYR, WEDGE
 *
 *		N(1)n(2)(3)[4][5] : (N)umber of (1) (n)odes of (2) type on elements which are (3) of order [4] belonging to
 *		                    element class [5]
 *		                    (1): (v)olume, (f)ace, (e)dge
 *		                    (2): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                    (3): (s)traight, (c)urved
 *		  Example: NfnIs[1][0] == (N)umber of (f)ace (n)odes for (I)ntegration on (straight) P1 faces of Eclass 0.
 *
 *		rst_(1)(2)(3) : coordinates on the reference element (rst) of (1) nodes of (2) type which are (3).
 *		                (1): (v)olume, (f)ace, (e)dge
 *		                (2): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution, (V)ertex
 *		                (3): (s)traight, (c)urved (optional)
 *
 *		(Grad)Chi(Ref)(1)(2)_(3)(4)(5)[6] : (Grad)ient (optional) (Ref)erence (optional) Basis functions (Chi) of type
 *		                                    (1) which are (2) evaluated at (3) nodes of (4) type which are (5) of order [6]
 *		                                    (1/4): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                                    (2/5): (s)traight, (c)urved
 *		                                    (3): (v)olume, (f)ace, (e)dge
 *
 *		_(1)(2) subscripts for operators indicate that this is a (1) operator, operating on VOLUME nodes and
 *		transferring to (2) nodes.
 *			Options (1)(2) : VV, FF
 *
 *	setup_ELEMENT_operators/setup_TP_operators:
 *		w_(1)(2)(3) : Cubature (w)eights (1) nodes of (2) type which are (3)
 *		             (1): (v)olume, (f)ace
 *		             (2): (I)ntegration
 *		             (3): (s)traight, (c)urved
 *		I(1)(2) : (I)dentity matrix of type (1) which is (2)
 *		          (1): (G)eometry, (C)ofactor, (S)olution
 *		          (2): (s)traight, (c)urved, (2): order 2
 *		T(1)(2)(3) : (T)ransformation matrix of type (1) which is (2). (Zwanenburg(2016): eq. 2.14)
 *		             Optional argument (3) for additional functionality.
 *		          (1): (G)eometry, (C)ofactor, (S)olution
 *		          (2): (s)traight, (c)urved
 *		          (3): _(v)olume(B)ezier
 *		I_(1)(2)(3)_(4)(5)(6) : (I)nterpolation operator from (1) nodes of type (2) which are (3) to (4) nodes of type
 *		                        (5) which are (6)
 *		                        (1/4): (v)olume, (f)ace, (e)dge
 *		                        (2/5): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                        (3/6): (s)traight, (c)urved
 *		  Note: Interpolation operators to FACEs may have up to 4 levels of dereferencing *I_x_x[1][2][3]:
 *		        [0] : The operator is a pointer to  matrix
 *		        [1] : Order of the element from which you are interpolating
 *		        [2] : Order of the element to which you are interpolating
 *		        [3] : Index of operator (e.g. there is a need for many operators interpolating to different parts of the
 *		                                 VOLUME FACEs)
 *		Ihat_(1)(2)(3)_(4)(5)(6) : Analogous to I (above) but where interpolation is of coefficients instead of values.
 *		  Note: values and coefficients are the same for the nodal scheme, but not for modal.
 *		D_(1)(2)(3)_(4)(5)(6) : (D)ifferentiation + interpolation operator from (1) nodes of type (2) which are (3) to
 *		                        (4) nodes of type (5) which are (6)
 *		                        (1/4): (v)olume, (f)ace, (e)dge
 *		                        (2/5): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                        (3/6): (s)traight, (c)urved
 *		  Note: Differentiation + interpolation operators may have up to 5 levels of dereferencing *D_x_x[1][2][3][4]:
 *		        [0] : The operator is a pointer to  matrix
 *		        [1] : Order of the element from which you are interpolating with differentiation
 *		        [2] : Order of the element to which you are interpolating with differentiation
 *		        [3] : Index of operator
 *		        [4] : Index of dimension of the operator
 *
 *		*_sp : Standard (*) operator stored in (sp)arse format ((C)ompressed (S)parse (R)ow), only used for TP ELEMENTs.
 *
 *	setup_L2_projection_operators:
 *		L2hat_vS_vS : L2 projection operator from coefficients of the VOLUMEs on level l+1 to the VOLUME on level l.
 *
 *
 *	Returned operators : [range used depending on Adapt] (ToBeModified: Add missing operators)
 *		- For maximum array size, see memory_constructors.c.
 *		- If no range is indicated, the range is given by that of the previous operator in the list.
 *
 *		- (*) indicates a special case for SF operators.
 *		- [^] is a placeholder.
 *		- Further notation simplifications for table below:
 *			*    : 0:NVREFSFMAX
 *			**   : 0:NREFMAX
 *			^    : 0:NFMAX
 *			^^   : 0:NFREFMAX*NFMAX
 *
 *			rP   : (r)ange P  == 0:PMax
 *			P(P) : P as a function of P (i.e. for a range in field 1 of [0:PMax] -> [0][0], [1][1], ..., [P][P])
 *			rPb  : (r)ange Pb (always a function of P) == PbMin:PbMax (PbMin = min(P-1,0), PbMax = max(P+1,PMax))
 *			rd   : (r)ange (d)imension == 0:d-1
 *
 *		                   Adapt
 *		Returned Op (SF) | ADAPT_0         ADAPT_P             ADAPT_H         ADAPT_HP
 *		                 |
 *		NvnGs            | [1]             [1]                 [1]             [1]
 *		NvnGc            | [P]             [rP]                [P]             [rP]
 *		NvnCs            |
 *		NvnCc            |
 *		NvnIs            |
 *		NvnIc            |
 *		NfnIs            | [P][0:1]        [rP][0:1]           [P][2]          [rP][0:1]
 *		NfnIc            |
 *		NenGc            | [P]             [rP]                [P]             [rP]
 *		                 |
 *		w_vIs            | [P]             [rP]                [P]             [rP]
 *		w_vIc            |
 *
 *		ChiS_vP          | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		ChiS_vS      (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		ChiS_vIs     (*) |
 *		ChiS_vIc     (*) |
 *		ChiInvS_vS       | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		                 |
 *		IGc              | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		ICs              |
 *		ICc              |
 *		                 |
 *		I_vGs_vP         | [1][P][0]       [1][rP][0]          [1][P][0]       [1][rP][0]
 *		I_vGs_vGs        | NOT_USED        NOT_USED            [1][1:2][**]    [1][1:2][**]
 *		I_vGs_vG2        | NOT_USED        NOT_USED            NOT_USED        [1][2][0]
 *		I_vGs_vGc        | [1][P][0]       [1][rP][0]          [1][P][0]       [1][rP][0]
 *		I_vGs_vCs        |
 *		I_vGs_vS     (*) | [1][P][0]       [1][rPb][0]         [1][P][*]       [1][rP][*]
 *		I_vGs_vIs    (*) |
 *		I_vGs_vIc    (*) |
 *		I_vGc_vP         | [P][P][0]       [P][P(P)][0]        [P][P][0]       [P][P(P)][0]
 *		I_vGc_vCc        |
 *		I_vGc_vS     (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		I_vGc_vIs    (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		I_vGc_vIc    (*) |
 *		I_vCs_vS     (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		I_vCs_vIs    (*) |
 *		I_vCs_vIc    (*) |
 *		I_vCc_vS     (*) |
 *		I_vCc_vIs    (*) |
 *		I_vCc_vIc    (*) |
 *		                 |
 *		Ihat_vS_vS       | NOT_USED        [rP][rPb][0]        [P][P][*]       [rP][rP][0] & [P][P][*]
 *		                 |
 *		D_vGs_vCs        | [1][P][0][rd]   [1][rP][0][rd]      [1][P][0][rd]   [1][rP][0][rd]
 *		D_vGs_vIs        |
 *		D_vGc_vCc        | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		D_vGc_vIc        |
 *		D_vCs_vCs        |
 *		D_vCc_vCc        |
 *		                 |
 *		ChiS_fS          | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rP][^^]
 *		ChiS_fIs         |
 *		ChiS_fIc         |
 *		                 |
 *		I_vGs_fS         | [1][P][^]       [P][rP][^]          [1][P][^^]      [1][rP][^^]
 *		I_vGs_fIs        |
 *		I_vGs_fIc        |
 *		I_vGc_fS         | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rP][^^]
 *		I_vGc_fIs        |
 *		I_vGc_fIc        |
 *		I_vCs_fIs        |
 *		I_vCs_fIc        |
 *		I_vCc_fIs        |
 *		I_vCc_fIc        |
 *		I_fGc_vGc        | [P][P][0]       [P][P(P)][0]        [P][P][0]       [P][P(P)][0]
 *		                 |
 *		Is_Weak_VV   (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		Ic_Weak_VV   (*) |
 *		Ds_Weak_VV       | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		Dc_Weak_VV       |
 *		                 |
 *		Is_Weak_FV       | [P][P][^]       [rP][rPb][^]       [P][P][^^]       [rP][rP][^^]
 *		Ic_Weak_FV       |
 *		                 |
 *		ChiS_fIs_sp      | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rP][^^]
 *		ChiS_fIc_sp      |
 *		Ds_Weak_VV_sp    | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		Dc_Weak_VV_sp    |
 *		Is_Weak_FV_sp    | [P][P][^]       [rP][rPb][^]       [P][P][^^]       [rP][rP][^^]
 *		Ic_Weak_FV_sp    |
 *
 *
 *	References:
 *		Zwanenburg(2016)_Equivalence between the Energy Stable Flux Reconstruction and Filtered Discontinuous Galerkin Schemes
 *		Karniadakis(1999)_Spectral hp Element Methods for CFD
 *		Hesthaven(2000)_Stable Spectral Methods on Tetrahedral Elements
 *		Stiller(2008)_Factorization Techniques for Nodal Spectral Elements in Curved Domains
 *		Alpert(2002)-Adaptive_Solution_of_Partial_Differential_Equations_in_Multiwavelet_Bases
 *		Kopriva(1996)-A_Conservative_Staggered-Grid_Chebyshev_Multidomain_Method_for_Compressible_Flows_II._A_Semi-Structured_Method
 *		Szabo(1991)-Finite_Element_Analysis
 *
 *		ToBeDeleted: Potentially add in MATH578 report instead of Alpert(2002) for the QMF discussion.
 */

static void setup_ELEMENT_operator_dependencies(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnIs, *NvnIc, *NvnS;
	double       ****ChiInvS_vS;

	// Initialize DB Parameters
	unsigned int PGs     = DB.PGs,
	             **PIvs  = DB.PIvs,
	             **PIvc  = DB.PIvc;

	char         *BasisType     = DB.BasisType,
	             **NodeTypeG    = DB.NodeTypeG,
	             ***NodeTypeIvs = DB.NodeTypeIvs,
	             ***NodeTypeIvc = DB.NodeTypeIvc,
	             ***NodeTypeS   = DB.NodeTypeS;

	// Standard datatypes
	unsigned int dE, Eclass, Nbf;
	double       *rst_vS, *IS, *ChiRefS_vS, *ChiS_vS;

	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	ChiS_vS = NULL;

	ELEMENT = get_ELEMENT_type(EType);
	Eclass  = get_Eclass(EType);

	dE = ELEMENT->d;

	select_functions(&basis,&grad_basis,&cubature,EType);

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NvnS  = ELEMENT->NvnS;

	ChiInvS_vS = ELEMENT->ChiInvS_vS;

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,PGs,cubature); // free
	free(CUBDATA->rst);
	NvnGs[1] = CUBDATA->Nn;

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	for (size_t P = PSMin; P <= PSMax; P++) {
		size_t PbMin, PbMax;
		get_Pb_range(P,&PbMin,&PbMax);
		for (size_t Pb = PbMax; Pb >= P; Pb--) {
			set_cubdata(CUBDATA,false,false,NodeTypeS[Pb][Eclass],dE,Pb,cubature); // free
			set_from_cubdata(CUBDATA,&NvnS[Pb],NULL,&rst_vS,NULL,NULL);

			IS         = identity_d(NvnS[Pb]);  // free
			ChiRefS_vS = basis(Pb,rst_vS,NvnS[Pb],&Nbf,dE); // free
			free(rst_vS);

			if (strstr(BasisType,"Modal")) {
				ChiS_vS = ChiRefS_vS;
			} else if (strstr(BasisType,"Nodal")) {
				ChiS_vS = IS;
			}

			if (ChiInvS_vS[Pb][Pb][0] == NULL)
				ChiInvS_vS[Pb][Pb][0] = inverse_d(NvnS[Pb],NvnS[Pb],ChiS_vS,IS); // keep

			free(IS);
			free(ChiRefS_vS);
			if (Pb == P) {
				// Avoid negative Pb for P = 0
				break;
			}
		}

		for (size_t Pb = PbMin; Pb <= PbMax; Pb++) {
			set_cubdata(CUBDATA,false,false,NodeTypeIvs[Pb][Eclass],dE,PIvs[Pb][Eclass],cubature); // free
			free(CUBDATA->rst);
			NvnIs[Pb] = CUBDATA->Nn;

			set_cubdata(CUBDATA,false,false,NodeTypeIvc[Pb][Eclass],dE,PIvc[Pb][Eclass],cubature); // free
			free(CUBDATA->rst);
			NvnIc[Pb] = CUBDATA->Nn;
		}
	}

	free(CUBDATA);
}

static void setup_ELEMENT_operators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnG2, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS, **NfnG2, **NfnGc, **NfnS, **NfnIs,
	             **NfnIc, *NenG2, *NenGc;
	double       **w_vIs, **w_vIc, ***w_fIs, ***w_fIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS, ****ChiInvGs_vGs, ****ChiBezInvS_vS,
	             ****TGs, ****TS, ****TS_vB, ****TInvS_vB,
	             ****IG2, ****IGc, ****ICs, ****ICc,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vG2, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vIs, ****I_vGs_vIc,
	             ****I_vGs_vS,
	             ****I_vGc_vP,                ****I_vGc_vCc, ****I_vGc_vIs, ****I_vGc_vIc, ****I_vGc_vS,
	             ****I_vCs_vS, ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vS, ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             *****GradChiS_vS, *****GradChiS_vIs, *****GradChiS_vIc,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
	             *****GradChiS_fIs, *****GradChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vG2_fG2, ****I_vG2_eG2,
	             ****I_vGc_fGc, ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc, ****I_vGc_eGc,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
	             *****D_vGs_fIs, *****D_vGs_fIc,
	             *****D_vGc_fIs, *****D_vGc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV, ****Is_Strong_VV, ****Ic_Strong_VV,
	             ****Is_Weak_FV, ****Ic_Weak_FV,
	             *****Ds_Weak_VV, *****Dc_Weak_VV,*****Ds_Strong_VV, *****Dc_Strong_VV;

	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               ****Is_Weak_FV_sp, ****Ic_Weak_FV_sp;

	// Initialize DB Parameters
	unsigned int PMax         = DB.PMax,
	             PGlobal      = DB.PGlobal,
	             NP           = DB.NP,
	             PGs          = DB.PGs,
	             *PGc         = DB.PGc,
	             **PCs        = DB.PCs,
	             **PCc        = DB.PCc,
	             **PIvs       = DB.PIvs,
	             **PIvc       = DB.PIvc,
	             Collocated   = DB.Collocated,
	             EFE          = DB.EFE,
	             *VFPartUnity = DB.VFPartUnity;

	char         *BasisType     = DB.BasisType,
	             **NodeTypeG    = DB.NodeTypeG,
	             ***NodeTypeIvs = DB.NodeTypeIvs,
	             ***NodeTypeIvc = DB.NodeTypeIvc,
	             ***NodeTypeS   = DB.NodeTypeS;

	// Standard datatypes
	unsigned int i, iMax, dim, dE, P, e, f, vh, fh, Ve, Vf, IndFType, fhMax,
	             Nve, Nve_h, Ne, Nf, Nbf, Eclass, NFTypes, Nvref, vref, vrefSF, NvrefSF, *Nfref, *ones_Nf, NEhref,
	             Indh, NvnP, u1 = 1,
	             BE_Nve, BF_Nve[2], Neve, *Nfve, *Nvve, *EType_h,
	             dummy_ui, *dummyPtr_ui[3];
	double       *E_rst_vV, *rst_vV, **VeE, **VeF, **VeV,
	             *rst_vP, *rst_vG2, *rst_vGc, *rst_vCs, *rst_vCc, **rst_vIs, **rst_vIc, **rst_vS,
	             *rst_fG2, *rst_fGc, *rst_fS, *rst_fIs, *rst_fIc,
	             *rst_eG2, *rst_eGc,
	             *wInv_vIs, *wInv_vIc,
	             *diag_w_vIs, *diag_w_vIc, *diag_wInv_vIs, *diag_wInv_vIc,
	             *diag_w_fIs, *diag_w_fIc,
	             *IGs, *IS,
	             *TG2, *TGc, *TCs, *TCc,
	             *ChiRefGc_vGc, *ChiRefCs_vCs, *ChiRefCc_vCc, *ChiRefS_vS, *ChiBezS_vS,
	             *ChiGc_vGc,    *ChiCs_vCs,    *ChiCc_vCc,
	             *ChiRefInvGs_vGs, *ChiRefInvG2_vG2, *ChiRefInvGc_vGc, *ChiRefInvCs_vCs, *ChiRefInvCc_vCc, *ChiRefInvS_vS,
	                               *ChiInvG2_vG2,    *ChiInvGc_vGc,    *ChiInvCs_vCs,    *ChiInvCc_vCc,
	             *ChiRefGs_vP, *ChiRefGs_vGs, *ChiRefGs_vG2, *ChiRefGs_vGc, *ChiRefGs_vCs, *ChiRefGs_vIs, *ChiRefGs_vIc, *ChiRefGs_vS,
	             *ChiRefG2_vG2,
	             *ChiRefGs_fS, *ChiRefGs_fIs, *ChiRefGs_fIc,
	             *ChiRefGc_vP,                *ChiRefGc_vCc, *ChiRefGc_vIs, *ChiRefGc_vIc, *ChiRefGc_vS,
	             *ChiRefG2_fG2,
	             *ChiRefGc_fGc, *ChiRefGc_fS, *ChiRefGc_fIs, *ChiRefGc_fIc, *ChiRefG2_eG2, *ChiRefGc_eGc,
	             *ChiRefCs_vS, *ChiRefCs_vIs, *ChiRefCs_vIc, *ChiRefCs_fS, *ChiRefCs_fIs, *ChiRefCs_fIc,
	             *ChiRefCc_vS, *ChiRefCc_vIs, *ChiRefCc_vIc, *ChiRefCc_fS, *ChiRefCc_fIs, *ChiRefCc_fIc,
	             *ChiRefS_vP, *ChiRefS_vIs, *ChiRefS_vIc, *ChiRefS_fS, *ChiRefS_fIs, *ChiRefS_fIc,
	             *ChiGs_vP, *ChiGs_vGs, *ChiGs_vG2, *ChiGs_vGc, *ChiGs_vCs, *ChiGs_vIs, *ChiGs_vIc, *ChiGs_vS,
	             *ChiG2_vG2,
	             *ChiGs_fS, *ChiGs_fIs, *ChiGs_fIc,
	             *ChiGc_vP,             *ChiGc_vCc, *ChiGc_vIs, *ChiGc_vIc, *ChiGc_vS,
	             *ChiG2_fG2,
	             *ChiGc_fGc, *ChiGc_fS, *ChiGc_fIs, *ChiGc_fIc, *ChiG2_eG2, *ChiGc_eGc,
	             *ChiCs_vS, *ChiCs_vIs, *ChiCs_vIc, *ChiCs_fS, *ChiCs_fIs, *ChiCs_fIc,
	             *ChiCc_vS, *ChiCc_vIs, *ChiCc_vIc, *ChiCc_fS, *ChiCc_fIs, *ChiCc_fIc,
	             **GradChiRefGs_vCs, **GradChiRefGs_vIs,
	             **GradChiRefGc_vCc, **GradChiRefGc_vIc,
	             **GradChiRefCs_vCs, **GradChiRefCc_vCc,
	             **GradChiRefS_vS, **GradChiRefS_vIs,  **GradChiRefS_vIc,
	             **GradChiRefGs_fIs, **GradChiRefGs_fIc,
	             **GradChiRefGc_fIs, **GradChiRefGc_fIc,
	             **GradChiRefS_fIs,  **GradChiRefS_fIc,
	             **GradChiGs_vCs, **GradChiGs_vIs,
	             **GradChiGc_vCc, **GradChiGc_vIc,
	             **GradChiCs_vCs, **GradChiCc_vCc,
	             **GradChiGs_fIs, **GradChiGs_fIc,
	             **GradChiGc_fIs, **GradChiGc_fIc,
	             *dummyPtr_d;

	struct BCoords {
		double **G2, **Gc, **Is, **Ic, **S;
	} *BCoords_E, *BCoords_F[2], **BCoords_V;

	struct S_BCOORDS *BCoords_dEm2, *BCoords_dEm1[2];
	struct S_ELEMENT *ELEMENT, *ELEMENT_h;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis, basis_Bezier;
	grad_basis_tdef grad_basis;

	// silence
	ChiGs_vGs = ChiG2_vG2 = ChiGc_vGc = NULL;
	ChiCs_vCs = NULL; ChiCc_vCc = NULL;
	ChiRefGs_vIs = ChiRefGs_vIc = NULL;
	BE_Nve = 0; BCoords_E = NULL; BCoords_dEm2 = NULL;

	ELEMENT   = get_ELEMENT_type(EType);
	ELEMENT_h = ELEMENT;

	// No need to consider the second Eclass as WEDGE basis functions will be built through a combination of lower
	// dimensional operators.
	Eclass = get_Eclass(EType);
	setup_ELEMENT_VeV(EType);

	dE      = ELEMENT->d;
	Nve     = ELEMENT->Nve;
	Ne      = ELEMENT->Ne;
	Nf      = ELEMENT->Nf;
	Nfref   = ELEMENT->Nfref;
	Neve    = ELEMENT->Neve;
	Nfve    = ELEMENT->Nfve;
	VeE     = ELEMENT->VeE;
	VeF     = ELEMENT->VeF;
	VeV     = ELEMENT->VeV;
	Nvve    = ELEMENT->Nvve;
	Nvref   = ELEMENT->Nvref;
	NvrefSF = ELEMENT->NvrefSF;
	NEhref  = ELEMENT->NEhref;
	EType_h = ELEMENT->type_h;

	select_functions(&basis,&grad_basis,&cubature,EType);
	if (EType != TET && EType != PYR) {
		select_functions_basis_Bezier(&basis_Bezier,EType);
	}

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnG2 = ELEMENT->NvnG2;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NvnS  = ELEMENT->NvnS;
	NfnG2 = ELEMENT->NfnG2;
	NfnGc = ELEMENT->NfnGc;
	NfnS  = ELEMENT->NfnS;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;
	NenG2 = ELEMENT->NenG2;
	NenGc = ELEMENT->NenGc;

	w_fIs = ELEMENT->w_fIs;
	w_fIc = ELEMENT->w_fIc;

	w_vIs = ELEMENT->w_vIs;
	w_vIc = ELEMENT->w_vIc;

	ChiS_vP      = ELEMENT->ChiS_vP;
	ChiS_vS      = ELEMENT->ChiS_vS;
	ChiS_vIs     = ELEMENT->ChiS_vIs;
	ChiS_vIc     = ELEMENT->ChiS_vIc;
	ChiInvS_vS   = ELEMENT->ChiInvS_vS;
	ChiInvGs_vGs = ELEMENT->ChiInvGs_vGs;

	ChiBezInvS_vS = ELEMENT->ChiBezInvS_vS;

	IG2 = ELEMENT->IG2;
	IGc = ELEMENT->IGc;
	ICs = ELEMENT->ICs;
	ICc = ELEMENT->ICc;

	TGs   = ELEMENT->TGs;
	TS    = ELEMENT->TS;
	TS_vB = ELEMENT->TS_vB;

	TInvS_vB = ELEMENT->TInvS_vB;

	GradChiS_vS  = ELEMENT->GradChiS_vS;
	GradChiS_vIs = ELEMENT->GradChiS_vIs;
	GradChiS_vIc = ELEMENT->GradChiS_vIc;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGs = ELEMENT->I_vGs_vGs;
	I_vGs_vG2 = ELEMENT->I_vGs_vG2;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vIs = ELEMENT->I_vGs_vIs;
	I_vGs_vIc = ELEMENT->I_vGs_vIc;
	I_vGs_vS  = ELEMENT->I_vGs_vS;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vIs = ELEMENT->I_vGc_vIs;
	I_vGc_vIc = ELEMENT->I_vGc_vIc;
	I_vGc_vS  = ELEMENT->I_vGc_vS;
	I_vCs_vS  = ELEMENT->I_vCs_vS;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vS  = ELEMENT->I_vCc_vS;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	Ihat_vS_vS  = ELEMENT->Ihat_vS_vS;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	ChiS_fS     = ELEMENT->ChiS_fS;
	ChiS_fIs    = ELEMENT->ChiS_fIs;
	ChiS_fIc    = ELEMENT->ChiS_fIc;
	ChiS_fIs_sp = ELEMENT->ChiS_fIs_sp;
	ChiS_fIc_sp = ELEMENT->ChiS_fIc_sp;

	GradChiS_fIs = ELEMENT->GradChiS_fIs;
	GradChiS_fIc = ELEMENT->GradChiS_fIc;

	I_vGs_fS  = ELEMENT->I_vGs_fS;
	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vG2_fG2 = ELEMENT->I_vG2_fG2;
	I_vGc_fGc = ELEMENT->I_vGc_fGc;
	I_vGc_fS  = ELEMENT->I_vGc_fS;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vGc_eGc = ELEMENT->I_vGc_eGc;
	I_vG2_eG2 = ELEMENT->I_vG2_eG2;
	I_vCs_fS  = ELEMENT->I_vCs_fS;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fS  = ELEMENT->I_vCc_fS;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	D_vGs_fIs = ELEMENT->D_vGs_fIs;
	D_vGs_fIc = ELEMENT->D_vGs_fIc;
	D_vGc_fIs = ELEMENT->D_vGc_fIs;
	D_vGc_fIc = ELEMENT->D_vGc_fIc;

	Is_Weak_VV    = ELEMENT->Is_Weak_VV;
	Ic_Weak_VV    = ELEMENT->Ic_Weak_VV;
	Is_Strong_VV  = ELEMENT->Is_Strong_VV;
	Ic_Strong_VV  = ELEMENT->Ic_Strong_VV;
	Is_Weak_FV    = ELEMENT->Is_Weak_FV;
	Ic_Weak_FV    = ELEMENT->Ic_Weak_FV;
	Is_Weak_FV_sp = ELEMENT->Is_Weak_FV_sp;
	Ic_Weak_FV_sp = ELEMENT->Ic_Weak_FV_sp;

	Ds_Weak_VV = ELEMENT->Ds_Weak_VV;
	Dc_Weak_VV = ELEMENT->Dc_Weak_VV;
	Ds_Strong_VV = ELEMENT->Ds_Strong_VV;
	Dc_Strong_VV = ELEMENT->Dc_Strong_VV;

	// Allocate memory for arrays with multiple levels of dereferencing
	rst_vS        = malloc(NVREFMAX * sizeof *rst_vS);  // free
	rst_vIs       = malloc(NVREFMAX * sizeof *rst_vIs); // free
	rst_vIc       = malloc(NVREFMAX * sizeof *rst_vIc); // free

	GradChiGs_vCs = malloc(dE * sizeof *GradChiGs_vCs); // free
	GradChiGs_vIs = malloc(dE * sizeof *GradChiGs_vIs); // free
	GradChiGc_vCc = malloc(dE * sizeof *GradChiGc_vCc); // free
	GradChiGc_vIc = malloc(dE * sizeof *GradChiGc_vIc); // free
	GradChiCs_vCs = malloc(dE * sizeof *GradChiCs_vCs); // free
	GradChiCc_vCc = malloc(dE * sizeof *GradChiCc_vCc); // free

	GradChiGs_fIs = malloc(dE * sizeof *GradChiGs_fIs); // free
	GradChiGs_fIc = malloc(dE * sizeof *GradChiGs_fIc); // free
	GradChiGc_fIs = malloc(dE * sizeof *GradChiGc_fIs); // free
	GradChiGc_fIc = malloc(dE * sizeof *GradChiGc_fIc); // free

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	// VOLUME Nodes (Order Independent)
	set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,PGs,cubature); // free
	free(CUBDATA->rst);
	NvnGs[1] = CUBDATA->Nn;

	// Use E_rst_vV instead of rst_vGs (!= for PYR ELEMENTs)
	E_rst_vV = get_rst_vV(ELEMENT); // free
	rst_vV = malloc((Nve+1)*dE * sizeof *rst_vV); // free (+1 for h-refined TET -> PYR) ToBeModified (remove +1)

	// Preliminary Operators
	IGs = identity_d(NvnGs[1]); // free

	ChiRefGs_vGs = basis(PGs,E_rst_vV,NvnGs[1],&Nbf,dE); // free

	if (strstr(BasisType,"Modal"))
		ChiGs_vGs = ChiRefGs_vGs;
	else if (strstr(BasisType,"Nodal"))
		ChiGs_vGs = IGs;

	ChiRefInvGs_vGs       = inverse_d(NvnGs[1],NvnGs[1],ChiRefGs_vGs,IGs); // free
	ChiInvGs_vGs[1][1][0] = inverse_d(NvnGs[1],NvnGs[1],ChiGs_vGs,IGs);    // keep

	TGs[1][1][0]          = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGs[1],NvnGs[1],NvnGs[1],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // keep

	free(IGs);
	free(ChiRefGs_vGs);
	free(ChiRefInvGs_vGs);


	// Vertex blending operators
	P = 2;
	set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,P,cubature); // free
	rst_vG2 = CUBDATA->rst;
	NvnG2[P] = CUBDATA->Nn;

	ChiRefGs_vG2 = basis(PGs,rst_vG2,NvnG2[P],&Nbf,dE); // free
	ChiGs_vG2    = mm_Alloc_d(CBRM,CBNT,CBNT,NvnG2[P],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vG2,TGs[1][1][0]); // free

	I_vGs_vG2[1][P][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnG2[P],NvnGs[1],NvnGs[1],1.0,ChiGs_vG2,ChiInvGs_vGs[1][1][0]); // keep

	IG2[P][P][0] = identity_d(NvnG2[P]); // keep

	ChiRefG2_vG2 = basis(P,rst_vG2,NvnG2[P],&Nbf,dE); // free
	free(rst_vG2);

	if (strstr(BasisType,"Modal"))
		ChiG2_vG2 = ChiRefG2_vG2;
	else if (strstr(BasisType,"Nodal"))
		ChiG2_vG2 = IG2[P][P][0];

	ChiRefInvG2_vG2 = inverse_d(NvnG2[P],NvnG2[P],ChiRefG2_vG2,IG2[P][P][0]); // free
	ChiInvG2_vG2    = inverse_d(NvnG2[P],NvnG2[P],ChiG2_vG2,IG2[P][P][0]);    // free
	TG2             = mm_Alloc_d(CBRM,CBNT,CBNT,NvnG2[P],NvnG2[P],NvnG2[P],1.0,ChiRefInvG2_vG2,ChiG2_vG2); // free

	free(ChiRefGs_vG2);
	free(ChiGs_vG2);
	free(ChiRefG2_vG2);
	free(ChiRefInvG2_vG2);


	// Get Barycentric coordinates for lower dimensional ELEMENTs

	// EDGEs
	if (dE == DMAX) {
		BCoords_E = malloc(sizeof *BCoords_E); // free
		BCoords_dEm2 = get_BCoords_dEm2(ELEMENT);

		BE_Nve        = BCoords_dEm2->Nve;
		BCoords_E->G2 = BCoords_dEm2->BCoords_G2;
		BCoords_E->Gc = BCoords_dEm2->BCoords_Gc;
		for (P = 0; P <= PMax; P++) {
			NenG2[P] = BCoords_dEm2->NenG2[P];
			NenGc[P] = BCoords_dEm2->NenGc[P];
		}
		free(BCoords_dEm2->NenG2);
		free(BCoords_dEm2->NenGc);
	}

	// FACEs
	NFTypes = 1;
	BCoords_dEm1[0] = get_BCoords_dEm1(ELEMENT,0); // keep/free
	if (EType == WEDGE || EType == PYR) {
		NFTypes = 2;
		BCoords_dEm1[1] = get_BCoords_dEm1(ELEMENT,1);
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++)
		BCoords_F[IndFType] = malloc(sizeof *BCoords_F[IndFType]); // free


	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		BF_Nve[IndFType]        = BCoords_dEm1[IndFType]->Nve;
		BCoords_F[IndFType]->G2 = BCoords_dEm1[IndFType]->BCoords_G2;
		BCoords_F[IndFType]->Gc = BCoords_dEm1[IndFType]->BCoords_Gc;
		BCoords_F[IndFType]->S  = BCoords_dEm1[IndFType]->BCoords_S;
		BCoords_F[IndFType]->Is = BCoords_dEm1[IndFType]->BCoords_Is;
		BCoords_F[IndFType]->Ic = BCoords_dEm1[IndFType]->BCoords_Ic;

		// Store Nfn* in a manner consistent with Nvn*
		for (P = 0; P <= PMax; P++) {
			NfnG2[P][IndFType] = BCoords_dEm1[IndFType]->NfnG2[P];
			NfnGc[P][IndFType] = BCoords_dEm1[IndFType]->NfnGc[P];
			NfnS[P][IndFType]  = BCoords_dEm1[IndFType]->NfnS[P];
			NfnIs[P][IndFType] = BCoords_dEm1[IndFType]->NfnIs[P];
			NfnIc[P][IndFType] = BCoords_dEm1[IndFType]->NfnIc[P];
			w_fIs[P][IndFType] = BCoords_dEm1[IndFType]->w_fIs[P];
			w_fIc[P][IndFType] = BCoords_dEm1[IndFType]->w_fIc[P];
		}
		free(BCoords_dEm1[IndFType]->NfnG2);
		free(BCoords_dEm1[IndFType]->NfnGc);
		free(BCoords_dEm1[IndFType]->NfnS);
		free(BCoords_dEm1[IndFType]->NfnIs);
		free(BCoords_dEm1[IndFType]->NfnIc);
		free(BCoords_dEm1[IndFType]->w_fIs);
		free(BCoords_dEm1[IndFType]->w_fIc);
	}

	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	BCoords_V = malloc(NEhref * sizeof *BCoords_V); // free
	for (i = 0, iMax = NEhref; i < iMax; i++) {
		BCoords_V[i]     = malloc(     sizeof *BCoords_V[i]);       // free
		BCoords_V[i]->S  = calloc(NP , sizeof *(BCoords_V[i]->S));  // free
		BCoords_V[i]->Is = calloc(NP , sizeof *(BCoords_V[i]->Is)); // free
		BCoords_V[i]->Ic = calloc(NP , sizeof *(BCoords_V[i]->Ic)); // free
	}

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	for (size_t P = PSMin; P <= PSMax; P++) {
		size_t PbMin, PbMax;
		get_Pb_range(P,&PbMin,&PbMax);

		// Build preliminary operators needed for P adaptation.
		IS = ChiRefS_vS = NULL;

		for (size_t Pb = PbMax; Pb >= P; Pb--) {
			set_cubdata(CUBDATA,false,false,NodeTypeS[Pb][Eclass],dE,Pb,cubature); // free
			set_from_cubdata(CUBDATA,&NvnS[Pb],NULL,&rst_vS[0],NULL,NULL);

			IS         = identity_d(NvnS[Pb]);  // free
			ChiRefS_vS = basis(Pb,rst_vS[0],NvnS[Pb],&Nbf,dE); // free
			free(rst_vS[0]);

			if (strstr(BasisType,"Modal")) {
				ChiS_vS[P][Pb][0] = ChiRefS_vS;
			} else if (strstr(BasisType,"Nodal")) {
				ChiS_vS[P][Pb][0] = IS;
			}

			if (ChiInvS_vS[Pb][Pb][0] == NULL)
				ChiInvS_vS[Pb][Pb][0] = inverse_d(NvnS[Pb],NvnS[Pb],ChiS_vS[P][Pb][0],IS); // keep

			if (Pb != P) {
				free(IS);
				free(ChiRefS_vS);
			} else {
				// Avoid negative Pb for P = 0
				break;
			}
		}
		set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,PGc[P],cubature); // free
		set_from_cubdata(CUBDATA,&NvnGc[P],NULL,&rst_vGc,NULL,NULL);

		set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,PCs[P][Eclass],cubature); // free
		set_from_cubdata(CUBDATA,&NvnCs[P],NULL,&rst_vCs,NULL,NULL);

		set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,PCc[P][Eclass],cubature); // free
		set_from_cubdata(CUBDATA,&NvnCc[P],NULL,&rst_vCc,NULL,NULL);

		// Preliminary Operators
		IGc[P][P][0] = identity_d(NvnGc[P]); // keep
		ICs[P][P][0] = identity_d(NvnCs[P]); // keep
		ICc[P][P][0] = identity_d(NvnCc[P]); // keep

		ChiRefGc_vGc = basis(PGc[P],rst_vGc,NvnGc[P],&Nbf,dE);         // free
		ChiRefCs_vCs = basis(PCs[P][Eclass],rst_vCs,NvnCs[P],&Nbf,dE); // free
		ChiRefCc_vCc = basis(PCc[P][Eclass],rst_vCc,NvnCc[P],&Nbf,dE); // free

		if (strstr(BasisType,"Modal")) {
			ChiGc_vGc = ChiRefGc_vGc;
			ChiCs_vCs = ChiRefCs_vCs;
			ChiCc_vCc = ChiRefCc_vCc;
		} else if (strstr(BasisType,"Nodal")) {
			ChiGc_vGc = IGc[P][P][0];
			ChiCs_vCs = ICs[P][P][0];
			ChiCc_vCc = ICc[P][P][0];
		}

		ChiRefInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiRefGc_vGc,IGc[P][P][0]); // free
		ChiRefInvCs_vCs = inverse_d(NvnCs[P],NvnCs[P],ChiRefCs_vCs,ICs[P][P][0]); // free
		ChiRefInvCc_vCc = inverse_d(NvnCc[P],NvnCc[P],ChiRefCc_vCc,ICc[P][P][0]); // free
		ChiRefInvS_vS   = inverse_d(NvnS[P], NvnS[P], ChiRefS_vS,IS);             // free

		ChiInvGc_vGc        = inverse_d(NvnGc[P],NvnGc[P],ChiGc_vGc,IGc[P][P][0]); // free
		ChiInvCs_vCs        = inverse_d(NvnCs[P],NvnCs[P],ChiCs_vCs,ICs[P][P][0]); // free
		ChiInvCc_vCc        = inverse_d(NvnCc[P],NvnCc[P],ChiCc_vCc,ICc[P][P][0]); // free

		TGc = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGc[P],NvnGc[P],1.0,ChiRefInvGc_vGc,ChiGc_vGc);      // free
		TCs = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,ChiRefInvCs_vCs,ChiCs_vCs);      // free
		TCc = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,ChiRefInvCc_vCc,ChiCc_vCc);      // free

		TS[P][P][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnS[P], NvnS[P], 1.0,ChiRefInvS_vS,ChiS_vS[P][P][0]); // keep

		if (EType != TET && EType != PYR) {
			set_cubdata(CUBDATA,false,false,NodeTypeS[P][Eclass],dE,P,cubature); // free
			set_from_cubdata(CUBDATA,&NvnS[P],NULL,&rst_vS[0],NULL,NULL);

			ChiBezS_vS             = basis_Bezier(P,rst_vS[0],NvnS[P],&Nbf,dE); // free
			ChiBezInvS_vS[P][P][0] = inverse_d(NvnS[P],NvnS[P],ChiBezS_vS,IS);  // keep
			TS_vB[P][P][0]         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P],NvnS[P],NvnS[P],1.0,ChiBezInvS_vS[P][P][0],ChiS_vS[P][P][0]); // keep

			TInvS_vB[P][P][0] = inverse_d(NvnS[P],NvnS[P],TS_vB[P][P][0],IS); // keep

			free(rst_vS[0]);
			free(ChiBezS_vS);
		}

		free(IS);

		free(ChiRefGc_vGc);
		free(ChiRefCs_vCs);
		free(ChiRefCc_vCc);
		free(ChiRefS_vS);

		free(ChiRefInvGc_vGc);
		free(ChiRefInvCs_vCs);
		free(ChiRefInvCc_vCc);
		free(ChiRefInvS_vS);

		for (size_t Pb = PbMin; Pb <= PbMax; Pb++) {
			if (w_vIs[Pb])
				free(w_vIs[Pb]);
			set_cubdata(CUBDATA,true,false,NodeTypeIvs[Pb][Eclass],dE,PIvs[Pb][Eclass],cubature); // free
			set_from_cubdata(CUBDATA,&NvnIs[Pb],NULL,&rst_vIs[0],&w_vIs[Pb],NULL);
			free(rst_vIs[0]);

			if (w_vIc[Pb])
				free(w_vIc[Pb]);
			set_cubdata(CUBDATA,true,false,NodeTypeIvc[Pb][Eclass],dE,PIvc[Pb][Eclass],cubature); // free
			set_from_cubdata(CUBDATA,&NvnIc[Pb],NULL,&rst_vIc[0],&w_vIc[Pb],NULL);
			free(rst_vIc[0]);

			wInv_vIs = malloc(NvnIs[P] * sizeof *wInv_vIs); // free
			wInv_vIc = malloc(NvnIc[P] * sizeof *wInv_vIc); // free

			for (i = 0, iMax = NvnIs[P]; i < iMax; i++)
				wInv_vIs[i] = 1./w_vIs[P][i];

			for (i = 0, iMax = NvnIc[P]; i < iMax; i++)
				wInv_vIc[i] = 1./w_vIc[P][i];

			diag_w_vIs    = diag_d(w_vIs[Pb],NvnIs[Pb]); // free
			diag_w_vIc    = diag_d(w_vIc[Pb],NvnIc[Pb]); // free
			diag_wInv_vIs = diag_d(wInv_vIs, NvnIs[P]);  // free
			diag_wInv_vIc = diag_d(wInv_vIc, NvnIc[P]);  // free

			free(wInv_vIs);
			free(wInv_vIc);

			// VOLUME Operators
			for (i = iMax = NEhref; i--; ) {
				ELEMENT_h = get_ELEMENT_type(EType_h[i]);

				Eclass = get_Eclass(EType_h[i]);
				select_functions(&basis,&grad_basis,&cubature,EType_h[i]);

				Nve_h = ELEMENT_h->Nve;
				NvnGs = ELEMENT_h->NvnGs;
				NvnS  = ELEMENT_h->NvnS;
				NvnIs = ELEMENT_h->NvnIs;
				NvnIc = ELEMENT_h->NvnIc;

				set_cubdata(CUBDATA,false,false,NodeTypeS[Pb][Eclass],dE,Pb,cubature); // free
				set_from_cubdata(CUBDATA,&NvnS[Pb],NULL,&rst_vS[0],NULL,NULL);

				set_cubdata(CUBDATA,false,false,NodeTypeIvs[Pb][Eclass],dE,PIvs[Pb][Eclass],cubature); // free
				set_from_cubdata(CUBDATA,&NvnIs[Pb],NULL,&rst_vIs[0],NULL,NULL);

				set_cubdata(CUBDATA,false,false,NodeTypeIvc[Pb][Eclass],dE,PIvc[Pb][Eclass],cubature); // free
				set_from_cubdata(CUBDATA,&NvnIc[Pb],NULL,&rst_vIc[0],NULL,NULL);

				ChiRefGs_vS  = basis(PGs,rst_vS[0], NvnS[Pb], &Nbf,dE); // free
				ChiRefGs_vIs = basis(PGs,rst_vIs[0],NvnIs[Pb],&Nbf,dE); // free
				ChiRefGs_vIc = basis(PGs,rst_vIc[0],NvnIc[Pb],&Nbf,dE); // free

				free(E_rst_vV);

				E_rst_vV        = get_rst_vV(ELEMENT_h);                   // free
				IGs             = identity_d(Nve_h);                       // free
				ChiRefGs_vGs    = basis(PGs,E_rst_vV,Nve_h,&Nbf,dE);       // free
				ChiRefInvGs_vGs = inverse_d(Nve_h,Nve_h,ChiRefGs_vGs,IGs); // free
				free(IGs);
				free(ChiRefGs_vGs);

				BCoords_V[i]->S[Pb]  = mm_Alloc_d(CBCM,CBT,CBT,NvnS[Pb], Nve_h,Nve_h,1.0,ChiRefGs_vS, ChiRefInvGs_vGs); // free
				BCoords_V[i]->Is[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIs[Pb],Nve_h,Nve_h,1.0,ChiRefGs_vIs,ChiRefInvGs_vGs); // free
				BCoords_V[i]->Ic[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIc[Pb],Nve_h,Nve_h,1.0,ChiRefGs_vIc,ChiRefInvGs_vGs); // free

				free(ChiRefInvGs_vGs);

				if (i) {
					free(rst_vS[0]);
					free(rst_vIs[0]);
					free(rst_vIc[0]);

					free(ChiRefGs_vIs);
					free(ChiRefGs_vIc);
				}

				free(ChiRefGs_vS);
			}

			// Note: Nvref >= NvrefSF
			for (vh = 0; vh < Nvref; vh++) {
				vrefSF = vh; // Used for operators defined only for setup_TP_operators
				mm_CTN_d(Nvve[vh],dE,Nve,VeV[vh],E_rst_vV,rst_vV);

				if (vh) {
					Indh = get_IndEhref(EType,vh);
					if (EType == TET || EType == PYR)
						ELEMENT_h = get_ELEMENT_type(EType_h[Indh]);
					else
						ELEMENT_h = ELEMENT;

					rst_vS[vh]  = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], dE,Nvve[vh],1.0,BCoords_V[Indh]->S[Pb], rst_vV); // free
					rst_vIs[vh] = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],dE,Nvve[vh],1.0,BCoords_V[Indh]->Is[Pb],rst_vV); // free
					rst_vIc[vh] = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],dE,Nvve[vh],1.0,BCoords_V[Indh]->Ic[Pb],rst_vV); // free

					ChiRefGs_vIs = basis(PGs,rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
					ChiRefGs_vIc = basis(PGs,rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				}

// Note: Nvve[vh] == ELEMENT_h->NvnGs[1]. Replace. (ToBeDeleted)
				ChiRefGs_vGs = basis(PGs,           rst_vV,     Nvve[vh], &Nbf,dE); // free
				ChiRefGs_vS  = basis(PGs,           rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefGc_vIs = basis(PGc[P],        rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefGc_vIc = basis(PGc[P],        rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefGc_vS  = basis(PGc[P],        rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefCs_vS  = basis(PCs[P][Eclass],rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefCs_vIs = basis(PCs[P][Eclass],rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefCs_vIc = basis(PCs[P][Eclass],rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefCc_vS  = basis(PCc[P][Eclass],rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefCc_vIs = basis(PCc[P][Eclass],rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefCc_vIc = basis(PCc[P][Eclass],rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefS_vS   = basis(P,             rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefS_vIs  = basis(P,             rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefS_vIc  = basis(P,             rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free

				ChiGs_vGs = mm_Alloc_d(CBRM,CBNT,CBNT,Nvve[vh], NvnGs[1],NvnGs[1],1.0,ChiRefGs_vGs,TGs[1][1][0]); // free
				ChiGs_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnGs[1],NvnGs[1],1.0,ChiRefGs_vS, TGs[1][1][0]); // free
				ChiGs_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vIs,TGs[1][1][0]); // free
				ChiGs_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vIc,TGs[1][1][0]); // free
				ChiGc_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnGc[P],NvnGc[P],1.0,ChiRefGc_vS, TGc); // free
				ChiGc_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIs,TGc); // free
				ChiGc_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIc,TGc); // free
				ChiCs_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCs[P],NvnCs[P],1.0,ChiRefCs_vS, TCs); // free
				ChiCs_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIs,TCs); // free
				ChiCs_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIc,TCs); // free
				ChiCc_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCc[P],NvnCc[P],1.0,ChiRefCc_vS, TCc); // free
				ChiCc_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIs,TCc); // free
				ChiCc_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIc,TCc); // free

				// Returned SF Operators
				ChiS_vIc[P][Pb][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIc,TS[P][P][0]); // keep
				if (vrefSF < NvrefSF) {
					ChiS_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnS[P],NvnS[P],1.0,ChiRefS_vS,TS[P][P][0]);  // keep
					ChiS_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIs,TS[P][P][0]); // keep

					if (P == Pb) {
						I_vGs_vS[1][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb] ,NvnGs[1],NvnGs[1],1.0,ChiGs_vS, ChiInvGs_vGs[1][1][0]); // keep
						I_vGs_vIs[1][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGs[1],NvnGs[1],1.0,ChiGs_vIs,ChiInvGs_vGs[1][1][0]); // keep
						I_vGs_vIc[1][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGs[1],NvnGs[1],1.0,ChiGs_vIc,ChiInvGs_vGs[1][1][0]); // keep
					}
					I_vGc_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnGc[P],NvnGc[P],1.0,ChiGc_vS, ChiInvGc_vGc); // keep
					I_vGc_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGc[P],NvnGc[P],1.0,ChiGc_vIs,ChiInvGc_vGc); // keep
					I_vGc_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGc[P],NvnGc[P],1.0,ChiGc_vIc,ChiInvGc_vGc); // keep
					I_vCs_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCs[P],NvnCs[P],1.0,ChiCs_vS, ChiInvCs_vCs); // keep
					I_vCs_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCs[P],NvnCs[P],1.0,ChiCs_vIs,ChiInvCs_vCs); // keep
					I_vCs_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCs[P],NvnCs[P],1.0,ChiCs_vIc,ChiInvCs_vCs); // keep
					I_vCc_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCc[P],NvnCc[P],1.0,ChiCc_vS, ChiInvCc_vCc); // keep
					I_vCc_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCc[P],NvnCc[P],1.0,ChiCc_vIs,ChiInvCc_vCc); // keep
					I_vCc_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCc[P],NvnCc[P],1.0,ChiCc_vIc,ChiInvCc_vCc); // keep

					if (EFE) {
						// Incorrect operator for PYR -> TET (vrefSF != 0) but this operator is not used.
						Is_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[Pb],NvnIs[Pb], 1.0,ChiS_vIs[P][Pb][vrefSF],diag_w_vIs); // keep
						Ic_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[Pb],NvnIc[Pb], 1.0,ChiS_vIc[P][Pb][vrefSF],diag_w_vIc); // keep

						Is_Strong_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[Pb],NvnS[P],NvnIs[Pb],1.0,diag_w_vIs,ChiS_vIs[P][Pb][vrefSF]); // keep
						Ic_Strong_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[Pb],NvnS[P],NvnIc[Pb],1.0,diag_w_vIc,ChiS_vIc[P][Pb][vrefSF]); // keep
					} else {
						printf("Error: Unsupported EFE.\n"), EXIT_MSG;
					}

					if (Collocated) {
						// Incorrect operator for PYR -> TET (vrefSF != 0) but this operator is not used.
						dummyPtr_d = Is_Weak_VV[P][Pb][vrefSF];
						Is_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[Pb],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Ic_Weak_VV[P][Pb][vrefSF];
						Ic_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnIc[Pb],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);

						dummyPtr_d = Is_Strong_VV[P][Pb][vrefSF];
						Is_Strong_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[Pb],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Ic_Strong_VV[P][Pb][vrefSF];
						Ic_Strong_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnIc[Pb],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);
					}
				}

				// Returned Adaptation Operators
				if (vh == 0 || P == Pb) {
					i = get_IndEhref(EType,vh);
					dummyPtr_d = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnS[P],NvnS[P],1.0,ChiRefS_vS,TS[P][P][0]);  // free
					Ihat_vS_vS[P][Pb][vh]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb],NvnS[P],ELEMENT_h->NvnS[Pb],
					                                    1.0,ELEMENT_h->ChiInvS_vS[Pb][Pb][0],dummyPtr_d); // keep
					free(dummyPtr_d);
					if (P == PGlobal && Pb == PGlobal) {
						I_vGs_vGs[1][1][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,Nvve[vh],NvnGs[1],NvnGs[1],1.0,ChiGs_vGs,ChiInvGs_vGs[1][1][0]); // keep
					}
				}

				free(ChiRefGs_vGs);
				free(ChiRefGs_vS);
				free(ChiRefGs_vIs);
				free(ChiRefGs_vIc);
				free(ChiRefGc_vS);
				free(ChiRefGc_vIs);
				free(ChiRefGc_vIc);
				free(ChiRefCs_vS);
				free(ChiRefCs_vIs);
				free(ChiRefCs_vIc);
				free(ChiRefCc_vS);
				free(ChiRefCc_vIs);
				free(ChiRefCc_vIc);
				free(ChiRefS_vS);
				free(ChiRefS_vIs);
				free(ChiRefS_vIc);

				free(ChiGs_vGs);
				free(ChiGs_vS);
				free(ChiGs_vIs);
				free(ChiGs_vIc);
				free(ChiGc_vS);
				free(ChiGc_vIs);
				free(ChiGc_vIc);
				free(ChiCs_vS);
				free(ChiCs_vIs);
				free(ChiCs_vIc);
				free(ChiCc_vS);
				free(ChiCc_vIs);
				free(ChiCc_vIc);
			}

			plotting_element_info(&rst_vP,&dummyPtr_ui[0],&dummyPtr_ui[1],&dummyPtr_ui[2],&NvnP,&dummy_ui,max(Pb,u1),EType); // free
			free(dummyPtr_ui[0]); free(dummyPtr_ui[1]); free(dummyPtr_ui[2]);

			ChiRefGc_vP  = basis(PGc[P],        rst_vP,    NvnP,    &Nbf,dE); // free
			ChiRefS_vP   = basis(P,             rst_vP,    NvnP,    &Nbf,dE); // free

			ChiGc_vP            = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP, TGc); // free

			// Returned Operators
			ChiS_vP[P][Pb][0]   = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnS[P], NvnS[P], 1.0,ChiRefS_vP,TS[P][P][0]);           // keep
			I_vGc_vP[P][Pb][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiGc_vP,  ChiInvGc_vGc); // keep

			free(ChiRefGc_vP);
			free(ChiRefS_vP);
			free(ChiGc_vP);

			// Operators (1/P -> 2/P only)
			if (P == Pb) {
				ChiRefGs_vP  = basis(PGs,           rst_vP,    NvnP,    &Nbf,dE); // free
				ChiRefGs_vGc = basis(PGs,           rst_vGc,   NvnGc[P],&Nbf,dE); // free
				ChiRefGs_vCs = basis(PGs,           rst_vCs,   NvnCs[P],&Nbf,dE); // free
				ChiRefGc_vCc = basis(PGc[P],        rst_vCc,   NvnCc[P],&Nbf,dE); // free
				ChiRefS_vIs  = basis(P,             rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				ChiRefS_vIc  = basis(P,             rst_vIc[0],NvnIc[P],&Nbf,dE); // free

				ChiGs_vP          = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGs[1],NvnGs[1],1.0,ChiRefGs_vP, TGs[1][1][0]); // free
				ChiGs_vGc         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vGc,TGs[1][1][0]); // free
				ChiGs_vCs         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vCs,TGs[1][1][0]); // free
				ChiGc_vCc         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vCc,TGc); // free

				GradChiRefGs_vCs = grad_basis(PGs,           rst_vCs,   NvnCs[P],&Nbf,dE); // free
				GradChiRefGs_vIs = grad_basis(PGs,           rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				GradChiRefGc_vCc = grad_basis(PGc[P],        rst_vCc,   NvnCc[P],&Nbf,dE); // free
				GradChiRefGc_vIc = grad_basis(PGc[P],        rst_vIc[0],NvnIc[P],&Nbf,dE); // free
				GradChiRefCs_vCs = grad_basis(PCs[P][Eclass],rst_vCs,   NvnCs[P],&Nbf,dE); // free
				GradChiRefCc_vCc = grad_basis(PCc[P][Eclass],rst_vCc,   NvnCc[P],&Nbf,dE); // free
				GradChiRefS_vS   = grad_basis(P,             rst_vS[0] ,NvnS[P], &Nbf,dE); // free
				GradChiRefS_vIs  = grad_basis(P,             rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				GradChiRefS_vIc  = grad_basis(P,             rst_vIc[0],NvnIc[P],&Nbf,dE); // free

				for (dim = 0; dim < dE; dim++) {
					GradChiGs_vCs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_vCs[dim],TGs[1][1][0]); // free
					GradChiGs_vIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_vIs[dim],TGs[1][1][0]); // free
					GradChiGc_vCc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vCc[dim],TGc); // free
					GradChiGc_vIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vIc[dim],TGc); // free
					GradChiCs_vCs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiRefCs_vCs[dim],TCs); // free
					GradChiCc_vCc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiRefCc_vCc[dim],TCc); // free

					GradChiS_vS[P][Pb][0][dim]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnS[P],NvnS[P],1.0,GradChiRefS_vS[dim],TS[P][P][0]);  // keep
					GradChiS_vIs[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnS[P],NvnS[P],1.0,GradChiRefS_vIs[dim],TS[P][P][0]); // keep
					GradChiS_vIc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnS[P],NvnS[P],1.0,GradChiRefS_vIc[dim],TS[P][P][0]); // keep
				}

				// Returned Operators
				// VOLUME related operators
				I_vGs_vP[1][Pb][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGs[1],NvnGs[1],1.0,ChiGs_vP, ChiInvGs_vGs[1][1][0]); // keep
				I_vGs_vGc[1][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGs[1],NvnGs[1],1.0,ChiGs_vGc,ChiInvGs_vGs[1][1][0]); // keep
				I_vGs_vCs[1][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,ChiGs_vCs,ChiInvGs_vGs[1][1][0]); // keep
				I_vGc_vCc[P][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vCc,ChiInvGc_vGc); // keep

				for (dim = 0; dim < dE; dim++) {
					D_vGs_vCs[1][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,GradChiGs_vCs[dim],ChiInvGs_vGs[1][1][0]); // keep
					D_vGs_vIs[1][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[1],NvnGs[1],1.0,GradChiGs_vIs[dim],ChiInvGs_vGs[1][1][0]); // keep
					D_vGc_vCc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vCc[dim],ChiInvGc_vGc); // keep
					D_vGc_vIc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vIc[dim],ChiInvGc_vGc); // keep
					D_vCs_vCs[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiCs_vCs[dim],ChiInvCs_vCs); // keep
					D_vCc_vCc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiCc_vCc[dim],ChiInvCc_vCc); // keep

					if (EFE) {
						Ds_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[P],NvnIs[P],1.0,GradChiS_vIs[P][Pb][0][dim],diag_w_vIs); // keep
						Dc_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[P],NvnIc[P],1.0,GradChiS_vIc[P][Pb][0][dim],diag_w_vIc); // keep
						Ds_Strong_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnS[P],NvnIs[P],1.0,diag_w_vIs,GradChiS_vIs[P][Pb][0][dim]); // keep
						Dc_Strong_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnS[P],NvnIc[P],1.0,diag_w_vIc,GradChiS_vIc[P][Pb][0][dim]); // keep
					} else {
						printf("Error: Unsupported EFE.\n"), EXIT_MSG;
					}
					if (Collocated) {
						dummyPtr_d = Ds_Weak_VV[P][Pb][0][dim];
						Ds_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[Pb],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Dc_Weak_VV[P][Pb][0][dim];
						Dc_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnIc[Pb],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);

						dummyPtr_d = Ds_Strong_VV[P][Pb][0][dim];
						Ds_Strong_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[P],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Dc_Strong_VV[P][Pb][0][dim];
						Dc_Strong_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[P],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
					}
				}

				free(ChiRefGs_vP);
				free(ChiRefGs_vGc);
				free(ChiRefGs_vCs);
				free(ChiRefGc_vCc);
				free(ChiRefS_vIs);
				free(ChiRefS_vIc);

				free(ChiGs_vP);
				free(ChiGs_vGc);
				free(ChiGs_vCs);
				free(ChiGc_vCc);

				array_free2_d(dE,GradChiRefGs_vCs);
				array_free2_d(dE,GradChiRefGs_vIs);
				array_free2_d(dE,GradChiRefGc_vCc);
				array_free2_d(dE,GradChiRefGc_vIc);
				array_free2_d(dE,GradChiRefCs_vCs);
				array_free2_d(dE,GradChiRefCc_vCc);
				array_free2_d(dE,GradChiRefS_vS);
				array_free2_d(dE,GradChiRefS_vIs);
				array_free2_d(dE,GradChiRefS_vIc);

				for (dim = 0; dim < dE; dim++) {
					free(GradChiGs_vCs[dim]);
					free(GradChiGs_vIs[dim]);
					free(GradChiGc_vCc[dim]);
					free(GradChiGc_vIc[dim]);
					free(GradChiCs_vCs[dim]);
					free(GradChiCc_vCc[dim]);
				}
			}

			// FACE related operators
			for (f = 0; f < Nf; f++) {
				IndFType = get_IndFType(Eclass,f);

				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;

					mm_CTN_d(Nfve[f],dE,Nve,VeF[Vf],E_rst_vV,rst_vV);

					rst_fGc = mm_Alloc_d(CBCM,CBNT,CBNT,NfnGc[Pb][IndFType],dE,BF_Nve[IndFType],1.0,BCoords_F[IndFType]->Gc[Pb],rst_vV); // free
					rst_fS  = mm_Alloc_d(CBCM,CBNT,CBNT,NfnS[Pb][IndFType], dE,BF_Nve[IndFType],1.0,BCoords_F[IndFType]->S[Pb], rst_vV); // free
					rst_fIs = mm_Alloc_d(CBCM,CBNT,CBNT,NfnIs[Pb][IndFType],dE,BF_Nve[IndFType],1.0,BCoords_F[IndFType]->Is[Pb],rst_vV); // free
					rst_fIc = mm_Alloc_d(CBCM,CBNT,CBNT,NfnIc[Pb][IndFType],dE,BF_Nve[IndFType],1.0,BCoords_F[IndFType]->Ic[Pb],rst_vV); // free

					diag_w_fIs = diag_d(w_fIs[Pb][IndFType],NfnIs[Pb][IndFType]); // free
					diag_w_fIc = diag_d(w_fIc[Pb][IndFType],NfnIc[Pb][IndFType]); // free

					ChiRefS_fS  = basis(P,rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefS_fIs = basis(P,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefS_fIc = basis(P,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					GradChiRefS_fIs = grad_basis(P,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					GradChiRefS_fIc = grad_basis(P,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					// Returned Operators
					ChiS_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnS[P],NvnS[P],1.0,ChiRefS_fS, TS[P][P][0]); // keep
					ChiS_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnS[P],NvnS[P],1.0,ChiRefS_fIs,TS[P][P][0]); // keep
					ChiS_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnS[P],NvnS[P],1.0,ChiRefS_fIc,TS[P][P][0]); // keep

					for (dim = 0; dim < dE; dim++) {
						GradChiS_fIs[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnS[P],NvnS[P],1.0,GradChiRefS_fIs[dim],TS[P][P][0]); // keep
						GradChiS_fIc[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnS[P],NvnS[P],1.0,GradChiRefS_fIc[dim],TS[P][P][0]); // keep
					}

					Is_Weak_FV[P][Pb][Vf] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NfnIs[Pb][IndFType],NfnIs[Pb][IndFType],-1.0,ChiS_fIs[P][Pb][Vf],diag_w_fIs); // keep
					Ic_Weak_FV[P][Pb][Vf] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NfnIc[Pb][IndFType],NfnIc[Pb][IndFType],-1.0,ChiS_fIc[P][Pb][Vf],diag_w_fIc); // keep

					if (Collocated) {
						dummyPtr_d = Is_Weak_FV[P][Pb][Vf];
						Is_Weak_FV[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NfnIs[Pb][IndFType],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Ic_Weak_FV[P][Pb][Vf];
						Ic_Weak_FV[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NfnIc[Pb][IndFType],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);
					}

					if (Collocated || VFPartUnity[Eclass]) {
						convert_to_CSR_d(NfnIs[Pb][IndFType],NvnS[P],ChiS_fIs[P][Pb][Vf],&ChiS_fIs_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NfnIc[Pb][IndFType],NvnS[P],ChiS_fIc[P][Pb][Vf],&ChiS_fIc_sp[P][Pb][Vf]); // keep

						convert_to_CSR_d(NvnS[P],NfnIs[Pb][IndFType],Is_Weak_FV[P][Pb][Vf],&Is_Weak_FV_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NvnS[P],NfnIc[Pb][IndFType],Ic_Weak_FV[P][Pb][Vf],&Ic_Weak_FV_sp[P][Pb][Vf]); // keep
					}

					ChiRefGs_fS  = basis(PGs           ,rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefGs_fIs = basis(PGs           ,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefGs_fIc = basis(PGs           ,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					ChiRefGc_fGc = basis(PGc[P]        ,rst_fGc,NfnGc[Pb][IndFType],&Nbf,dE); // free
					ChiRefGc_fS  = basis(PGc[P]        ,rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefGc_fIs = basis(PGc[P]        ,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefGc_fIc = basis(PGc[P]        ,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					ChiRefCs_fS  = basis(PCs[P][Eclass],rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefCs_fIs = basis(PCs[P][Eclass],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefCs_fIc = basis(PCs[P][Eclass],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					ChiRefCc_fS  = basis(PCc[P][Eclass],rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefCc_fIs = basis(PCc[P][Eclass],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefCc_fIc = basis(PCc[P][Eclass],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					GradChiRefGs_fIs = grad_basis(PGs,   rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					GradChiRefGs_fIc = grad_basis(PGs,   rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					GradChiRefGc_fIs = grad_basis(PGc[P],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					GradChiRefGc_fIc = grad_basis(PGc[P],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					ChiGs_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGs[1],NvnGs[1],1.0,ChiRefGs_fS, TGs[1][1][0]); // free
					ChiGs_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiRefGs_fIs,TGs[1][1][0]); // free
					ChiGs_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiRefGs_fIc,TGs[1][1][0]); // free
					ChiGc_fGc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnGc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fGc,TGc); // free
					ChiGc_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGc[P],NvnGc[P],1.0,ChiRefGc_fS, TGc); // free
					ChiGc_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIs,TGc); // free
					ChiGc_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIc,TGc); // free
					ChiCs_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCs[P],NvnCs[P],1.0,ChiRefCs_fS, TCs); // free
					ChiCs_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIs,TCs); // free
					ChiCs_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIc,TCs); // free
					ChiCc_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCc[P],NvnCc[P],1.0,ChiRefCc_fS, TCc); // free
					ChiCc_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIs,TCc); // free
					ChiCc_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIc,TCc); // free

					for (dim = 0; dim < dE; dim++) {
						GradChiGs_fIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_fIs[dim],TGs[1][1][0]); // free
						GradChiGs_fIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_fIc[dim],TGs[1][1][0]); // free
						GradChiGc_fIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_fIs[dim],TGc); // free
						GradChiGc_fIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_fIc[dim],TGc); // free
					}

					// Returned Operators
					if (P == Pb) {
						// Outside of fh == 0 if condition as operator is used for testing h-adaptation connectivity.
						I_vGs_fS[1][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGs[1],NvnGs[1],1.0,ChiGs_fS, ChiInvGs_vGs[1][1][0]); // keep
					}

					if (fh == 0) {
						if (P == Pb) {
							I_vGs_fIs[1][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiGs_fIs,ChiInvGs_vGs[1][1][0]); // keep
							I_vGs_fIc[1][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiGs_fIc,ChiInvGs_vGs[1][1][0]); // keep

							// Vertex blending operators
							if (P == 2) {
								rst_fG2      = mm_Alloc_d(CBCM,CBNT,CBNT,NfnG2[Pb][IndFType],dE,BF_Nve[IndFType],1.0,BCoords_F[IndFType]->G2[Pb],rst_vV); // free
								ChiRefG2_fG2 = basis(P,rst_fG2,NfnG2[Pb][IndFType],&Nbf,dE); // free
								ChiG2_fG2    = mm_Alloc_d(CBRM,CBNT,CBNT,NfnG2[Pb][IndFType],NvnG2[P],NvnG2[P],1.0,ChiRefG2_fG2,TG2); // free
								I_vG2_fG2[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnG2[Pb][IndFType],NvnG2[P],NvnG2[P],1.0,ChiG2_fG2,ChiInvG2_vG2); // keep

								free(rst_fG2);
								free(ChiRefG2_fG2);
								free(ChiG2_fG2);
							}
						}
						I_vGc_fGc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnGc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fGc,ChiInvGc_vGc); // keep
						I_vGc_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGc[P],NvnGc[P],1.0,ChiGc_fS, ChiInvGc_vGc); // keep
						I_vGc_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIs,ChiInvGc_vGc); // keep
						I_vGc_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIc,ChiInvGc_vGc); // keep
						I_vCs_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCs[P],NvnCs[P],1.0,ChiCs_fS, ChiInvCs_vCs); // keep
						I_vCs_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIs,ChiInvCs_vCs); // keep
						I_vCs_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIc,ChiInvCs_vCs); // keep
						I_vCc_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCc[P],NvnCc[P],1.0,ChiCc_fS, ChiInvCc_vCc); // keep
						I_vCc_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIs,ChiInvCc_vCc); // keep
						I_vCc_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIc,ChiInvCc_vCc); // keep
					}

					for (dim = 0; dim < dE; dim++) {
						if (P == Pb) {
							D_vGs_fIs[1][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiGs_fIs[dim],ChiInvGs_vGs[1][1][0]); // keep
							D_vGs_fIc[1][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiGs_fIc[dim],ChiInvGs_vGs[1][1][0]); // keep
						}
						D_vGc_fIs[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiGc_fIs[dim],ChiInvGc_vGc); // keep
						D_vGc_fIc[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiGc_fIc[dim],ChiInvGc_vGc); // keep
					}

					free(ChiRefGs_fS);
					free(ChiRefGs_fIs);
					free(ChiRefGs_fIc);
					free(ChiRefGc_fGc);
					free(ChiRefGc_fS);
					free(ChiRefGc_fIs);
					free(ChiRefGc_fIc);
					free(ChiRefCs_fS);
					free(ChiRefCs_fIs);
					free(ChiRefCs_fIc);
					free(ChiRefCc_fS);
					free(ChiRefCc_fIs);
					free(ChiRefCc_fIc);

					free(ChiGs_fS);
					free(ChiGs_fIs);
					free(ChiGs_fIc);
					free(ChiGc_fGc);
					free(ChiGc_fS);
					free(ChiGc_fIs);
					free(ChiGc_fIc);
					free(ChiCs_fS);
					free(ChiCs_fIs);
					free(ChiCs_fIc);
					free(ChiCc_fS);
					free(ChiCc_fIs);
					free(ChiCc_fIc);

					array_free2_d(dE,GradChiRefGs_fIs);
					array_free2_d(dE,GradChiRefGs_fIc);
					array_free2_d(dE,GradChiRefGc_fIs);
					array_free2_d(dE,GradChiRefGc_fIc);

					for (dim = 0; dim < dE; dim++) {
						free(GradChiGs_fIs[dim]);
						free(GradChiGs_fIc[dim]);
						free(GradChiGc_fIs[dim]);
						free(GradChiGc_fIc[dim]);
					}

					free(rst_fGc);
					free(rst_fS);
					free(rst_fIs);
					free(rst_fIc);

					free(diag_w_fIs);
					free(diag_w_fIc);

					free(ChiRefS_fS);
					free(ChiRefS_fIs);
					free(ChiRefS_fIc);

					array_free2_d(dE,GradChiRefS_fIs);
					array_free2_d(dE,GradChiRefS_fIc);
				}
			}

			// EDGE related operators (d == 3 only)
			for (e = 0; dE == DMAX && e < Ne; e++) {
				Ve = e*NEREFMAX;

				mm_CTN_d(Neve,dE,Nve,VeE[Ve],E_rst_vV,rst_vV);

				rst_eGc = mm_Alloc_d(CBCM,CBNT,CBNT,NenGc[Pb],dE,BE_Nve,1.0,BCoords_E->Gc[Pb],rst_vV); // free

				ChiRefGc_eGc = basis(PGc[P]        ,rst_eGc,NenGc[Pb],&Nbf,dE);                             // free
				ChiGc_eGc    = mm_Alloc_d(CBRM,CBNT,CBNT,NenGc[Pb],NvnGc[P],NvnGc[P],1.0,ChiRefGc_eGc,TGc); // free

				I_vGc_eGc[P][Pb][Ve] = mm_Alloc_d(CBRM,CBNT,CBNT,NenGc[Pb],NvnGc[P],NvnGc[P],1.0,ChiGc_eGc,ChiInvGc_vGc); // keep

				if (P == Pb && P == 2) {
					rst_eG2 = mm_Alloc_d(CBCM,CBNT,CBNT,NenG2[Pb],dE,BE_Nve,1.0,BCoords_E->G2[Pb],rst_vV); // free

					ChiRefG2_eG2 = basis(P             ,rst_eG2,NenG2[Pb],&Nbf,dE);                             // free
					ChiG2_eG2    = mm_Alloc_d(CBRM,CBNT,CBNT,NenG2[Pb],NvnG2[P],NvnG2[P],1.0,ChiRefG2_eG2,TG2); // free

					I_vG2_eG2[P][Pb][Ve] = mm_Alloc_d(CBRM,CBNT,CBNT,NenG2[Pb],NvnG2[P],NvnG2[P],1.0,ChiG2_eG2,ChiInvG2_vG2); // keep

					free(rst_eG2);
					free(ChiRefG2_eG2);
					free(ChiG2_eG2);
				}

				free(rst_eGc);
				free(ChiRefGc_eGc);
				free(ChiGc_eGc);
			}

			for (vref = 0; vref < Nvref; vref++) {
				free(rst_vS[vref]);
				free(rst_vIs[vref]);
				free(rst_vIc[vref]);
			}
			for (i = 0; i < NEhref; i++) {
				free(BCoords_V[i]->S[Pb]);
				free(BCoords_V[i]->Is[Pb]);
				free(BCoords_V[i]->Ic[Pb]);
			}

			free(rst_vP);

			free(diag_w_vIs);
			free(diag_w_vIc);
			free(diag_wInv_vIs);
			free(diag_wInv_vIc);
		}
		free(rst_vGc);
		free(rst_vCs);
		free(rst_vCc);

		free(ChiInvGc_vGc);
		free(ChiInvCs_vCs);
		free(ChiInvCc_vCc);

		free(TGc);
		free(TCs);
		free(TCc);
	}
	free(rst_vS);
	free(rst_vIs);
	free(rst_vIc);

	free(ChiInvG2_vG2);

	free(TG2);

	if (dE == DMAX) {
		array_free2_d(PMax+1,BCoords_E->G2);
		array_free2_d(PMax+1,BCoords_E->Gc);
		free(BCoords_E);
		free(BCoords_dEm2);
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		array_free2_d(PMax+1,BCoords_F[IndFType]->G2);
		array_free2_d(PMax+1,BCoords_F[IndFType]->Gc);
		array_free2_d(PMax+1,BCoords_F[IndFType]->S);
		array_free2_d(PMax+1,BCoords_F[IndFType]->Is);
		array_free2_d(PMax+1,BCoords_F[IndFType]->Ic);
		free(BCoords_F[IndFType]);
		free(BCoords_dEm1[IndFType]);
	}

	for (i = 0; i < NEhref; i++) {
		free(BCoords_V[i]->S);
		free(BCoords_V[i]->Is);
		free(BCoords_V[i]->Ic);
		free(BCoords_V[i]);
	}
	free(BCoords_V);

	free(E_rst_vV);
	free(rst_vV);

	free(GradChiGs_vCs);
	free(GradChiGs_vIs);
	free(GradChiGc_vCc);
	free(GradChiGc_vIc);
	free(GradChiCs_vCs);
	free(GradChiCc_vCc);

	free(GradChiGs_fIs);
	free(GradChiGs_fIc);
	free(GradChiGc_fIs);
	free(GradChiGc_fIc);

	free(ones_Nf);

	free(CUBDATA);

	compute_ELEMENT_Volume(EType);
}

static void setup_TP_operators(const unsigned int EType)
{
	/*
	 *	Purpose:
	 *		Compute operators for elements which are tensor-products of lower dimensional elements.
	 *
	 *	Comments:
	 *		This functions works by taking tensor-products of the operators of lower-dimensional ELEMENT operators. This
	 *		allows for a reduction in the computational complexity of computing these operators and is very intuitive
	 *		once the concept is understood.
	 *
	 *		The casts to "const *const ..." are required because the pointers have more than one level of dereferencing.
	 *		It would be beneficial to cast the stored matrices of setup_ELEMENT_operators to their appropriate fully
	 *		const types in the future, making the casts here redundant. This was not done immediately as it would likely
	 *		result in errors throughout the entire code, requiring changes to many functions.
	 *
	 *	References:
	 *		Thesis (ToBeModified)
	 */

	// Returned operators
	unsigned int *NvnGs, *NvnG2, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS, **NfnGc, **NfnS, **NfnIs,
	             **NfnIc;
	double       **w_vIs, **w_vIc, ***w_fIs, ***w_fIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS, ****ChiInvGs_vGs, ****ChiBezInvS_vS,
				 ****TGs, ****TS, ****TS_vB, ****TInvS_vB,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vG2, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vS, ****I_vGs_vIs,
	             ****I_vGc_vP,                               ****I_vGc_vCc,                ****I_vGc_vS, ****I_vGc_vIc,
	             ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             *****GradChiS_vIs, *****GradChiS_vIc,
	             ****L2hat_vS_vS,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vG2_fG2, ****I_vG2_eG2,
	             ****I_vGc_fGc, ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc, ****I_vGc_eGc,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
				 *****D_vGs_fIs, *****D_vGs_fIc,
				 *****D_vGc_fIs, *****D_vGc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FV, ****Ic_Weak_FV,
	             *****Ds_Weak_VV, *****Dc_Weak_VV, *****Ds_Strong_VV, *****Dc_Strong_VV;

	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               *****Ds_Weak_VV_sp, *****Dc_Weak_VV_sp,
	               ****Is_Weak_FV_sp, ****Ic_Weak_FV_sp;

	// Initialize DB Parameters
	unsigned int EFE          = DB.EFE,
	             PGlobal      = DB.PGlobal,
	             Collocated   = DB.Collocated,
	             *VFPartUnity = DB.VFPartUnity;

	// Standard datatypes
	unsigned int dim, vh, f, fh, e, fhMax, Nvref, IndClass,
	             Eclass, dE, Nf, Ne, Vf, Ve, *Nfref, *ones_Nf,
	             NIn[3], NOut[3];
	double const      * OP[3];

	struct S_ELEMENT *ELEMENT, *ELEMENTclass[2];

	// silence
	ELEMENTclass[1] = NULL;

	Eclass = get_Eclass(EType);
	ELEMENT = get_ELEMENT_type(EType);

	ELEMENTclass[0] = ELEMENT->ELEMENTclass[0];
	if (EType == WEDGE)
		ELEMENTclass[1] = ELEMENT->ELEMENTclass[1];

	dE    = ELEMENT->d;
	Nf    = ELEMENT->Nf;
	Ne    = ELEMENT->Ne;
	Nfref = ELEMENT->Nfref;
	Nvref = ELEMENT->Nvref;

	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnG2 = ELEMENT->NvnG2;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NvnS  = ELEMENT->NvnS;
	NfnGc = ELEMENT->NfnGc;
	NfnS  = ELEMENT->NfnS;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;

	w_fIs = ELEMENT->w_fIs;
	w_fIc = ELEMENT->w_fIc;

	w_vIs = ELEMENT->w_vIs;
	w_vIc = ELEMENT->w_vIc;

	ChiS_vP      = ELEMENT->ChiS_vP;
	ChiS_vS      = ELEMENT->ChiS_vS;
	ChiS_vIs     = ELEMENT->ChiS_vIs;
	ChiS_vIc     = ELEMENT->ChiS_vIc;
	ChiInvS_vS   = ELEMENT->ChiInvS_vS;
	ChiInvGs_vGs = ELEMENT->ChiInvGs_vGs;

	ChiBezInvS_vS = ELEMENT->ChiBezInvS_vS;

	GradChiS_vIs = ELEMENT->GradChiS_vIs;
	GradChiS_vIc = ELEMENT->GradChiS_vIc;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGs = ELEMENT->I_vGs_vGs;
	I_vGs_vG2 = ELEMENT->I_vGs_vG2;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vS  = ELEMENT->I_vGs_vS;
	I_vGs_vIs = ELEMENT->I_vGs_vIs;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vS  = ELEMENT->I_vGc_vS;
	I_vGc_vIc = ELEMENT->I_vGc_vIc;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	TGs   = ELEMENT->TGs;
	TS    = ELEMENT->TS;
	TS_vB = ELEMENT->TS_vB;

	TInvS_vB = ELEMENT->TInvS_vB;

	Ihat_vS_vS  = ELEMENT->Ihat_vS_vS;
	L2hat_vS_vS = ELEMENT->L2hat_vS_vS;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	ChiS_fS     = ELEMENT->ChiS_fS;
	ChiS_fIs    = ELEMENT->ChiS_fIs;
	ChiS_fIc    = ELEMENT->ChiS_fIc;
	ChiS_fIs_sp = ELEMENT->ChiS_fIs_sp;
	ChiS_fIc_sp = ELEMENT->ChiS_fIc_sp;

	I_vGs_fS  = ELEMENT->I_vGs_fS;
	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vG2_fG2 = ELEMENT->I_vG2_fG2;
	I_vGc_fGc = ELEMENT->I_vGc_fGc;
	I_vGc_fS  = ELEMENT->I_vGc_fS;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vG2_eG2 = ELEMENT->I_vG2_eG2;
	I_vGc_eGc = ELEMENT->I_vGc_eGc;
	I_vCs_fS  = ELEMENT->I_vCs_fS;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fS  = ELEMENT->I_vCc_fS;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	D_vGs_fIs = ELEMENT->D_vGs_fIs;
	D_vGs_fIc = ELEMENT->D_vGs_fIc;
	D_vGc_fIs = ELEMENT->D_vGc_fIs;
	D_vGc_fIc = ELEMENT->D_vGc_fIc;

	Is_Weak_VV    = ELEMENT->Is_Weak_VV;
	Ic_Weak_VV    = ELEMENT->Ic_Weak_VV;
	Is_Weak_FV    = ELEMENT->Is_Weak_FV;
	Ic_Weak_FV    = ELEMENT->Ic_Weak_FV;
	Is_Weak_FV_sp = ELEMENT->Is_Weak_FV_sp;
	Ic_Weak_FV_sp = ELEMENT->Ic_Weak_FV_sp;

	Ds_Weak_VV    = ELEMENT->Ds_Weak_VV;
	Dc_Weak_VV    = ELEMENT->Dc_Weak_VV;
	Ds_Strong_VV  = ELEMENT->Ds_Strong_VV;
	Dc_Strong_VV  = ELEMENT->Dc_Strong_VV;
	Ds_Weak_VV_sp = ELEMENT->Ds_Weak_VV_sp;
	Dc_Weak_VV_sp = ELEMENT->Dc_Weak_VV_sp;

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	if (Eclass == C_TP) {
		NvnGs[1] = pow(ELEMENTclass[0]->NvnGs[1],dE);
		NvnGs[2] = pow(ELEMENTclass[0]->NvnGs[2],dE);
		NvnG2[2] = pow(ELEMENTclass[0]->NvnG2[2],dE);

		get_sf_parametersV(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGs[1],(double const *const *const) ELEMENTclass[0]->TGs[1][1],
		                   0,0,NULL,NIn,NOut,OP,dE,0,Eclass);
		TGs[1][1][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
		get_sf_parametersV(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGs[1],(double const *const *const) ELEMENTclass[0]->ChiInvGs_vGs[1][1],
		                   0,0,NULL,NIn,NOut,OP,dE,0,Eclass);
		ChiInvGs_vGs[1][1][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
		get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnG2[2],ELEMENTclass[0]->I_vGs_vG2[1][2][0],
						  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
		I_vGs_vG2[1][2][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

		for (size_t P = PSMin; P <= PSMax; P++) {
			NvnGc[P] = pow(ELEMENTclass[0]->NvnGc[P],dE);
			NvnCs[P] = pow(ELEMENTclass[0]->NvnCs[P],dE);
			NvnCc[P] = pow(ELEMENTclass[0]->NvnCc[P],dE);
			NvnS[P]  = pow(ELEMENTclass[0]->NvnS[P],dE);

			size_t PbMin, PbMax;
			get_Pb_range(P,&PbMin,&PbMax);
			for (size_t Pb = PbMin; Pb <= PbMax; Pb++) {
				NvnIs[Pb] = pow(ELEMENTclass[0]->NvnIs[Pb],dE);
				NvnIc[Pb] = pow(ELEMENTclass[0]->NvnIc[Pb],dE);

				NfnGc[Pb][0] = pow(ELEMENTclass[0]->NvnGc[Pb],dE-1)*(ELEMENTclass[0]->NfnGc[Pb][0]);
				NfnS[Pb][0]  = pow(ELEMENTclass[0]->NvnS[Pb], dE-1)*(ELEMENTclass[0]->NfnS[Pb][0]);
				NfnIs[Pb][0] = pow(ELEMENTclass[0]->NvnIs[Pb],dE-1)*(ELEMENTclass[0]->NfnIs[Pb][0]);
				NfnIc[Pb][0] = pow(ELEMENTclass[0]->NvnIc[Pb],dE-1)*(ELEMENTclass[0]->NfnIc[Pb][0]);

				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->w_vIs[Pb],
				                  1,ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->w_fIs[Pb][0],NIn,NOut,OP,dE,0,Eclass);
				if (w_fIs[Pb][0] == NULL)
					w_fIs[Pb][0] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->w_vIc[Pb],
				                  1,ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->w_fIc[Pb][0],NIn,NOut,OP,dE,0,Eclass);
				if (w_fIc[Pb][0] == NULL)
					w_fIc[Pb][0] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->w_vIs[Pb],0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				if (w_vIs[Pb] == NULL)
					w_vIs[Pb] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->w_vIc[Pb],0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				if (w_vIc[Pb] == NULL)
					w_vIc[Pb] = sf_assemble_d(NIn,NOut,dE,OP);

				for (vh = 0; vh < Nvref; vh++) {
					if (vh == 0 || P == Pb) {
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],(double const *const *const) ELEMENTclass[0]->Ihat_vS_vS[P][Pb],
						                   0,0,NULL,NIn,NOut,OP,dE,vh,Eclass);
						Ihat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],(double const *const *const) ELEMENTclass[0]->L2hat_vS_vS[P][Pb],
						                   0,0,NULL,NIn,NOut,OP,dE,vh,Eclass);
						L2hat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (P == PGlobal && Pb == PGlobal) {
							get_sf_parametersV(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGs[1],(double const *const *const) ELEMENTclass[0]->I_vGs_vGs[1][1],
							                   0,0,NULL,NIn,NOut,OP,dE,vh,Eclass);
							I_vGs_vGs[1][1][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
					}
				}

				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->ChiS_vP[P][Pb][0],
				                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				ChiS_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
				                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				ChiS_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->I_vGc_vP[P][Pb][0],
				                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				I_vGc_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

				// Note: Most VOLUME operators need not interpolate between different orders
				if (P == Pb) {
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->ChiS_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnP[P],ELEMENTclass[0]->I_vGs_vP[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vP[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->ChiInvS_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->ChiBezInvS_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiBezInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGc[Pb],ELEMENTclass[0]->I_vGs_vGc[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vGc[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vCs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGs_vS[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vS[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vIs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vCc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGc_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCs_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCs_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCc_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCc_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (EFE) {
						get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
						                  ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
						                  NIn,NOut,OP,dE,3,Eclass);
						Is_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
						                  ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
						                  NIn,NOut,OP,dE,3,Eclass);
						Ic_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					} else {
						EXIT_UNSUPPORTED;
					}

					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->TS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					TS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->TS_vB[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					TS_vB[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->TInvS_vB[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					TInvS_vB[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					for (dim = 0; dim < dE; dim++) {
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],
						                  ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->D_vGs_vCs[1][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGs_vCs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],
						                  ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->D_vGs_vIs[1][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGs_vIs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
						                  ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->D_vGc_vCc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],
						                  ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->D_vGc_vIc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGc_vIc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->ICs[P][Pb][0],
						                  ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->D_vCs_vCs[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vCs_vCs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->ICc[P][Pb][0],
						                  ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->D_vCc_vCc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vCc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
						                  ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->GradChiS_vIs[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						GradChiS_vIs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
						                  ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->GradChiS_vIc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						GradChiS_vIc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (EFE) {
							get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
							                  ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ds_Weak_VV[P][Pb][0][0],
							                  NIn,NOut,OP,dE,dim,Eclass);
							Ds_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
							                  ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Dc_Weak_VV[P][Pb][0][0],
							                  NIn,NOut,OP,dE,dim,Eclass);
							Dc_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->Is_Strong_VV[P][Pb][0],
							                  ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->Ds_Strong_VV[P][Pb][0][0],
							                  NIn,NOut,OP,dE,dim,Eclass);
							Ds_Strong_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->Ic_Strong_VV[P][Pb][0],
							                  ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->Dc_Strong_VV[P][Pb][0][0],
							                  NIn,NOut,OP,dE,dim,Eclass);
							Dc_Strong_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							if (Collocated) {
								convert_to_CSR_d(NvnS[P],NvnIs[Pb],Ds_Weak_VV[P][Pb][0][dim],&Ds_Weak_VV_sp[P][Pb][0][dim]); // keep
								convert_to_CSR_d(NvnS[P],NvnIc[Pb],Dc_Weak_VV[P][Pb][0][dim],&Dc_Weak_VV_sp[P][Pb][0][dim]); // keep
							}
						} else {
							printf("Error: Unsupported EFE.\n"), EXIT_MSG;
						}
					}
				}

				for (f = 0; f < Nf; f++) {
				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],   (double const *const *const) ELEMENTclass[0]->ChiS_vS[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnS[Pb][0],(double const *const *const) ELEMENTclass[0]->ChiS_fS[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->ChiS_vIs[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const) ELEMENTclass[0]->ChiS_fIs[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->ChiS_vIc[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const) ELEMENTclass[0]->ChiS_fIc[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						convert_to_CSR_d(NfnIs[Pb][0],NvnS[P],ChiS_fIs[P][Pb][Vf],&ChiS_fIs_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NfnIc[Pb][0],NvnS[P],ChiS_fIc[P][Pb][Vf],&ChiS_fIc_sp[P][Pb][Vf]); // keep
					}

					if (P == Pb) {
						// Outside of fh == 0 if condition as operator is used for testing h-adaptation connectivity.
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnS[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGs_vS[1][Pb],
						                   ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnS[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGs_fS[1][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fS[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}

					get_sf_parametersF(ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->NvnS[P],(double const *const *const) ELEMENTclass[0]->Is_Weak_VV[P][Pb],
					                   ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->NvnS[P],(double const *const *const) ELEMENTclass[0]->Is_Weak_FV[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					Is_Weak_FV[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->NvnS[P],(double const *const *const) ELEMENTclass[0]->Ic_Weak_VV[P][Pb],
					                   ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->NvnS[P],(double const *const *const) ELEMENTclass[0]->Ic_Weak_FV[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					Ic_Weak_FV[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						convert_to_CSR_d(NvnS[P],NfnIs[Pb][0],Is_Weak_FV[P][Pb][Vf],&Is_Weak_FV_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NvnS[P],NfnIc[Pb][0],Ic_Weak_FV[P][Pb][Vf],&Ic_Weak_FV_sp[P][Pb][Vf]); // keep
					}

					if (fh == 0) {
						if (P == Pb) {
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGs_vIs[1][Pb],
							                   ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGs_fIs[1][Pb],
							                   NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIs[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGs_vIc[1][Pb],
							                   ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGs_fIc[1][Pb],
							                   NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIc[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnGc[Pb],   (double const *const *const) ELEMENTclass[0]->IGc[P][Pb],
							                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnGc[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGc_fGc[P][Pb],
							                   NIn,NOut,OP,dE,Vf,Eclass);
							I_vGc_fGc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							if (P == 2) {
								get_sf_parametersF(ELEMENTclass[0]->NvnG2[P],ELEMENTclass[0]->NvnG2[Pb],   (double const *const *const) ELEMENTclass[0]->IG2[P][Pb],
								                   ELEMENTclass[0]->NvnG2[P],ELEMENTclass[0]->NfnG2[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vG2_fG2[P][Pb],
								                   NIn,NOut,OP,dE,Vf,Eclass);
								I_vG2_fG2[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							}
						}
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGc_vS[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnS[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGc_fS[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGc_vIs[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGc_fIs[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGc_vIc[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGc_fIc[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnS[Pb],   (double const *const *const) ELEMENTclass[0]->I_vCs_vS[P][Pb],
						                   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnS[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vCs_fS[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->I_vCs_vIs[P][Pb],
						                   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vCs_fIs[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->I_vCs_vIc[P][Pb],
						                   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vCs_fIc[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnS[Pb],   (double const *const *const) ELEMENTclass[0]->I_vCc_vS[P][Pb],
						                   ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NfnS[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vCc_fS[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->I_vCc_vIs[P][Pb],
						                   ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vCc_fIs[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->I_vCc_vIc[P][Pb],
						                   ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vCc_fIc[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}

					for (dim = 0; dim < dE; dim++) {
						if (P == Pb) {
							get_sf_parametersFd(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGs_vIs[1][Pb],
						                        ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const *const) ELEMENTclass[0]->D_vGs_fIs[1][Pb],
						                        NIn,NOut,OP,dE,Vf,Eclass,dim,0);
							D_vGs_fIs[1][Pb][Vf][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							get_sf_parametersFd(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGs_vIc[1][Pb],
						                        ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const *const) ELEMENTclass[0]->D_vGs_fIc[1][Pb],
						                        NIn,NOut,OP,dE,Vf,Eclass,dim,0);
							D_vGs_fIc[1][Pb][Vf][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
						get_sf_parametersFd(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIs[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGc_vIs[P][Pb],
						                    ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIs[Pb][0],(double const *const *const *const) ELEMENTclass[0]->D_vGc_fIs[P][Pb],
						                    NIn,NOut,OP,dE,Vf,Eclass,dim,0);
						D_vGc_fIs[P][Pb][Vf][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersFd(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],   (double const *const *const) ELEMENTclass[0]->I_vGc_vIc[P][Pb],
						                    ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIc[Pb][0],(double const *const *const *const) ELEMENTclass[0]->D_vGc_fIc[P][Pb],
						                    NIn,NOut,OP,dE,Vf,Eclass,dim,0);
						D_vGc_fIc[P][Pb][Vf][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}
				}}

				for (e = 0; dE == DMAX && e < Ne; e++) {
					Ve = e*NEREFMAX;
					if (P == Pb) {
						get_sf_parametersE(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnGc[Pb],   (double const *const *const) ELEMENTclass[0]->IGc[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnGc[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vGc_fGc[P][Pb],
						                   NIn,NOut,OP,dE,Ve,Eclass);
						I_vGc_eGc[P][Pb][Ve] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (P == 2) {
							get_sf_parametersE(ELEMENTclass[0]->NvnG2[P],ELEMENTclass[0]->NvnG2[Pb],   (double const *const *const) ELEMENTclass[0]->IG2[P][Pb],
							                   ELEMENTclass[0]->NvnG2[P],ELEMENTclass[0]->NfnG2[Pb][0],(double const *const *const) ELEMENTclass[0]->I_vG2_fG2[P][Pb],
							                   NIn,NOut,OP,dE,Ve,Eclass);
							I_vG2_eG2[P][Pb][Ve] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
					}
				}

			}
		}
	} else if (Eclass == C_WEDGE) {
		unsigned int NIn0, NIn1, NOut0, NOut1;
		double       *OP0, *OP1, **OPF0, **OPF1;

		NvnGs[1] = (ELEMENTclass[0]->NvnGs[1])*(ELEMENTclass[1]->NvnGs[1]);

		for (size_t P = PSMin; P <= PSMax; P++) {
			NvnGc[P] = (ELEMENTclass[0]->NvnGc[P])*(ELEMENTclass[1]->NvnGc[P]);
			NvnCs[P] = (ELEMENTclass[0]->NvnCs[P])*(ELEMENTclass[1]->NvnCs[P]);
			NvnCc[P] = (ELEMENTclass[0]->NvnCc[P])*(ELEMENTclass[1]->NvnCc[P]);
			NvnS[P]  = (ELEMENTclass[0]->NvnS[P])*(ELEMENTclass[1]->NvnS[P]);

			size_t PbMin, PbMax;
			get_Pb_range(P,&PbMin,&PbMax);
			for (size_t Pb = PbMin; Pb <= PbMax; Pb++) {
				NvnIs[Pb] = (ELEMENTclass[0]->NvnIs[Pb])*(ELEMENTclass[1]->NvnIs[Pb]);
				NvnIc[Pb] = (ELEMENTclass[0]->NvnIc[Pb])*(ELEMENTclass[1]->NvnIc[Pb]);

				NfnS[Pb][0]  = (ELEMENTclass[0]->NfnS[Pb][0])* (ELEMENTclass[1]->NvnS[Pb]);
				NfnS[Pb][1]  = (ELEMENTclass[0]->NvnS[Pb])*    (ELEMENTclass[1]->NfnS[Pb][0]);
				NfnIs[Pb][0] = (ELEMENTclass[0]->NfnIs[Pb][0])*(ELEMENTclass[1]->NvnIs[Pb]);
				NfnIs[Pb][1] = (ELEMENTclass[0]->NvnIs[Pb])*   (ELEMENTclass[1]->NfnIs[Pb][0]);
				NfnIc[Pb][0] = (ELEMENTclass[0]->NfnIc[Pb][0])*(ELEMENTclass[1]->NvnIc[Pb]);
				NfnIc[Pb][1] = (ELEMENTclass[0]->NvnIc[Pb])*   (ELEMENTclass[1]->NfnIc[Pb][0]);

				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->w_vIs[Pb],
				                  1,ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->w_vIs[Pb],NIn,NOut,OP,dE,3,Eclass);
				if (w_vIs[Pb] == NULL)
					w_vIs[Pb] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->w_vIc[Pb],
				                  1,ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->w_vIc[Pb],NIn,NOut,OP,dE,3,Eclass);
				if (w_vIc[Pb] == NULL)
					w_vIc[Pb] = sf_assemble_d(NIn,NOut,dE,OP);

				for (vh = 0; vh < Nvref; vh++) {
					if (vh == 0 || P == Pb) {
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],(double const *const *const) ELEMENTclass[0]->Ihat_vS_vS[P][Pb],
						                   ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[Pb],(double const *const *const) ELEMENTclass[1]->Ihat_vS_vS[P][Pb],
						                   NIn,NOut,OP,dE,vh,Eclass);
						Ihat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],(double const *const *const) ELEMENTclass[0]->L2hat_vS_vS[P][Pb],
						                   ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[Pb],(double const *const *const) ELEMENTclass[1]->L2hat_vS_vS[P][Pb],
						                   NIn,NOut,OP,dE,vh,Eclass);
						L2hat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (P == PGlobal && Pb == PGlobal) {
							get_sf_parametersV(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGs[1],(double const *const *const) ELEMENTclass[0]->I_vGs_vGs[1][1],
							                   ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnGs[1],(double const *const *const) ELEMENTclass[1]->I_vGs_vGs[1][1],
							                   NIn,NOut,OP,dE,vh,Eclass);
							I_vGs_vGs[1][1][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
					}
				}

				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->ChiS_vP[P][Pb][0],
				                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnP[Pb],ELEMENTclass[1]->ChiS_vP[P][Pb][0],
				                  NIn,NOut,OP,dE,3,Eclass);
				ChiS_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
				                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->ChiS_vIc[P][Pb][0],
				                  NIn,NOut,OP,dE,3,Eclass);
				ChiS_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->I_vGc_vP[P][Pb][0],
				                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnP[Pb],ELEMENTclass[1]->I_vGc_vP[P][Pb][0],
				                  NIn,NOut,OP,dE,3,Eclass);
				I_vGc_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

				// Note: Most VOLUME operators need not interpolate between different orders
				if (P == Pb) {
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->ChiS_vS[P][Pb][0],
					                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[Pb],ELEMENTclass[1]->ChiS_vS[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					ChiS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnP[P],ELEMENTclass[0]->I_vGs_vP[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnP[P],ELEMENTclass[1]->I_vGs_vP[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vP[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->ChiInvS_vS[P][Pb][0],
					                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->ChiInvS_vS[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					ChiInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
					                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->ChiS_vIs[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGc[Pb],ELEMENTclass[0]->I_vGs_vGc[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnGc[Pb],ELEMENTclass[1]->I_vGs_vGc[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vGc[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnCs[Pb],ELEMENTclass[1]->I_vGs_vCs[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vCs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGs_vS[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->I_vGs_vS[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vS[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->I_vGs_vIs[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vIs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
					                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnCc[Pb],ELEMENTclass[1]->I_vGc_vCc[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vCc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGc_vS[P][Pb][0],
					                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->I_vGc_vS[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],
					                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->I_vGc_vIc[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCs_vIs[P][Pb][0],
					                  ELEMENTclass[1]->NvnCs[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->I_vCs_vIs[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCs_vIc[P][Pb][0],
					                  ELEMENTclass[1]->NvnCs[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->I_vCs_vIc[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCc_vIs[P][Pb][0],
					                  ELEMENTclass[1]->NvnCc[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->I_vCc_vIs[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCc_vIc[P][Pb][0],
					                  ELEMENTclass[1]->NvnCc[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->I_vCc_vIc[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (EFE) {
						get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
						                  ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->Is_Weak_VV[P][Pb][0],
						                  NIn,NOut,OP,dE,3,Eclass);
						Is_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
						                  ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->Ic_Weak_VV[P][Pb][0],
						                  NIn,NOut,OP,dE,3,Eclass);
						Ic_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					} else {
						printf("Error: Unsupported EFE.\n"), EXIT_MSG;
					}

					for (dim = 0; dim < dE; dim++) {
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGs_vCs[1][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGs_vCs[1][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],      OP1 = ELEMENTclass[1]->D_vGs_vCs[1][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],OP0,
						                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnCs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGs_vCs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGs_vIs[1][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGs_vIs[1][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],      OP1 = ELEMENTclass[1]->D_vGs_vIs[1][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],OP0,
						                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnIs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGs_vIs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGc_vCc[P][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGc_vCc[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],      OP1 = ELEMENTclass[1]->D_vGc_vCc[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],OP0,
						                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnCc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGc_vIc[P][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGc_vIc[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],      OP1 = ELEMENTclass[1]->D_vGc_vIc[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],OP0,
						                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnIc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGc_vIc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vCs_vCs[P][Pb][0][dim], OP1 = ELEMENTclass[1]->ICs[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->ICs[P][Pb][0],            OP1 = ELEMENTclass[1]->D_vCs_vCs[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnCs[Pb],OP0,
						                  ELEMENTclass[1]->NvnCs[P],ELEMENTclass[1]->NvnCs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vCs_vCs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vCc_vCc[P][Pb][0][dim], OP1 = ELEMENTclass[1]->ICc[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->ICc[P][Pb][0],            OP1 = ELEMENTclass[1]->D_vCc_vCc[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnCc[Pb],OP0,
						                  ELEMENTclass[1]->NvnCc[P],ELEMENTclass[1]->NvnCc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vCc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (EFE) {
							if (dim < 2) OP0 = ELEMENTclass[0]->Ds_Weak_VV[P][Pb][0][dim], OP1 = ELEMENTclass[1]->Is_Weak_VV[P][Pb][0];
							else         OP0 = ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],      OP1 = ELEMENTclass[1]->Ds_Weak_VV[P][Pb][0][0];
							get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],OP0,
							                  ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->NvnS[P],OP1,NIn,NOut,OP,dE,3,Eclass);
							Ds_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							if (dim < 2) OP0 = ELEMENTclass[0]->Dc_Weak_VV[P][Pb][0][dim], OP1 = ELEMENTclass[1]->Ic_Weak_VV[P][Pb][0];
							else         OP0 = ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],      OP1 = ELEMENTclass[1]->Dc_Weak_VV[P][Pb][0][0];
							get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],OP0,
							                  ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->NvnS[P],OP1,NIn,NOut,OP,dE,3,Eclass);
							Dc_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							if (Collocated) {
								convert_to_CSR_d(NvnS[P],NvnIs[Pb],Ds_Weak_VV[P][Pb][0][dim],&Ds_Weak_VV_sp[P][Pb][0][dim]); // keep
								convert_to_CSR_d(NvnS[P],NvnIc[Pb],Dc_Weak_VV[P][Pb][0][dim],&Dc_Weak_VV_sp[P][Pb][0][dim]); // keep
							}
						} else {
							;
						}
					}
				}

				for (f = 0; f < Nf; f++) {
				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;

					if (f < 3) { OPF0  = ELEMENTclass[0]->ChiS_fS[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_vS[P][Pb];
					             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],    NOut1 = ELEMENTclass[1]->NvnS[Pb];
					} else {     OPF0  = ELEMENTclass[0]->ChiS_vS[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_fS[P][Pb];
					             NOut0 = ELEMENTclass[0]->NvnS[Pb],       NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],NOut0,(double const *const *const) OPF0,
					                   ELEMENTclass[1]->NvnS[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					if (f < 3) { OPF0  = ELEMENTclass[0]->ChiS_fIs[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_vIs[P][Pb];
					             NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],    NOut1 = ELEMENTclass[1]->NvnIs[Pb];
					} else {     OPF0  = ELEMENTclass[0]->ChiS_vIs[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_fIs[P][Pb];
					             NOut0 = ELEMENTclass[0]->NvnIs[Pb],       NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],NOut0,(double const *const *const) OPF0,
					                   ELEMENTclass[1]->NvnS[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					if (f < 3) { OPF0  = ELEMENTclass[0]->ChiS_fIc[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_vIc[P][Pb];
					             NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],    NOut1 = ELEMENTclass[1]->NvnIc[Pb];
					} else {     OPF0  = ELEMENTclass[0]->ChiS_vIc[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_fIc[P][Pb];
					             NOut0 = ELEMENTclass[0]->NvnIc[Pb],       NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],NOut0,(double const *const *const) OPF0,
					                   ELEMENTclass[1]->NvnS[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						if (f < 3) IndClass = 0;
						else       IndClass = 1;

						convert_to_CSR_d(NfnIs[Pb][IndClass],NvnS[P],ChiS_fIs[P][Pb][Vf],&ChiS_fIs_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NfnIc[Pb][IndClass],NvnS[P],ChiS_fIc[P][Pb][Vf],&ChiS_fIc_sp[P][Pb][Vf]); // keep
					}

					if (P == Pb) {
						// Outside of fh == 0 if condition as operator is used for testing h-adaptation connectivity.
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fS[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vS[1][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGs_vS[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fS[1][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnGs[1],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fS[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}

					if (f < 3) { OPF0 = ELEMENTclass[0]->Is_Weak_FV[P][Pb], OPF1 = ELEMENTclass[1]->Is_Weak_VV[P][Pb];
					             NIn0 = ELEMENTclass[0]->NfnIs[Pb][0],      NIn1 = ELEMENTclass[1]->NvnIs[Pb];
					} else {     OPF0 = ELEMENTclass[0]->Is_Weak_VV[P][Pb], OPF1 = ELEMENTclass[1]->Is_Weak_FV[P][Pb];
					             NIn0 = ELEMENTclass[0]->NvnIs[Pb],         NIn1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
					get_sf_parametersF(NIn0,ELEMENTclass[0]->NvnS[P],(double const *const *const) OPF0,
					                   NIn1,ELEMENTclass[1]->NvnS[P],(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					Is_Weak_FV[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					if (f < 3) { OPF0 = ELEMENTclass[0]->Ic_Weak_FV[P][Pb], OPF1 = ELEMENTclass[1]->Ic_Weak_VV[P][Pb];
					             NIn0 = ELEMENTclass[0]->NfnIc[Pb][0],      NIn1 = ELEMENTclass[1]->NvnIc[Pb];
					} else {     OPF0 = ELEMENTclass[0]->Ic_Weak_VV[P][Pb], OPF1 = ELEMENTclass[1]->Ic_Weak_FV[P][Pb];
					             NIn0 = ELEMENTclass[0]->NvnIc[Pb],         NIn1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
					get_sf_parametersF(NIn0,ELEMENTclass[0]->NvnS[P],(double const *const *const) OPF0,
					                   NIn1,ELEMENTclass[1]->NvnS[P],(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					Ic_Weak_FV[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						if (f < 3) IndClass = 0;
						else       IndClass = 1;

						convert_to_CSR_d(NvnS[P],NfnIs[Pb][IndClass],Is_Weak_FV[P][Pb][Vf],&Is_Weak_FV_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NvnS[P],NfnIc[Pb][IndClass],Ic_Weak_FV[P][Pb][Vf],&Ic_Weak_FV_sp[P][Pb][Vf]); // keep
					}

					if (fh == 0) {
						if (P == Pb) {
							if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fIs[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vIs[1][Pb];
							             NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
							} else {     OPF0  = ELEMENTclass[0]->I_vGs_vIs[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fIs[1][Pb];
							             NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],NOut0,(double const *const *const) OPF0,
							                   ELEMENTclass[1]->NvnGs[1],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIs[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fIc[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vIc[1][Pb];
							             NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
							} else {     OPF0  = ELEMENTclass[0]->I_vGs_vIc[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fIc[1][Pb];
							             NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],NOut0,(double const *const *const) OPF0,
							                   ELEMENTclass[1]->NvnGs[1],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIc[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGc_fS[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_vS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGc_vS[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_fS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnGc[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGc_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_vIs[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGc_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_fIs[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnGc[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGc_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_vIc[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGc_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_fIc[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnGc[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCs_fS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_vS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCs_vS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_fS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnCs[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCs_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_vIs[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCs_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_fIs[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnCs[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCs_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_vIc[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCs_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_fIc[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnCs[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCc_fS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_vS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCc_vS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_fS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnCc[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCc_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_vIs[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCc_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_fIs[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnCc[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCc_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_vIc[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCc_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_fIc[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],NOut0,(double const *const *const) OPF0,
						                   ELEMENTclass[1]->NvnCc[P],NOut1,(double const *const *const) OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}
				}}
			}
		}
	} else {
		printf("Error: Unsupported Eclass.\n"), EXIT_MSG;
	}
	free(ones_Nf);

	compute_ELEMENT_Volume(EType);
}

static void setup_L2_projection_preoperators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnIc, *NvnS;
	double       **w_vIc, ****ChiS_vIc;

	// Initialize DB Parameters
	unsigned int PGs    = DB.PGs,
	             NP     = DB.NP;

	char         *BasisType     = DB.BasisType,
	             ***NodeTypeIvc = DB.NodeTypeIvc,
	             ***NodeTypeS   = DB.NodeTypeS;

	// Standard datatypes
	unsigned int i, iMax, dE, vh, PIvc[NEC],
	             Nve, Nbf, Eclass, Nvref, NEhref, Indh, *Nvve, *EType_h;
	double       *E_rst_vV, *rst_vV, **VeV, **rst_vIc, *rst_vS,
	             *IGs, *IS, *TS, *ChiRefS_vS, *ChiRefInvGs_vGs, *ChiRefInvS_vS,
	             *ChiRefGs_vGs, *ChiRefGs_vIc, *ChiRefS_vIc, *ChiS_vS;

	struct BCoords {
		double **Ic;
	} **BCoords_V;
	struct S_ELEMENT *ELEMENT, *ELEMENT_h;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	E_rst_vV = ChiS_vS = NULL;

	if (!(EType == TET || EType == PYR))
		printf("Error: Unsupported EType.\n"), EXIT_MSG;

	PIvc[1] = PIvcMaxTET;
	PIvc[2] = PIvcMaxPYR;

	ELEMENT   = get_ELEMENT_type(EType);
	ELEMENT_h = ELEMENT;

	Eclass = get_Eclass(EType);
	setup_ELEMENT_VeV(EType);

	dE      = ELEMENT->d;
	Nve     = ELEMENT->Nve;
	VeV     = ELEMENT->VeV;
	Nvve    = ELEMENT->Nvve;
	Nvref   = ELEMENT->Nvref;
	NEhref  = ELEMENT->NEhref;
	EType_h = ELEMENT->type_h;

	select_functions(&basis,&grad_basis,&cubature,EType);

	// Stored operators
	NvnIc    = ELEMENT->NvnIc;
	NvnS     = ELEMENT->NvnS;
	w_vIc    = ELEMENT->w_vIc;
	ChiS_vIc = ELEMENT->ChiS_vIc;

	// Allocate memory for arrays with multiple levels of dereferencing
	rst_vIc = malloc(NVREFMAX * sizeof *rst_vIc); // free

	BCoords_V = malloc(NEhref * sizeof *BCoords_V); // free
	for (i = 0, iMax = NEhref; i < iMax; i++) {
		BCoords_V[i]     = malloc(     sizeof *BCoords_V[i]);       // free
		BCoords_V[i]->Ic = calloc(NP , sizeof *(BCoords_V[i]->Ic)); // free
	}

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	// Nve + 1 for TETrefineType == TET6
	rst_vV = malloc((Nve+1)*dE * sizeof *rst_vV); // free

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	for (size_t P = PSMin; P <= PSMax; P++) {
		set_cubdata(CUBDATA,false,false,NodeTypeS[P][Eclass],dE,P,cubature); // free
		set_from_cubdata(CUBDATA,&NvnS[P],NULL,&rst_vS,NULL,NULL);

		IS         = identity_d(NvnS[P]);             // free
		ChiRefS_vS = basis(P,rst_vS,NvnS[P],&Nbf,dE); // free
		free(rst_vS);
		if (strstr(BasisType,"Modal")) {
			ChiS_vS = ChiRefS_vS;
		} else if (strstr(BasisType,"Nodal")) {
			ChiS_vS = IS;
		}

		ChiRefInvS_vS = inverse_d(NvnS[P],NvnS[P],ChiRefS_vS,IS);                                     // free
		TS            = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P],NvnS[P],NvnS[P],1.0,ChiRefInvS_vS,ChiS_vS); // free

		free(IS);
		free(ChiRefS_vS);
		free(ChiRefInvS_vS);

		size_t PbMin, PbMax;
		get_Pb_range(P,&PbMin,&PbMax);
		for (size_t Pb = PbMin; Pb <= PbMax; Pb++) {
			if (w_vIc[Pb])
				free(w_vIc[Pb]);
			set_cubdata(CUBDATA,true,false,NodeTypeIvc[Pb][Eclass],dE,PIvc[Eclass],cubature); // free
			set_from_cubdata(CUBDATA,&NvnIc[Pb],NULL,&rst_vIc[0],&w_vIc[Pb],NULL);
			free(rst_vIc[0]);

			for (i = iMax = NEhref; i--; ) {
				ELEMENT_h = get_ELEMENT_type(EType_h[i]);

				Eclass = get_Eclass(EType_h[i]);
				select_functions(&basis,&grad_basis,&cubature,EType_h[i]);

				Nve   = ELEMENT_h->Nve;
				NvnIc = ELEMENT_h->NvnIc;

				set_cubdata(CUBDATA,false,false,NodeTypeIvc[Pb][Eclass],dE,PIvc[Eclass],cubature); // free
				set_from_cubdata(CUBDATA,&NvnIc[Pb],NULL,&rst_vIc[0],NULL,NULL);

				E_rst_vV        = get_rst_vV(ELEMENT_h);                   // free
				IGs             = identity_d(Nve);                         // free
				ChiRefGs_vGs    = basis(PGs,E_rst_vV,Nve,&Nbf,dE);         // free
				ChiRefInvGs_vGs = inverse_d(Nve,Nve,ChiRefGs_vGs,IGs);     // free
				ChiRefGs_vIc    = basis(PGs,rst_vIc[0],NvnIc[Pb],&Nbf,dE); // free

				BCoords_V[i]->Ic[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIc[Pb],Nve,Nve,1.0,ChiRefGs_vIc,ChiRefInvGs_vGs); // free

				free(IGs);
				free(ChiRefGs_vGs);
				free(ChiRefInvGs_vGs);
				free(ChiRefGs_vIc);

				if (i) {
					free(E_rst_vV);
					free(rst_vIc[0]);
				}
			}

			for (vh = 0; vh < Nvref; vh++) {
				mm_CTN_d(Nvve[vh],dE,Nve,VeV[vh],E_rst_vV,rst_vV);
				if (vh) {
					Indh = get_IndEhref(EType,vh);
					if (EType == TET || EType == PYR)
						ELEMENT_h = get_ELEMENT_type(EType_h[Indh]);
					else
						ELEMENT_h = ELEMENT;

					rst_vIc[vh] = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],dE,Nvve[vh],1.0,BCoords_V[Indh]->Ic[Pb],rst_vV); // free
				}

				ChiRefS_vIc         = basis(P,rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE);                                  // free
				ChiS_vIc[P][Pb][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIc,TS); // keep

				free(ChiRefS_vIc);
			}

			free(E_rst_vV);
			for (vh = 0; vh < Nvref; vh++)
				free(rst_vIc[vh]);

			for (i = 0; i < NEhref; i++)
				free(BCoords_V[i]->Ic[Pb]);
		}
		free(TS);
	}
	free(rst_vV);
	free(rst_vIc);
	for (i = 0; i < NEhref; i++) {
		free(BCoords_V[i]->Ic);
		free(BCoords_V[i]);
	}
	free(BCoords_V);

	free(CUBDATA);
}

static void setup_L2_projection_operators(const unsigned int EType)
{
	/*
	 *	Comments:
	 *		This function is separated from the setup_ELEMENT_operators function as the order P operators defined here
	 *		rely on order Pb operators for their assembly.
	 *		L2 projection operators are only needed for coarsening operations. During refinement, these operators reduce
	 *		to standard interpolations.
	 *		For PYR ELEMENTs, the order used for the cubature nodes has a signifiant impact on the accuracy of the L2
	 *		projections (several orders of magnitude), which remain quite poor even when using the highest order
	 *		cubature nodes supported (TET: P10, PYR: P6). As the cubature order is not seen in the L2 projection
	 *		operator, it is advantageous to compute it using the highest possible cubature order. However, this order is
	 *		not required for other operators in the code and is thus redefined appropriately in setup_operators.
	 */

	// Returned Operators
	double ****L2hat_vS_vS;

	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int vh,
	             Indh, Nvref, *NvnS, *NvnIc, *EType_h;
	double       L2_scale, *IS, *ChiTW, *M, *MInv, *S, *diag_w_vIc,
	             ****ChiS_vIc, **w_vIc;

	struct S_ELEMENT *ELEMENT, *ELEMENT_h;

	// silence
	ELEMENT_h = NULL;

	if (Adapt == ADAPT_0)
		return;

	ELEMENT = get_ELEMENT_type(EType);

	// Stored operators
	L2hat_vS_vS = ELEMENT->L2hat_vS_vS;

	EType_h = ELEMENT->type_h;

	Nvref = ELEMENT->Nvref;
	NvnS  = ELEMENT->NvnS;
	NvnIc = ELEMENT->NvnIc;

	ChiS_vIc = ELEMENT->ChiS_vIc;
	w_vIc    = ELEMENT->w_vIc;

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	for (size_t P = PSMin; P <= PSMax; P++) {
		size_t PbMin, PbMax;
		get_Pb_range(P,&PbMin,&PbMax);
		for (size_t Pb = PbMin; Pb <= PbMax; Pb++) {
			diag_w_vIc = diag_d(w_vIc[Pb],NvnIc[Pb]); // free

			IS    = identity_d(NvnS[Pb]);                                                                       // free
			ChiTW = mm_Alloc_d(CBRM,CBT, CBNT,NvnS[Pb],NvnIc[Pb],NvnIc[Pb],1.0,ChiS_vIc[Pb][Pb][0],diag_w_vIc); // free
			M     = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[Pb],NvnS[Pb], NvnIc[Pb],1.0,ChiTW,ChiS_vIc[Pb][Pb][0]);      // free
			MInv  = inverse_d(NvnS[Pb],NvnS[Pb],M,IS);                                                          // free

			free(diag_w_vIc);
			free(IS);
			free(ChiTW);
			free(M);
			for (vh = 0; vh < Nvref; vh++) {
				if (vh == 0 || P == Pb) {
					Indh = get_IndEhref(EType,vh);
					if (EType == TET || EType == PYR)
						ELEMENT_h = get_ELEMENT_type(EType_h[Indh]);
					else
						ELEMENT_h = ELEMENT;

					diag_w_vIc = diag_d(ELEMENT_h->w_vIc[Pb],ELEMENT_h->NvnIc[Pb]); // free

					ChiTW = mm_Alloc_d(CBRM,CBT, CBNT,NvnS[Pb],ELEMENT_h->NvnIc[Pb],ELEMENT_h->NvnIc[Pb],1.0,ChiS_vIc[Pb][Pb][vh],diag_w_vIc);     // free
					S     = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[Pb],ELEMENT_h->NvnS[P],  ELEMENT_h->NvnIc[Pb],1.0,ChiTW,ELEMENT_h->ChiS_vIc[P][Pb][0]); // free

					L2_scale = get_L2_scaling(EType,vh);

					// Returned Operators
					L2hat_vS_vS[P][Pb][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[Pb],ELEMENT_h->NvnS[P],NvnS[Pb],L2_scale,MInv,S); // keep

					free(diag_w_vIc);
					free(ChiTW);
					free(S);
				}
			}
			free(MInv);
		}
	}


/*
// ToBeDeleted
	// Returned Operators
	double ****L2hat_vS_vS,
	       ****GfS_fIs, ****GfS_fIc; // Change G to L2 if these operators are kept (ToBeDeleted)

	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int P, Pb, vref, f, fh, PSMin, PSMax, PbMin, PbMax, Nf, fhMax,
	             Vf, PM, *Nfref, Eclass, IndFType, Nvref,
	             NvnS, NvnSN, NvnI, NfnM, NfnF, NfnIs, NfnIc,
	             *ones_Nf;
	double       L2_scale, *ChiSN_vI, *wSN_vI, *diag_wSN_vI, *ISN, *ChiTW_SN, *MSN, *MSNInv, *ChiS_vI, *SSNS,
	             *one_d, *ChiS_vF, *ChiM_vI, *ChiF_vIs, *ChiF_vIc,
	             *ChiInvF_vF, *ChiInvM_vM, *wM_vI, *diag_wM_vI, *IF, *IM, *ChiTW_M, *ChiTW_F, *ChiTW_SNr,
	             *SMF, *SFM, *MM, *MF, *MMInv, *MFInv, *MInvSMF, *MInvSFM,
	             *IhatS_fS, *L2hatvS_fS, *L2hatfS_fI;

	struct S_ELEMENT *ELEMENT, *ELEMENT_F;

	if (Adapt == ADAPT_0)
		return;

	// silence
	IF = IM = NULL;
	diag_wM_vI = ChiTW_M = NULL;
	MM = MMInv = NULL;

	one_d = malloc(1 * sizeof *one_d); // free
	one_d[0] = 1.0;

	ELEMENT = get_ELEMENT_type(EType);
	Eclass = get_Eclass(EType);

	Nf    = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;
	Nvref = ELEMENT->Nvref;

	// Stored operators
	L2hat_vS_vS = ELEMENT->L2hat_vS_vS;
	GfS_fIs     = ELEMENT->GfS_fIs;
	GfS_fIc     = ELEMENT->GfS_fIc;


	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		NvnS = ELEMENT->NvnS[P];

		get_Pb_range(P,&PbMin,&PbMax);
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			PM = max(P,Pb);

			NvnSN = ELEMENT->NvnS[Pb]; // N: New
			NvnI  = ELEMENT->NvnIc[Pb];

			ChiSN_vI = ELEMENT->ChiS_vIc[Pb][Pb][0];

			wSN_vI      = ELEMENT->w_vIc[Pb];
			diag_wSN_vI = diag_d(wSN_vI,NvnI); // free

			ISN      = identity_d(NvnSN); // free
			ChiTW_SN = mm_Alloc_d(CBRM,CBT, CBNT,NvnSN,NvnI, NvnI,1.0,ChiSN_vI,diag_wSN_vI); // free
			MSN      = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSN,NvnSN,NvnI,1.0,ChiTW_SN,ChiSN_vI);    // free
			MSNInv   = inverse_d(NvnSN,NvnSN,MSN,ISN);                                       // free

			ChiS_vI  = ELEMENT->ChiS_vIc[P][Pb][0];
			for (vref = 0; vref < Nvref; vref++) {
				ChiSN_vI = ELEMENT->ChiS_vIc[Pb][Pb][vref];

				ChiTW_SNr = mm_Alloc_d(CBRM,CBT, CBNT,NvnSN,NvnI,NvnI,1.0,ChiSN_vI, diag_wSN_vI); // free
				SSNS      = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSN,NvnS,NvnI,1.0,ChiTW_SNr,ChiS_vI);     // free

				L2_scale = get_L2_scaling(EType,vref);

				// Returned Operators
				L2hat_vS_vS[P][Pb][vref] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSN,NvnS,NvnSN,L2_scale,MSNInv,SSNS); // keep

				free(ChiTW_SNr);
				free(SSNS);
			}

			free(diag_wSN_vI);
			free(ISN);
			free(ChiTW_SN);
			free(MSN);
			free(MSNInv);
		}

		for (Pb = P; Pb <= PbMax; Pb++) {
			PM = max(P,Pb);
			for (f = 0; f < Nf; f++) {
				IndFType = get_IndFType(Eclass,f);
				ELEMENT_F = get_ELEMENT_FACE(EType,IndFType);

				if (EType == LINE) {
					NfnM = 1; // (M)ortar
					NfnF = 1; // (F)ace
					NfnIs = 1;
					NfnIc = 1;

//					ChiM_vS   = one_d;
				} else {
					NfnM = ELEMENT_F->NvnS[PM];
					NfnF = ELEMENT_F->NvnS[Pb];
					NfnIs = ELEMENT_F->NvnIs[PM];
					NfnIc = ELEMENT_F->NvnIc[PM];

//					ChiM_vS  = ELEMENT_F->ChiS_vS[PM][PM][0];
					ChiM_vI  = ELEMENT_F->ChiS_vIc[PM][PM][0];
					wM_vI    = ELEMENT_F->w_vIc[PM];

					ChiInvF_vF = ELEMENT_F->ChiInvS_vS[Pb][Pb][0];
					ChiInvM_vM = ELEMENT_F->ChiInvS_vS[PM][PM][0];

					IF   = identity_d(NfnF); // free
					IM   = identity_d(NfnM); // free

					diag_wM_vI = diag_d(wM_vI,NfnIc); // free
					ChiTW_M    = mm_Alloc_d(CBRM,CBT,CBNT,NfnM,NfnIc,NfnIc,1.0,ChiM_vI,diag_wM_vI); // free
					MM         = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NfnM,NfnIc,1.0,ChiTW_M,ChiM_vI);    // free
					MMInv      = inverse_d(NfnM,NfnM,MM,IM);                                        // free
				}

				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;
//printf("P, Pb, f, fh, Vf: %d %d %d %d %d\n",P,Pb,f,fh,Vf);

					ChiS_vF = ELEMENT->ChiS_fS[P][Pb][f*NFREFMAX];

					if (EType == LINE) {
						ChiF_vIs  = one_d;
						ChiF_vIc  = one_d;
						L2hatvS_fS = mm_Alloc_d(CBRM,CBNT,CBNT,1,NvnS,1,1.0,one_d,ChiS_vF); // free
						L2hatfS_fI = mm_Alloc_d(CBRM,CBNT,CBNT,1,1,1,1.0,one_d,one_d);      // free
					} else {
						ChiF_vIs = ELEMENT_F->ChiS_vIs[Pb][PM][fh];
						ChiF_vIc = ELEMENT_F->ChiS_vIc[Pb][PM][fh];

						ChiTW_F  = mm_Alloc_d(CBRM,CBT,CBNT,NfnF,NfnIc,NfnIc,1.0,ChiF_vIc,diag_wM_vI); // free

						SMF      = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NfnF,NfnIc,1.0,ChiTW_M,ChiF_vIc);  // free
						SFM      = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnM,NfnIc,1.0,ChiTW_M,ChiF_vIc);  // free
						MF       = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnF,NfnIc,1.0,ChiTW_F,ChiF_vIc);  // free
						MFInv    = inverse_d(NfnF,NfnF,MF,IF);                                       // free
						MInvSMF  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NfnF,NfnM,1.0,MMInv,SMF);          // free
						MInvSFM  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnM,NfnF,1.0,MFInv,SFM);          // free
						IhatS_fS = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NvnS,NfnF,1.0,ChiInvF_vF,ChiS_vF); // free

						L2hatvS_fS = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NvnS,NfnF,1.0,MInvSMF,IhatS_fS);   // free
						L2hatfS_fI = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnM,NfnM,1.0,MInvSFM,ChiInvM_vM); // free

						free(ChiTW_F);
						free(SMF), free(SFM), free(MF);
						free(MFInv), free(MInvSMF), free(MInvSFM);
						free(IhatS_fS);
					}

					// Returned Operators
					GfS_fIs[P][Pb][Vf]   = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs,NfnM,NfnF,1.0,ChiF_vIs,L2hatfS_fI); // keep
					GfS_fIc[P][Pb][Vf]   = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc,NfnM,NfnF,1.0,ChiF_vIc,L2hatfS_fI); // keep

					free(L2hatvS_fS);
					free(L2hatfS_fI);
				}
				free(IF), free(IM);
				free(diag_wM_vI), free(ChiTW_M);
				free(MM);
				free(MMInv);
			}
		}
	}
	free(one_d);
	free(ones_Nf);
*/
}

static void setup_blending(const unsigned int EType)
{
	// Returned operators
	double ****I_fGs_vG2, ****I_fG2_vG2,
	       ****I_fGs_vGc, ****I_fGc_vGc, ****I_eGs_vGc, ****I_eGc_vGc;

	// Initialize DB Parameters
	unsigned int PGs      = DB.PGs,
	             *PGc     = DB.PGc,
	             Blending = DB.Blending;

	char         **NodeTypeG = DB.NodeTypeG;

	// Standard datatypes
	unsigned int f, e, n, ve, P, Pb, Vf, Ve, dim, PbMin, PbMax, Nf, Ne, Nve, *Nfve, Neve, *VeFcon, *VeEcon,
	             *NvnGc, *Nv0nGs, *Nv0nGc, *NvnG2, *Nv0nG2,
	             Eclass, EclassF, EclassE, dE, Nbf, dimF[2], dimE;
	double       *BCoords_V, *BCoords_V2, *BCoords_F, *BCoords_E,
	             *rst_v0Gs, *rst_proj, *rst_proj2, *rst_v0G2, *rst_v0Gc, *rst_vGs, *rst_vGsF, *rst_vGsE, rst_Num, rst_Den,
	             *ChiRefGs_vGs, *ChiRefGc_vGc, *ChiRefGs_vProj, *ChiRefGc_vProj, *ChiRefG2_vG2, *ChiRefGs_vProj2, *ChiRefG2_vProj2,
	             *ChiRefInvGs_vGs, *IGs, *ChiRefInvGc_vGc, *IGc, *IG2, *ChiRefInvG2_vG2;

	struct S_ELEMENT *ELEMENT, *ELEMENT_F, *ELEMENT_E;

	// Function pointers
	cubature_tdef cubature;
	basis_tdef    basis;

	ELEMENT = get_ELEMENT_type(EType);

	dE = ELEMENT->d;
	Eclass = get_Eclass(EType);

	Ne     = ELEMENT->Ne;
	Nf     = ELEMENT->Nf;
	Nve    = ELEMENT->Nve;
	Neve   = ELEMENT->Neve;
	Nfve   = ELEMENT->Nfve;
	VeEcon = ELEMENT->VeEcon;
	VeFcon = ELEMENT->VeFcon;
	NvnG2  = ELEMENT->NvnG2;
	NvnGc  = ELEMENT->NvnGc;

	I_fGs_vG2 = ELEMENT->I_fGs_vG2;
	I_fG2_vG2 = ELEMENT->I_fG2_vG2;
	I_fGs_vGc = ELEMENT->I_fGs_vGc;
	I_fGc_vGc = ELEMENT->I_fGc_vGc;

	I_eGs_vGc = ELEMENT->I_eGs_vGc;
	I_eGc_vGc = ELEMENT->I_eGc_vGc;

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		BCoords_V  = ELEMENT->I_vGs_vGc[1][P][0];
		BCoords_V2 = ELEMENT->I_vGs_vG2[1][2][0];

		// FACE related operators
		for (f = 0; f < Nf; f++) {
			Vf = f*NFREFMAX;

			ELEMENT_F = get_ELEMENT_F_type(EType,f);
			EclassF   = get_Eclass(ELEMENT_F->type);

			rst_v0Gs  = get_rst_vV(ELEMENT_F); // free
			rst_proj  = malloc(NvnGc[P]*(dE-1) * sizeof *rst_proj);  // free
			rst_proj2 = malloc(NvnG2[2]*(dE-1) * sizeof *rst_proj2); // free

			if (Blending == SCOTT && EType == TRI) {
				// dE == 2
				for (dim = 0; dim < dE-1; dim++) {
					for (n = 0; n < NvnGc[P]; n++)
						rst_proj[n+dim*NvnGc[P]] = 2.0*BCoords_V[n*Nve+VeFcon[f*NFVEMAX+1]]-1.0;
					for (n = 0; n < NvnG2[2]; n++)
						rst_proj2[n+dim*NvnG2[2]] = 2.0*BCoords_V2[n*Nve+VeFcon[f*NFVEMAX+1]]-1.0;
				}
			} else if (Blending == NIELSON && EType == TRI) {
				// dE == 2
				for (dim = 0; dim < dE-1; dim++) {
					for (n = 0; n < NvnGc[P]; n++) {
						rst_Num = BCoords_V[n*Nve+VeFcon[f*NFVEMAX+1]]-BCoords_V[n*Nve+VeFcon[f*NFVEMAX+0]];

						if (fabs(rst_Num) < EPS) {
							rst_proj[n+dim*NvnGc[P]] = 0.0;
							continue;
						}

						rst_Den = 1.0-BCoords_V[n*Nve+f];
						rst_proj[n+dim*NvnGc[P]] = rst_Num/rst_Den;
					}
					for (n = 0; n < NvnG2[2]; n++) {
						rst_Num = BCoords_V2[n*Nve+VeFcon[f*NFVEMAX+1]]-BCoords_V2[n*Nve+VeFcon[f*NFVEMAX+0]];

						if (fabs(rst_Num) < EPS) {
							rst_proj2[n+dim*NvnG2[2]] = 0.0;
							continue;
						}

						rst_Den = 1.0-BCoords_V2[n*Nve+f];
						rst_proj2[n+dim*NvnG2[2]] = rst_Num/rst_Den;
					}
				}
			} else if (Eclass == C_SI) { // Default to SZABO_BABUSKA blending for simplex elements
				// These BCoords_F determine rst_proj according to Szabo(1991, p. 111) in 2D  (i.e. xi = L2-L1).
				BCoords_F = malloc(NvnGc[P]*Nfve[f] * sizeof *BCoords_F); // free

				for (n = 0; n < NvnGc[P]; n++) {
				for (ve = 0; ve < Nfve[f]; ve++) {
					BCoords_F[n*Nfve[f]+ve] = BCoords_V[n*Nve+VeFcon[f*NFVEMAX+ve]];
				}}

				// Compute projection of VOLUME nodes to FACE
				mm_CTN_d(NvnGc[P],dE-1,Nfve[f],BCoords_F,rst_v0Gs,rst_proj);
				free(BCoords_F);

				BCoords_F = malloc(NvnG2[2]*Nfve[f] * sizeof *BCoords_F); // free

				for (n = 0; n < NvnG2[2]; n++) {
				for (ve = 0; ve < Nfve[f]; ve++) {
					BCoords_F[n*Nfve[f]+ve] = BCoords_V2[n*Nve+VeFcon[f*NFVEMAX+ve]];
				}}

				// Compute projection of VOLUME nodes to FACE
				mm_CTN_d(NvnG2[2],dE-1,Nfve[f],BCoords_F,rst_v0Gs,rst_proj2);
				free(BCoords_F);
			} else if (Eclass == C_TP) { // GORDON_HALL blending
				rst_vGs = get_rst_vV(ELEMENT); // free
				rst_vGsF = malloc(Nve*(dE-1) * sizeof *rst_vGsF); // free

				if      (f/2 == 0) { dimF[0] = 1; dimF[1] = 2; }
				else if (f/2 == 1) { dimF[0] = 0; dimF[1] = 2; }
				else if (f/2 == 2) { dimF[0] = 0; dimF[1] = 1; }

				for (dim = 0; dim < dE-1; dim++) {
					for (n = 0; n < Nve; n++)
						rst_vGsF[n+dim*Nve] = rst_vGs[n+dimF[dim]*Nve];
				}
				free(rst_vGs);

				mm_CTN_d(NvnGc[P],dE-1,Nve,BCoords_V,rst_vGsF,rst_proj);
				mm_CTN_d(NvnG2[2],dE-1,Nve,BCoords_V2,rst_vGsF,rst_proj2);
				free(rst_vGsF);
			} else {
				// Allow correct_ToBeCurved in setup_ToBeCurved for WEDGE and PYR elements. (ToBeDeleted)
				printf("Add Support.\n"), EXIT_MSG;
			}
if (P == 4 && EType == TRI && f == 2) {
//	printf("rst_proj:\n");
//	array_print_d(NvnGc[P],dE-1,rst_proj,'C');
}

			// Compute blending operators
			Nv0nGs = ELEMENT_F->NvnGs;
			Nv0nG2 = ELEMENT_F->NvnG2;
			Nv0nGc = ELEMENT_F->NvnGc;

			select_functions_cubature(&cubature,ELEMENT_F->type);
			select_functions_basis(&basis,ELEMENT_F->type);

			if (DB.Blending_HO) {
				PbMin = 1;
				PbMax = P;
				// Needs modification for PGc[P] != P (ToBeDeleted)
				if (P > 0 && PGc[P] != P)
					printf("Error: Add support. %d %d\n",P,PGc[P]), EXIT_MSG;
			} else {
				PbMin = P;
				PbMax = P;
			}

			for (Pb = PbMin; Pb <= PbMax; Pb++) {
				set_cubdata(CUBDATA,false,false,NodeTypeG[EclassF],dE-1,PGc[Pb],cubature); // free
				set_from_cubdata(CUBDATA,&Nv0nGc[Pb],NULL,&rst_v0Gc,NULL,NULL);

				ChiRefGs_vGs    = basis(PGs,    rst_v0Gs, Nv0nGs[1],&Nbf,dE-1);           // free
				ChiRefGc_vGc    = basis(PGc[Pb],rst_v0Gc, Nv0nGc[Pb],&Nbf,dE-1);          // free
				ChiRefGs_vProj  = basis(PGs,    rst_proj, NvnGc[P],&Nbf,dE-1);            // free
				ChiRefGc_vProj  = basis(PGc[Pb],rst_proj, NvnGc[P],&Nbf,dE-1);            // free

				IGs             = identity_d(Nv0nGs[1]);                             // free
				IGc             = identity_d(Nv0nGc[Pb]);                            // free
				ChiRefInvGs_vGs = inverse_d(Nv0nGs[1], Nv0nGs[1], ChiRefGs_vGs,IGs); // free
				ChiRefInvGc_vGc = inverse_d(Nv0nGc[Pb],Nv0nGc[Pb],ChiRefGc_vGc,IGc); // free

				if (Pb == P) {
					I_fGs_vGc[1][P][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],Nv0nGs[1],Nv0nGs[1],1.0,ChiRefGs_vProj,ChiRefInvGs_vGs); // keep
					if (P == 2) {
						set_cubdata(CUBDATA,false,false,NodeTypeG[EclassF],dE-1,P,cubature); // free
						set_from_cubdata(CUBDATA,&Nv0nG2[P],NULL,&rst_v0G2,NULL,NULL);

						ChiRefGs_vProj2 = basis(PGs,rst_proj2,NvnG2[P], &Nbf,dE-1); // free
						ChiRefG2_vG2    = basis(P,  rst_v0G2, Nv0nG2[P],&Nbf,dE-1); // free
						ChiRefG2_vProj2 = basis(P,  rst_proj2,NvnG2[P], &Nbf,dE-1); // free

						IG2             = identity_d(Nv0nG2[P]);                           // free
						ChiRefInvG2_vG2 = inverse_d(Nv0nG2[P],Nv0nG2[P],ChiRefG2_vG2,IG2); // free

						I_fGs_vG2[1][P][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnG2[P],Nv0nGs[1],Nv0nGs[1],1.0,ChiRefGs_vProj2,ChiRefInvGs_vGs); // keep
						I_fG2_vG2[Pb][P][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnG2[P],Nv0nG2[P],Nv0nG2[P],1.0,ChiRefG2_vProj2,ChiRefInvG2_vG2); // keep

						free(rst_v0G2);

						free(ChiRefGs_vProj2);
						free(ChiRefG2_vG2);
						free(ChiRefG2_vProj2);

						free(IG2);
						free(ChiRefInvG2_vG2);
					}
				}
				I_fGc_vGc[Pb][P][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],Nv0nGc[Pb],Nv0nGc[Pb],1.0,ChiRefGc_vProj,ChiRefInvGc_vGc); // keep

				free(rst_v0Gc);

				free(ChiRefGs_vGs);
				free(ChiRefGc_vGc);
				free(ChiRefGs_vProj);
				free(ChiRefGc_vProj);

				free(IGs);
				free(IGc);
				free(ChiRefInvGs_vGs);
				free(ChiRefInvGc_vGc);
			}

			free(rst_v0Gs);
			free(rst_proj);
			free(rst_proj2);
		}

		// EDGE related operators
		for (e = 0; dE == DMAX && e < Ne; e++) {
			Ve = e*NEREFMAX;

			ELEMENT_E = get_ELEMENT_type(LINE);
			EclassE   = get_Eclass(ELEMENT_E->type);

			rst_v0Gs = get_rst_vV(ELEMENT_E); // free
			rst_proj = malloc(NvnGc[P]*(dE-2) * sizeof *rst_proj); // free

			if (Blending == SZABO_BABUSKA) { // EclassE = C_TP (always)
				// These BCoords_E determine rst_proj according to Szabo(1991, p. 111) in 2D  (i.e. xi = L2-L1).
				BCoords_E = malloc(NvnGc[P]*Neve * sizeof *BCoords_E); // free

				for (n = 0; n < NvnGc[P]; n++) {
				for (ve = 0; ve < Neve; ve++) {
					BCoords_E[n*Neve+ve] = BCoords_V[n*Nve+VeEcon[e*NEVEMAX+ve]];
				}}

				// Compute projection of VOLUME nodes to EDGE
				mm_CTN_d(NvnGc[P],dE-2,Neve,BCoords_E,rst_v0Gs,rst_proj);
				free(BCoords_E);
			} else if (Blending == GORDON_HALL || EclassE == C_TP) {
				rst_vGs = get_rst_vV(ELEMENT); // free
				rst_vGsE = malloc(Nve * sizeof *rst_vGsE); // free

				if      (e/4 == 0) { dimE = 0; }
				else if (e/4 == 1) { dimE = 1; }
				else if (e/4 == 2) { dimE = 2; }

				for (n = 0; n < Nve; n++)
					rst_vGsE[n] = rst_vGs[n+dimE*Nve];
				free(rst_vGs);

				mm_CTN_d(NvnGc[P],dE-2,Nve,BCoords_V,rst_vGsE,rst_proj);
				free(rst_vGsE);
			} else {
				printf("Add Support.\n"), EXIT_MSG;
			}

			// Compute blending operators
			Nv0nGs = ELEMENT_E->NvnGs;
			Nv0nGc = ELEMENT_E->NvnGc;

			select_functions_cubature(&cubature,ELEMENT_E->type);
			select_functions_basis(&basis,ELEMENT_E->type);

			set_cubdata(CUBDATA,false,false,NodeTypeG[EclassE],dE-2,PGc[P],cubature); // free
			set_from_cubdata(CUBDATA,&Nv0nGc[P],NULL,&rst_v0Gc,NULL,NULL);

			ChiRefGs_vGs    = basis(PGs,   rst_v0Gs,Nv0nGs[1],&Nbf,dE-2);        // free
			ChiRefGc_vGc    = basis(PGc[P],rst_v0Gc,Nv0nGc[P],&Nbf,dE-2);        // free
			ChiRefGs_vProj  = basis(PGs,   rst_proj,NvnGc[P],&Nbf,dE-2);         // free
			ChiRefGc_vProj  = basis(PGc[P],rst_proj,NvnGc[P],&Nbf,dE-2);         // free

			IGs             = identity_d(Nv0nGs[1]);                           // free
			IGc             = identity_d(Nv0nGc[P]);                           // free
			ChiRefInvGs_vGs = inverse_d(Nv0nGs[1],Nv0nGs[1],ChiRefGs_vGs,IGs); // free
			ChiRefInvGc_vGc = inverse_d(Nv0nGc[P],Nv0nGc[P],ChiRefGc_vGc,IGc); // free

			I_eGs_vGc[1][P][Ve] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],Nv0nGs[1],Nv0nGs[1],1.0,ChiRefGs_vProj,ChiRefInvGs_vGs); // keep
			I_eGc_vGc[P][P][Ve] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],Nv0nGc[P],Nv0nGc[P],1.0,ChiRefGc_vProj,ChiRefInvGc_vGc); // keep


			free(ChiRefGs_vGs);
			free(ChiRefGc_vGc);
			free(ChiRefGs_vProj);
			free(ChiRefGc_vProj);

			free(IGs);
			free(IGc);
			free(ChiRefInvGs_vGs);
			free(ChiRefInvGc_vGc);

			free(rst_v0Gs);
			free(rst_v0Gc);
			free(rst_proj);
		}
	}
	free(CUBDATA);
}

static void setup_vertex_projection(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, **VeMask;
	double       ****I_vGs_vGs;

	// Initialize DB Parameters
	char         **NodeTypeG = DB.NodeTypeG;

	// Standard datatypes
	unsigned int vh, ve, veP2, Eclass, dE, Nve, NveP2, Nbf, Nvref, *Nvve;
	double       *rst_vGs, *TGs, *ChiInvGs_vGs, *ChiRefGs_vGs, *ChiGs_vGs;

	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef cubature;
	basis_tdef    basis;

	select_functions_cubature(&cubature,EType);
	select_functions_basis(&basis,EType);


	ELEMENT = get_ELEMENT_type(EType);

	Eclass = get_Eclass(EType);

	dE    = ELEMENT->d;
	Nve   = ELEMENT->Nve;
	NveP2 = ELEMENT->NveP2;
	Nvref = ELEMENT->Nvref;
	Nvve  = ELEMENT->Nvve;

	NvnGs        = ELEMENT->NvnGs;
	I_vGs_vGs    = ELEMENT->I_vGs_vGs;
	VeMask       = ELEMENT->VeMask[1][2];
	TGs          = ELEMENT->TGs[1][1][0];
	ChiInvGs_vGs = ELEMENT->ChiInvGs_vGs[1][1][0];


	// P2 geometry nodes
	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free
	set_cubdata(CUBDATA,false,false,NodeTypeG[Eclass],dE,2,cubature); // free
	set_from_cubdata(CUBDATA,&NvnGs[2],NULL,&rst_vGs,NULL,NULL);

	if (NvnGs[2] != NveP2)
		printf("Error: Incorrect.\n"), EXIT_MSG;

	ChiRefGs_vGs = basis(1,rst_vGs,NvnGs[2],&Nbf,dE); // free
	ChiGs_vGs    = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGs[2],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vGs,TGs); // free

	I_vGs_vGs[1][2][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGs[2],NvnGs[1],NvnGs[1],1.0,ChiGs_vGs,ChiInvGs_vGs); // keep

	free(rst_vGs);
	free(CUBDATA);

	free(ChiRefGs_vGs);
	free(ChiGs_vGs);

	for (vh = 0; vh < Nvref; vh++) {
		// Set Nvve[vh] for QUADs, HEXs, WEDGEs
		if (Nvve[vh] == 0)
			Nvve[vh] = Nve;

		VeMask[vh] = malloc(Nvve[vh] * sizeof *VeMask[vh]); // keep

		for (ve = 0; ve < Nvve[vh]; ve++) {
		for (veP2 = 0; veP2 < NveP2; veP2++) {
			if (array_norm_diff_d(Nve,&I_vGs_vGs[1][1][vh][ve*Nve],&I_vGs_vGs[1][2][0][veP2*Nve],"Inf") < EPS) {
				VeMask[vh][ve] = veP2;
				break;
			}
			if (veP2 == NveP2-1 && DB.TETrefineType != TET12) {
				array_print_d(Nvve[ve],Nve,I_vGs_vGs[1][1][vh],'R');
				array_print_d(NveP2,Nve,I_vGs_vGs[1][2][0],'R');
				printf("Error: Did not find VeMask.\n"), EXIT_MSG;
			}
		}}
	}
}

void setup_operators(void)
{
	/*
	 *	Comments:
	 *		QUADs/PYRs are set up if TETs are present and adaptation is enabled as TETs have the potential to refine
	 *		into PYRs when using TET6 refinement.
	 *		L2 projection operators must be set up before standard operators for TET and PYR ELEMENTs as several
	 *		preoperators must be subsequently overwritten.
	 */

	// Initialize DB Parameters
	unsigned int d     = DB.d,
	             Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int EType;

	// POINT
	EType = POINT;
	if (d == 1)
		setup_ELEMENT_FACE_ordering(EType);

	// LINE (Includes TP Class)
	EType = LINE;
	if (!DB.MPIrank && !DB.Testing)
		printf("    LINE\n");

	setup_ELEMENT_VeF(EType);
	setup_ELEMENT_plotting(EType);
	setup_ELEMENT_normals(EType);
	setup_ELEMENT_operators(EType);
	setup_L2_projection_operators(EType);
	if (d == 2)
		setup_ELEMENT_FACE_ordering(EType);

	// TRI
	EType = TRI;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    TRI\n");

		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
		setup_L2_projection_operators(EType);
		setup_blending(EType);
		setup_vertex_projection(EType);
		if (d == 3)
			setup_ELEMENT_FACE_ordering(EType);
	}

	// QUAD
	EType = QUAD;
	if (is_ELEMENT_present(EType) || (Adapt != ADAPT_0 && is_ELEMENT_present(TET))) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    QUAD\n");

		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
		setup_blending(EType);
		setup_vertex_projection(EType);
		if (d == 3)
			setup_ELEMENT_FACE_ordering(EType);
	}

	// HEX
	EType = HEX;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    HEX\n");

		setup_ELEMENT_VeE(EType);
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
		setup_blending(EType);
//		setup_vertex_projection(EType);
	}

	// TET/PYR
	if (Adapt != ADAPT_0 && is_ELEMENT_present(TET)) {
		setup_L2_projection_preoperators(TET);
		setup_L2_projection_preoperators(PYR);

		setup_L2_projection_operators(TET);
		setup_L2_projection_operators(PYR);

		memory_destructor_L2_projection(TET);
		memory_destructor_L2_projection(PYR);
	}

	// TET
	EType = TET;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    TET\n");

		setup_ELEMENT_operator_dependencies(PYR);

		setup_ELEMENT_VeE(EType);
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
		setup_blending(EType);
		setup_vertex_projection(EType);
	}

	// PYR
	EType = PYR;
	if (is_ELEMENT_present(EType) || (Adapt != ADAPT_0 && is_ELEMENT_present(TET) && DB.TETrefineType == TET6)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    PYR\n");

		setup_ELEMENT_VeE(EType);
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
//		setup_blending(EType);
		setup_vertex_projection(EType);
	}

	// WEDGE
	EType = WEDGE;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    WEDGE\n");

		setup_ELEMENT_VeE(EType);
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
//		setup_blending(EType);
//		setup_vertex_projection(EType);
	}

	// Free unused operators (including unused lower dimensional operators (if applicable) and sum factorized operators) (ToBeDeleted)
}
