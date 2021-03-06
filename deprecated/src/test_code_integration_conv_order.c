// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration_conv_order.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"

#include "adaptation.h"
#include "output_to_paraview.h"
#include "test_code_output_to_paraview.h"
#include "initialize_test_case.h"
#include "solver_Advection.h"
#include "solver_Poisson.h"
#include "solver_explicit.h"
#include "solver_implicit.h"
#include "compute_errors.h"
#include "element_functions.h"
#include "finalize_LHS.h"
#include "array_norm.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for (conv)ergence (order) integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

// Order for which to output the solution and mesh edges
#define OUTPUT_ORDER 5

static void h_adapt_test (void);

static void set_test_convorder_data(struct S_convorder *const data, char const *const TestName)
{
	/*
	 *	Comments:
	 *		As nodes having integration strength of order 2*P form a basis for the polynomial space of order P for TP
	 *		elements and SI elements (P <= 2), non-trivial L2 projection error requires:
	 *			TP         : InOrder_add > 1
	 *			SI (P <= 2): InOrder_add > 0
	 *
	 *		For the Peterson mesh suboptimal orders, the current meshes were generated with horizontal spacing such that
	 *		the converge is suboptimal for P = 1. For suboptimality of other orders, see Richter(2008) Figure 3 for the
	 *		required mesh spacing.
	 *
	 *	References:
	 *		Richter(2008)-On the Order of Convergence of the Discontinuous Galerkin Method for Hyperbolic Equations
	 */

	// default values
	data->PrintEnabled   = 0;
	data->Compute_L2proj = 0;
	data->SolveExplicit  = 1;
	data->SolveImplicit  = 1;
	data->AdaptiveRefine = 0;
	data->Adapt = ADAPT_HP;

	data->PMin  = 1;
	data->PMax  = 3;
	data->MLMin = 0;
	data->MLMax = 3;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;

	strcpy(data->argvNew[1],TestName);
	if (strstr(TestName,"Advection")) {
		data->IntOrder_add  = 2; // For certain cases, the exact solution is obtained if this is omitted.
		data->AdaptiveRefine = 0;
		data->MLMax  = 4;
		if (strstr(TestName,"Peterson_n-Cube")) {
			// Meshes are not nested and the spacing is selected for suboptimality for P = 1
			data->Adapt = ADAPT_P;
			data->PMin = 1;
			data->PMax = 1;
		}
	} else if (strstr(TestName,"Poisson")) {
		data->SolveExplicit = 0;
		data->AdaptiveRefine = 0;
		data->MLMax = 4;
		data->PG_add        = 0;
		data->IntOrder_add  = 2; // See comments
		if (strstr(TestName,"n-Cube")) {
			if (strstr(TestName,"LINE")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Cube_StraightLINE");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Ellipsoid_HollowSection")) {
			if (strstr(TestName,"ToBeCurvedTRI")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_ToBeCurvedTRI");
			} else if (strstr(TestName,"ToBeCurvedQUAD")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_ToBeCurvedQUAD");
			} else if (strstr(TestName,"TRI")) {
				if (strstr(TestName,"extended"))
					strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_CurvedTRI_extended");
				else
					strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_CurvedTRI");
			} else if (strstr(TestName,"QUAD")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_CurvedQUAD");
			} else if (strstr(TestName,"MIXED2D")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_CurvedMIXED2D");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			// ToBeModified (Include all relevant tests)
//	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_mixed");
//	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_Ellipsoidal_Section_2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_Ellipsoidal_Section_2D_QUAD");
//	strcpy(argvNew[1],"test/Test_Poisson_Ringleb2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_Ringleb2D_QUAD");
//	strcpy(argvNew[1],"test/Test_Poisson_HoldenRamp2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_GaussianBump2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_GaussianBump2D_mixed");
//	strcpy(argvNew[1],"test/Test_Poisson_3D_TET");
//	strcpy(argvNew[1],"test/Test_Poisson_3D_HEX");
//	strcpy(argvNew[1],"test/Test_Poisson_mixed3D_TP");
//	strcpy(argvNew[1],"test/Test_Poisson_mixed3D_HW");
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"Euler")) {
		if (strstr(TestName,"n-Cylinder_HollowSection")) {
			data->SolveExplicit = 0;
			if (strstr(TestName,"ToBeCurved")) {
				if (strstr(TestName,"MIXED2D")) {
//					data->PrintEnabled = 1;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
				} else if (strstr(TestName,"TET")) {
					// Blowing up in solver_implicit. Revisit once DPG is implemented. (ToBeDeleted)
					EXIT_UNSUPPORTED;

					data->PrintEnabled = 1;
					data->MLMax = 1;
					data->PMax  = 3;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedTET");
				} else if (strstr(TestName,"HEX")) {
					data->PrintEnabled = 1;
					data->MLMax = 2;
					data->PMax  = 2;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedHEX");
				} else if (strstr(TestName,"WEDGE")) {
					data->PrintEnabled = 1;
					data->MLMax = 2;
					data->PMax  = 2;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedWEDGE");
				} else if (strstr(TestName,"MIXED_TP")) {
					data->MLMax = 2;
					data->PMax  = 2;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_TP");
				} else if (strstr(TestName,"MIXED_HW")) {
					data->MLMax = 2;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_HW");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else if (strstr(TestName,"CurvedMIXED2D")) {
//				data->PrintEnabled = 1;
//test_print_warning("Modified Parameters"); PRINT_FILELINE;
//data->PMax = 2;
//data->MLMax = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Cube")) {
			data->SolveImplicit = 0;
			if (strstr(TestName,"Curved")) {
				EXIT_UNSUPPORTED;
			} else {
				data->PMin = 1;
				data->PMax = 3;
				data->MLMax = 5;
				if (strstr(TestName,"QUAD")) {
					if (strstr(TestName,"Stationary"))
						strcpy(data->argvNew[1],"test/Euler/Test_Euler_PeriodicVortex_Stationary_QUAD");
					else
						strcpy(data->argvNew[1],"test/Euler/Test_Euler_PeriodicVortex_QUAD");
				} else if (strstr(TestName,"TRI")) {
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_PeriodicVortex_TRI");
				} else {
					EXIT_UNSUPPORTED;
				}
			}
		} else if (strstr(TestName,"n-Elliptic_Pipe")) {
//			data->Compute_L2proj = 1;
			data->IntOrder_add  = 2;
			data->SolveExplicit = 0;
			data->PMin  = 2;
			data->PMax  = 5;
			data->MLMax = 4;
			if (strstr(TestName,"ToBeCurved")) {
				if (strstr(TestName,"TRI")) {
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_EllipticPipe_ToBeCurvedTRI");
				} else if (strstr(TestName,"QUAD")) {
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_EllipticPipe_ToBeCurvedQUAD");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
			EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Parabolic_Pipe")) {
//			data->Compute_L2proj = 1;
			data->IntOrder_add  = 2;
			data->SolveExplicit = 0;
			data->PMin  = 2;
			data->PMax  = 4;
			data->MLMax = 4;
			if (strstr(TestName,"ToBeCurved")) {
				if (strstr(TestName,"TRI")) {
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_ParabolicPipe_ToBeCurvedTRI");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Sinusoidal_Pipe")) {
//			data->SolveExplicit = 0;
			if (strstr(TestName,"ToBeCurved")) {
				if (strstr(TestName,"TRI")) {
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_SinusoidalPipe_ToBeCurvedTRI");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-GaussianBump")) {
			if (strstr(TestName,"BezierToBeCurved")) {
				if (strstr(TestName, "QUAD")){
data->PrintEnabled = 1;
data->PMin  = 3;
data->Adapt = ADAPT_0;
					// The ToBeCurved QUAD mesh for the gaussian bump case
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_GaussianBump_BezierToBeCurvedQUAD");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else if (strstr(TestName,"BezierCurvedQUAD")) {
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_GaussianBump_BezierCurvedQUAD");
			} else if (strstr(TestName,"Curved")) {
				if (strstr(TestName,"CurvedQUAD")) {
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_GaussianBump_CurvedQUAD");
				}else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"NavierStokes")) {
		data->PG_add = 0;
//		data->PrintEnabled = 1;
//test_print_warning("Modified Parameters"); PRINT_FILELINE;
//data->MLMax = 0;
//data->PMax  = 1;

		strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_");
		if (strstr(TestName,"n-Cylinder_Hollow")) {
			strcat(data->argvNew[1],"TaylorCouette_");
			if (strstr(TestName,"ToBeCurved")) {
				strcat(data->argvNew[1],"ToBeCurved");
				if (strstr(TestName,"TRI")) {
					strcat(data->argvNew[1],"TRI");
				} else if (strstr(TestName,"QUAD")) {
					strcat(data->argvNew[1],"QUAD");
				} else if (strstr(TestName,"MIXED2D")) {
					data->MLMax = 4;
					strcat(data->argvNew[1],"MIXED2D");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Cube")) {
			strcat(data->argvNew[1],"PlaneCouette_");
			if (strstr(TestName,"Straight")) {
				if (strstr(TestName,"QUAD")) { strcat(data->argvNew[1],"QUAD"); }
				else { EXIT_UNSUPPORTED; }
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		printf("%s\n",TestName); EXIT_UNSUPPORTED;
	}
}

void test_conv_order(struct S_convorder *const data, char const *const TestName)
{
	/*
	 *	Expected Output:
	 *		Optimal convergence orders in L2 for the solution (P+1) and its gradients (P).
	 */

	set_test_convorder_data(data,TestName);

	int  const               nargc   = data->nargc;
	char const *const *const argvNew = (char const *const *const) data->argvNew;

	bool const PrintEnabled   = data->PrintEnabled,
	           Compute_L2proj = data->Compute_L2proj,
	           SolveExplicit  = data->SolveExplicit,
	           SolveImplicit  = data->SolveImplicit,
	           AdaptiveRefine = data->AdaptiveRefine;

	unsigned int const Adapt = data->Adapt,
	                   PMin  = data->PMin,
	                   PMax  = data->PMax,
	                   MLMin = data->MLMin,
	                   MLMax = data->MLMax;

	TestDB.PGlobal       = 1; // ToBeModified (Not important for Adapt != ADAPT_0)
//	TestDB.PGlobal       = 4; // ToBeModified (Not important for Adapt != ADAPT_0)
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	double *mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

	if (!SolveExplicit && !SolveImplicit)
		EXIT_UNSUPPORTED; // Enabled at least one

	unsigned int pass = 0;
	if (Adapt != ADAPT_0 && Adapt != ADAPT_P) {
		TestDB.ML = DB.ML;
		code_startup(nargc,argvNew,0,2);
	}

	for (size_t P = PMin; P <= PMax; P++) {
	for (size_t ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		if (Adapt == ADAPT_0) {
			code_startup(nargc,argvNew,0,1);
		} else if (Adapt == ADAPT_P) {
			code_startup(nargc,argvNew,0,1);
			// mesh_to_order(P);
		} else if (Adapt == ADAPT_H || Adapt == ADAPT_HP) {
			mesh_to_level(ML);
			if (ML == MLMin && AdaptiveRefine)
				h_adapt_test();
			mesh_to_order(P);
		}

		if (!SolveImplicit)
			initialize_test_case(0);

		// Output mesh edges to paraview
		if (DB.d > 1 && TestDB.PGlobal == OUTPUT_ORDER && TestDB.ML <= 2) {
			char *const fNameOut = get_fNameOut("MeshEdges_");
			output_to_paraview(fNameOut);
			free(fNameOut);
		}

		if (Compute_L2proj) {
			if (!(strstr(TestName,"Poisson") || strstr(TestName,"Euler")))
				EXIT_UNSUPPORTED;

			DB.init_with_L2 = true;
			initialize_test_case(0);
			DB.init_with_L2 = false;
		} else {
			if (strstr(TestName,"Advection")) {
				solver_Advection(PrintEnabled);
			} else if (strstr(TestName,"Poisson")) {
				solver_Poisson(PrintEnabled);
			} else if (strstr(TestName,"Euler") || strstr(TestName,"NavierStokes")) {
//				if (SolveExplicit)
//				if (!SolveImplicit || (ML == MLMin && P == PMin))
				if (SolveExplicit && (!SolveImplicit || (ML <= MLMin+1))) {
					bool const update_SolverType = (strstr(DB.SolverType,"Implicit") ? 1 : 0);
					if (update_SolverType)
						strcpy(DB.SolverType,"Explicit");
					solver_explicit(PrintEnabled);
					if (update_SolverType)
						strcpy(DB.SolverType,"Implicit");
				}

				if (SolveImplicit)
					solver_implicit(PrintEnabled);
			} else {
				EXIT_UNSUPPORTED;
			}
		}

		// Output to paraview
		if (TestDB.ML <= 1 || (TestDB.PGlobal == 1) || (TestDB.PGlobal == OUTPUT_ORDER && TestDB.ML <= 4)) {
			char *const fNameOut = get_fNameOut("SolFinal_");
			output_to_paraview(fNameOut);
			free(fNameOut);
		}

		compute_errors_global();

		if (PrintEnabled) {
			compute_dof();
			printf("ML, P, dof: %zu %zu %d\n",ML,P,DB.dof);
		}

// ERROR FOR BEZIER WHEN WORKING WITH MESH REGULARITY
		if (P == PMin)
			evaluate_mesh_regularity(&mesh_quality[ML-MLMin]);

		if (P == PMax && ML == MLMax) {
			check_convergence_orders(MLMin,MLMax,PMin,PMax,&pass,PrintEnabled);
			check_mesh_regularity(mesh_quality,MLMax-MLMin+1,&pass,PrintEnabled);
		}

		if (Adapt == ADAPT_0 || Adapt == ADAPT_P) {
			set_PrintName("conv_orders",data->PrintName,&data->omit_root);
			code_cleanup();
		}


	}}
	if (Adapt == ADAPT_H || Adapt == ADAPT_HP) {
		set_PrintName("conv_orders",data->PrintName,&data->omit_root);
		code_cleanup();
	}
	free(mesh_quality);

	test_print2(pass,data->PrintName);
}

void h_adapt_test(void)
{
	/*
	 * Purpose:
	 *		Perform local mesh refinement around specified XYZref coordinate locations.
	 */

	// Initialize DB Parameters
	char         *Geometry = DB.Geometry;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int NrefMax = 5, MLMax = 5;

	unsigned int Nref, NML[NrefMax], CurvedOnly[NrefMax];
	double       *XYZref;

	XYZref = malloc(NrefMax*DMAX * sizeof *XYZref); // free

	if (TestDB.ML > 0)
		printf("Error: Only enter for ML == 0.\n"), EXIT_MSG;

	if (strstr(Geometry,"n-Cylinder") || strstr(Geometry,"n-Cube") ||
	    strstr(Geometry,"GaussianBump")) {
		Nref = 0;
	} else if (strstr(Geometry,"JoukowskiSymmetric")) {
		double a  = DB.JSa,
		       xL = DB.JSxL;

		Nref = 2;

		NML[0] = 0;
		NML[1] = 3;

		XYZref[0+0*DMAX] = xL;  XYZref[1+0*DMAX] = 0.0; XYZref[2+0*DMAX] = 0.0;
		XYZref[0+1*DMAX] = 2*a; XYZref[1+1*DMAX] = 0.0; XYZref[2+1*DMAX] = 0.0;
	} else if (strstr(Geometry,"n-Ellipsoid")) {
		if (strstr(DB.PDE,"Poisson")) {
			Nref = 1;
			MLMax = 2;

			unsigned int i = 0;
			NML[i] = 2; CurvedOnly[i] = 0; i++;

			i = 0;
			if (strstr(DB.MeshType,"ToBeCurved")) {
				// May need to add a an override for the refinement of one of the volumes to get the desired
				// non-conforming mesh.
				XYZref[0+i*DMAX] =  0.0;     XYZref[1+i*DMAX] = DB.aOut; XYZref[2+i*DMAX] = 0.0; i++;
			} else if (strstr(DB.MeshType,"Curved")) {
				XYZref[0+i*DMAX] =  0.0;     XYZref[1+i*DMAX] = DB.bOut; XYZref[2+i*DMAX] = 0.0; i++;
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(DB.PDE,"Euler")) {
			Nref = 5;

			unsigned int i = 0;
			NML[i] = 2; CurvedOnly[i] = 0; i++;
			NML[i] = 1; CurvedOnly[i] = 0; i++;
			NML[i] = 2; CurvedOnly[i] = 0; i++;
			NML[i] = 1; CurvedOnly[i] = 0; i++;
			NML[i] = 1; CurvedOnly[i] = 0; i++;

			i = 0;
			XYZref[0+i*DMAX] =  DB.aIn;  XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
			XYZref[0+i*DMAX] =  DB.aOut; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
			XYZref[0+i*DMAX] = -DB.aIn;  XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
			XYZref[0+i*DMAX] = -DB.aOut; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
			XYZref[0+i*DMAX] =  0.0;     XYZref[1+i*DMAX] = DB.bIn; XYZref[2+i*DMAX] = 0.0; i++;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(Geometry,"EllipsoidalBump")) {
		Nref = 2;

		unsigned int i = 0;
		NML[i] = 2; CurvedOnly[i] = 1; i++;
		NML[i] = 2; CurvedOnly[i] = 1; i++;

		i = 0;
		XYZref[0+i*DMAX] =  DB.aIn; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] = -DB.aIn; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	// Store indices of VOLUMEs to be updated and the adaptation type (HREFINE here) in hp_update
	unsigned int ML;
	for (ML = 0; ML < MLMax; ML++) {
		unsigned int NVglobal = DB.NVglobal;

		struct S_VOLUME *VOLUME;

		unsigned int *hp_update;
		hp_update = calloc(NVglobal , sizeof *hp_update); // free

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int n, ve, dim, Nve, indexg;
			double       *XYZ_vV;

			struct S_ELEMENT *ELEMENT;

			ELEMENT = get_ELEMENT_type(VOLUME->type);

			Nve = ELEMENT->Nve;

			indexg = VOLUME->indexg;

			// Store XYZ_vV as a row vector
			XYZ_vV = malloc(Nve*d * sizeof *XYZ_vV); // free
			for (ve = 0; ve < Nve; ve++) {
			for (dim = 0; dim < d; dim++) {
				XYZ_vV[ve*d+dim] = VOLUME->XYZ_vV[ve+dim*Nve];
			}}

			for (n = 0; n < Nref; n++) {
				if (VOLUME->level < NML[n]) {
					if (CurvedOnly[n] && !VOLUME->curved)
						continue;
					// Check if one of the XYZ_vV matches any of the specified XYZref
					for (ve = 0; ve < Nve; ve++) {
						if (array_norm_diff_d(d,&XYZref[n*DMAX],&XYZ_vV[ve*d],"Inf") < EPS) {
							hp_update[indexg] = HREFINE;
							break;
						}
					}
				}
			}
			free(XYZ_vV);
		}

		// Add additional indices to the list of hp_update to ensure that the resulting mesh will not be more than
		// 1-irregular.
		ensure_1irregular(hp_update);

		// Mark VOLUMEs for refinement and perform the mesh update
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (hp_update[VOLUME->indexg] == 0) {
				// Do nothing
			} else if (hp_update[VOLUME->indexg] == HREFINE) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = HREFINE;
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}
		}
		free(hp_update);
		mesh_update();
	}
	free(XYZref);
}
