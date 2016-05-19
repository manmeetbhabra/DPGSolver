// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of cubature_PYR.
 *
 *	Comments:
 *		Other than the grouping based on 1/4 symmetries, the ordering of the pyramidal nodes is arbitrary.
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_cubature_PYR(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		P, NodeType.
	 *
	 *	Expected output:
	 *
	 *		P = 2, GL:
	 *			rst = [ See below ]
	 *			w   = N/A
	 *			symms = [ 4 4 4 1 1 ]
	 *			Ns  = 5
	 *			Nn  = 14
	 *
	 *		P = 3, GLL:
	 *			rst = [ See below ]
	 *			w   = N/A
	 *			symms = [ 8 4 4 4 4 1 1 ]
	 *			Ns  = 8
	 *			Nn  = 30
	 *
	 *		P = 4, WV:
	 *			rst = [ See below ]
	 *			w   = [ See below ]
	 *			symms = [ 4 4 1 1 ]
	 *			Ns  = 4
	 *			Nn  = 10
	 */

	unsigned int Nn, Ns, P, d;
	unsigned int *symms;
	double *rst, *w;

	d = 3;

	// GL (P = 2)
	P = 2;
	unsigned int Nn2_GL = 14, Ns2_GL = 5;
	unsigned int symms2_GL[5] = { 4, 4, 4, 1, 1 };
	double rst2_GL[42] = { -0.687298334620742, -0.687298334620742, -0.123458488793238,
	                        0.687298334620742, -0.687298334620742, -0.123458488793238,
	                        0.687298334620742,  0.687298334620742, -0.123458488793238,
	                       -0.687298334620742,  0.687298334620742, -0.123458488793238,
	                        0.000000000000000, -0.687298334620742, -0.123458488793238,
	                        0.687298334620742, -0.000000000000000, -0.123458488793238,
	                        0.000000000000000,  0.687298334620742, -0.123458488793238,
	                       -0.687298334620742,  0.000000000000000, -0.123458488793238,
	                       -0.288675134594813, -0.288675134594813,  0.424264068711929,
	                        0.288675134594813, -0.288675134594813,  0.424264068711929,
	                        0.288675134594813,  0.288675134594813,  0.424264068711929,
	                       -0.288675134594813,  0.288675134594813,  0.424264068711929,
	                        0.000000000000000,  0.000000000000000, -0.123458488793238,
	                        0.000000000000000,  0.000000000000000,  0.971986626217095};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn2_GL,d,1.0,rst2_GL,d,Nn2_GL);

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,P,d,"GL");

	pass = 0;
	if (array_norm_diff_d(Nn2_GL*d,rst2_GL,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns2_GL,symms2_GL,symms,"Inf") < EPS &&
	    Nn2_GL == Nn && Ns2_GL == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_PYR (P2, GL) :                          ");
	test_print(pass);

	free(rst);
	free(symms);

	// GLL (P = 3)
	P = 3;
	unsigned int Nn3_GLL = 30, Ns3_GLL = 9;
	unsigned int symms3_GLL[9] = { 4, 4, 4, 4, 4, 4, 4, 1, 1 };
	double rst3_GLL[90] = { -0.447213595499958, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000, -0.447213595499958, -0.282842712474619,
	                         0.447213595499958,  1.000000000000000, -0.282842712474619,
	                        -1.000000000000000,  0.447213595499958, -0.282842712474619,
	                         0.447213595499958, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000,  0.447213595499958, -0.282842712474619,
	                        -0.447213595499958,  1.000000000000000, -0.282842712474619,
	                        -1.000000000000000, -0.447213595499958, -0.282842712474619,
	                        -1.000000000000000, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000,  1.000000000000000, -0.282842712474619,
	                        -1.000000000000000,  1.000000000000000, -0.282842712474619,
	                        -0.723606797749979, -0.723606797749979,  0.108036302695091,
	                         0.723606797749979, -0.723606797749979,  0.108036302695091,
	                         0.723606797749979,  0.723606797749979,  0.108036302695091,
	                        -0.723606797749979,  0.723606797749979,  0.108036302695091,
	                        -0.447213595499958, -0.447213595499958, -0.282842712474619,
	                         0.447213595499958, -0.447213595499958, -0.282842712474619,
	                         0.447213595499958,  0.447213595499958, -0.282842712474619,
	                        -0.447213595499958,  0.447213595499958, -0.282842712474619,
	                         0.000000000000000, -0.723606797749979,  0.108036302695091,
	                         0.723606797749979,  0.000000000000000,  0.108036302695091,
	                         0.000000000000000,  0.723606797749979,  0.108036302695091,
	                        -0.723606797749979,  0.000000000000000,  0.108036302695091,
	                        -0.276393202250021, -0.276393202250021,  0.740491834728766,
	                         0.276393202250021, -0.276393202250021,  0.740491834728766,
	                         0.276393202250021,  0.276393202250021,  0.740491834728766,
	                        -0.276393202250021,  0.276393202250021,  0.740491834728766,
	                        -0.000000000000000, -0.000000000000000,  0.108036302695091,
	                         0.000000000000000,  0.000000000000000,  1.131370849898476};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn3_GLL,d,1.0,rst3_GLL,d,Nn3_GLL);

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,P,d,"GLL");

	pass = 0;
	if (array_norm_diff_d(Nn3_GLL*d,rst3_GLL,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns3_GLL,symms3_GLL,symms,"Inf") < EPS &&
	    Nn3_GLL == Nn && Ns3_GLL == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (P3, GLL):                          ");
	test_print(pass);

	free(rst);
	free(symms);

	// WV (P = 4)
	P = 4;
	unsigned int Nn4_WV = 10, Ns4_WV = 4;
	unsigned int symms4_WV[4] = { 4, 4, 1, 1 };
	double rst4_WV[30] = { -0.657966997121690, -0.657966997121690, -0.227337257085045,
	                        0.657966997121690, -0.657966997121690, -0.227337257085045,
	                        0.657966997121690,  0.657966997121690, -0.227337257085045,
	                       -0.657966997121690,  0.657966997121690, -0.227337257085045,
	                       -0.000000000000000, -0.650581556398233,  0.173077324153007,
	                        0.650581556398233,  0.000000000000000,  0.173077324153007,
	                        0.000000000000000,  0.650581556398233,  0.173077324153007,
	                       -0.650581556398233,  0.000000000000000,  0.173077324153007,
	                        0.000000000000000,  0.000000000000000, -0.105872336234184,
	                        0.000000000000000,  0.000000000000000,  0.674909082451912};

	double w4_WV[10] = { 0.119772745994674,
	                     0.119772745994674,
	                     0.119772745994674,
	                     0.119772745994674,
	                     0.200487565609086,
	                     0.200487565609086,
	                     0.200487565609086,
	                     0.200487565609086,
	                     0.390103085029384,
	                     0.214473751719704};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn4_WV,d,1.0,rst4_WV,d,Nn4_WV);

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WV");

	pass = 0;
	if (array_norm_diff_d(Nn4_WV*d,rst4_WV,rst,"Inf")    < EPS*10 &&
	    array_norm_diff_d(Nn4_WV,w4_WV,w,"Inf")          < EPS    &&
	    array_norm_diff_ui(Ns4_WV,symms4_WV,symms,"Inf") < EPS    &&
	    Nn4_WV == Nn && Ns4_WV == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (P4, WV) :                          ");
	test_print(pass);

	free(rst);
	free(w);
	free(symms);
}
