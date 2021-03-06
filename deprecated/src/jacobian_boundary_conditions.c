// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "jacobian_boundary_conditions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "variable_functions.h"
#include "boundary_conditions.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Compute the jacobian of boundary conditions for given input parameters: XYZ, WL, WOut, nL
 *
 *	Comments:
 *		The jacobian_boundary_* functions compute the transpose of the standard dWdW matrices. (i.e. the terms are
 *		ordered in memory as: dWB0/dWL0 dWB1/dWL0 ... dWB(Neq)/dWL0 dWB0/dWL1 ...)
 *		If this is found to be slow after profiling, separate cases for different dimensions to avoid unnecessary
 *		operations (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

static void jacobian_boundary_Advection         (struct S_BC *const BCdata);
static void jacobian_boundary_Poisson           (struct S_BC *const BCdata);
static void jacobian_boundary_Riemann           (struct S_BC *const BCdata);
static void jacobian_boundary_SlipWall          (struct S_BC *const BCdata);
static void jacobian_boundary_BackPressure      (struct S_BC *const BCdata);
static void jacobian_boundary_Total_TP          (struct S_BC *const BCdata);
static void jacobian_boundary_SupersonicInflow  (struct S_BC *const BCdata);
static void jacobian_boundary_SupersonicOutflow (struct S_BC *const BCdata);
static void jacobian_boundary_NoSlip_Dirichlet  (struct S_BC *const BCdata);
static void jacobian_boundary_NoSlip_Adiabatic  (struct S_BC *const BCdata);

void compute_jacobian_boundary_values(struct S_BC *const BCdata)
{
	switch((BCdata->BC) % BC_STEP_SC) {
		case BC_RIEMANN:          jacobian_boundary_Riemann(BCdata);           break;
		case BC_SLIPWALL: {
			if (EXACT_SLIPWALL) {
				unsigned int const d    = BCdata->d,
				                   Nvar = d+2,
				                   Neq  = d+2,
				                   Nn   = BCdata->Nn,
				                   Nel  = BCdata->Nel;

				unsigned int const NnTotal = Nn*Nel;

				double *const dWBdWL = BCdata->dWBdWL;
				for (size_t i = 0; i < NnTotal*Nvar*Neq; i++)
					dWBdWL[i] = 0.0;
			} else if (EXACT_NORMAL) {
				double const *const nL = BCdata->nL;

				BCdata->nL = compute_exact_boundary_normal(BCdata);
				jacobian_boundary_SlipWall(BCdata);

				free((double *) BCdata->nL);
				BCdata->nL = nL;
			} else {
				jacobian_boundary_SlipWall(BCdata);
			}
			break;
		}
		case BC_BACKPRESSURE:     jacobian_boundary_BackPressure(BCdata);      break;
		case BC_TOTAL_TP:         jacobian_boundary_Total_TP(BCdata);          break;
		case BC_SUPERSONIC_IN:    jacobian_boundary_SupersonicInflow(BCdata);  break;
		case BC_SUPERSONIC_OUT:   jacobian_boundary_SupersonicOutflow(BCdata); break;
		case BC_NOSLIP_T:         jacobian_boundary_NoSlip_Dirichlet(BCdata);  break;
		case BC_NOSLIP_ADIABATIC: jacobian_boundary_NoSlip_Adiabatic(BCdata);  break;
		case BC_DIRICHLET:        // fallthrough
		case BC_NEUMANN:          jacobian_boundary_Poisson(BCdata);           break;
		case BC_INFLOW:           // fallthrough
		case BC_OUTFLOW:          jacobian_boundary_Advection(BCdata);         break;
		default:                  EXIT_UNSUPPORTED;                            break;
	}
}

static void jacobian_boundary_Advection(struct S_BC *const BCdata)
{
	unsigned int const BC_index = (BCdata->BC) % BC_STEP_SC;
	if (!(BC_index == BC_INFLOW || BC_index == BC_OUTFLOW))
		EXIT_UNSUPPORTED;

	/*
	 *	Comments:
	 *		It is implicitly assumed that Nvar = Neq = 1.
	 */

	unsigned int const Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double *const dWBdWL = BCdata->dWBdWL;

	if (BC_index == BC_INFLOW) {
		for (size_t n = 0; n < NnTotal; n++)
			dWBdWL[n] = 0.0;
	} else if (BC_index == BC_OUTFLOW) {
		for (size_t n = 0; n < NnTotal; n++)
			dWBdWL[n] = 1.0;
	}
}

static void jacobian_boundary_Poisson(struct S_BC *const BCdata)
{
	unsigned int const BC_index = (BCdata->BC) % BC_STEP_SC;
	if (!(BC_index == BC_DIRICHLET || BC_index == BC_NEUMANN))
		EXIT_UNSUPPORTED;

	/*
	 *	Comments:
	 *		The assumption that Nvar = 1 was used below.
	 *		Currently does not support computing dQBd*.
	 */

	unsigned int const Nn  = BCdata->Nn,
	                   Nel = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double *const dWBdWL = BCdata->dWBdWL;

	if (BC_index == BC_DIRICHLET) {
		for (size_t n = 0; n < NnTotal; n++)
			dWBdWL[n] = -1.0;
	} else if (BC_index == BC_NEUMANN) {
		for (size_t n = 0; n < NnTotal; n++)
			dWBdWL[n] =  1.0;
	}
}

static void jacobian_boundary_Riemann(struct S_BC *const BCdata)
{
	/*
	 *	Jacobian Matrices [var * eq]
	 *
	 *	Supersonic Inlet
	 *	dWBdWL = [  0  0  0  0  0
	 *	            0  0  0  0  0
	 *	            0  0  0  0  0
	 *	            0  0  0  0  0
	 *	            0  0  0  0  0 ]
	 *
	 *	Supersonic Outlet
	 *	dWBdWL = [  1  0  0  0  0
	 *	            0  1  0  0  0
	 *	            0  0  1  0  0
	 *	            0  0  0  1  0
	 *	            0  0  0  0  1 ]
	 *
	 *	Subsonic Inlet/Outlet
	 *		See below.
	 *
	 */

	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL  = BCdata->nL,
	             *const XYZ = BCdata->XYZ,
	             *const WL  = BCdata->WL;

	double *const dWdW = BCdata->dWBdWL;

	// Standard datatypes
	unsigned int i, iMax, n, eq, var, InddWdW;
	double       rhoL, rhoL_inv, uL, vL, wL, EL, V2L, pL, rhoR, uR, vR, wR, pR,
	             cL, RL, VnL, cR, RR, VnR, c, Vn,
	             X, Y, n1, n2, n3, *dWdW_ptr[Neq*Neq];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *n_ptr, *X_ptr, *Y_ptr;

	// silence
	rhovL_ptr = rhowL_ptr = NULL;
	n2 = n3 = 0.0;

	double zeros[NnTotal];

	for (i = 0; i < NnTotal; i++)
		zeros[i] = 0.0;

	rhoL_ptr  = &WL[        0];
	rhouL_ptr = &WL[NnTotal*1];
	EL_ptr    = &WL[NnTotal*(d+1)];

	n_ptr = nL;

	if (d == 3) {
		rhovL_ptr = &WL[NnTotal*2];
		rhowL_ptr = &WL[NnTotal*3];
	} else if (d == 2) {
		rhovL_ptr = &WL[NnTotal*2];
		rhowL_ptr = zeros;
	} else if (d == 1) {
		rhovL_ptr = zeros;
		rhowL_ptr = zeros;
	}

	X_ptr = &XYZ[NnTotal*0];
	Y_ptr = &XYZ[NnTotal*1];

/// Alternate order in updated code (Fastest on var, then eq (where eq is actually var_out)).
	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
		dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
	}}

	for (n = 0; n < NnTotal; n++) {
		InddWdW = 0;

		// Inner VOLUME
		rhoL     = *rhoL_ptr++;
		rhoL_inv = 1.0/rhoL;

		uL   = (*rhouL_ptr++)*rhoL_inv;
		vL   = (*rhovL_ptr++)*rhoL_inv;
		wL   = (*rhowL_ptr++)*rhoL_inv;
		EL   = *EL_ptr++;

		V2L = uL*uL+vL*vL+wL*wL;

		pL  = GM1*(EL-0.5*rhoL*V2L);

		n1 = *n_ptr++;
		if      (d == 3) { n2 = *n_ptr++; n3 = *n_ptr++; }
		else if (d == 2) { n2 = *n_ptr++; n3 = 0.0;      }
		else if (d == 1) { n2 = 0.0;      n3 = 0.0;      }

		VnL = uL*n1+vL*n2+wL*n3;

		X = *X_ptr++;
		Y = *Y_ptr++;

		// Outer VOLUME
		get_boundary_values(X,Y,&rhoR,&uR,&vR,&wR,&pR);
		VnR = n1*uR+n2*vR+n3*wR;

		cL = sqrt(GAMMA*pL/rhoL);
		cR = sqrt(GAMMA*pR/rhoR);

		// Riemann invariants
		RL = VnL + 2.0/GM1*cL;
		RR = VnR - 2.0/GM1*cR;

		Vn = 0.5*(RL+RR);
		c  = 0.25*GM1*(RL-RR);

		if (fabs(Vn) >= c) { // Supersonic
			if (Vn < 0.0) { // Inlet
//printf("j: Sup Inlet\n");
				for (var = 0; var < Nvar; var++) {
				for (eq = 0; eq < Neq; eq++) {
					*dWdW_ptr[InddWdW++] = 0.0;
				}}
			} else { // Outlet
//printf("j: Sup Outlet\n");
				for (var = 0; var < Nvar; var++) {
				for (eq = 0; eq < Neq; eq++) {
					if (var != eq)
						*dWdW_ptr[InddWdW++] = 0.0;
					else
						*dWdW_ptr[InddWdW++] = 1.0;
				}}
			}
		} else { // Subsonic
			double dcLdW, rho, u, v, w, V2, ut, vt, wt, un, vn, wn, cnst1, drhodW, dudW, dvdW, dwdW, dpdW,
			       drhoLdW[Nvar], duLdW[Nvar], dvLdW[Nvar], dwLdW[Nvar], dpLdW[Nvar],
			       dVnLdW[Nvar], dRLdW[Nvar], dcdW[Nvar];

			// silence
			un = vn = wn = 0.0;

			if (d == 3) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = -wL; dpLdW[4]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
				dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;

				for (var = 0; var < Nvar; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1 + dvLdW[var]*n2 + dwLdW[var]*n3;
				}

				un = Vn*n1;
				vn = Vn*n2;
				wn = Vn*n3;
			} else if (d == 2) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;

				for (var = 0; var < Nvar; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1 + dvLdW[var]*n2;
				}

				un = Vn*n1;
				vn = Vn*n2;
			} else if (d == 1) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;

				for (var = 0; var < Nvar; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1;
				}

				un = Vn*n1;
			}

			for (var = 0; var < Nvar; var++) {
				dcLdW      = 0.5*GAMMA/(cL*rhoL*rhoL)*(dpLdW[var]*rhoL-pL*drhoLdW[var]);
				dRLdW[var] = dVnLdW[var] + 2.0/GM1*dcLdW;
				dcdW[var]  = 0.25*GM1*dRLdW[var];
			}

			if (Vn < 0.0) { // Inlet
//printf("j: Sub Inlet\n");
				double sR;

				sR  = sqrt(pR/pow(rhoR,GAMMA));
				if (d == 3) {
					for (var = 0; var < Nvar; var++) {
						drhodW = pow(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						vt = vR - VnR*n2;
						wt = wR - VnR*n3;

						u = un + ut;
						v = vn + vt;
						w = wn + wt;
						V2 = u*u+v*v+w*w;

						dudW = 0.5*dRLdW[var]*n1;
						dvdW = 0.5*dRLdW[var]*n2;
						dwdW = 0.5*dRLdW[var]*n3;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = drhodW*w + rho*dwdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
					}
				} else if (d == 2) {
					for (var = 0; var < Nvar; var++) {
						drhodW = pow(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						vt = vR - VnR*n2;

						u = un + ut;
						v = vn + vt;
						V2 = u*u+v*v;

						dudW = 0.5*dRLdW[var]*n1;
						dvdW = 0.5*dRLdW[var]*n2;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW));
					}
				} else if (d == 1) {
					for (var = 0; var < Nvar; var++) {
						drhodW = pow(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						u  = un + ut;
						V2 = u*u;

						dudW = 0.5*dRLdW[var]*n1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW));
					}
				}
			} else { // Outlet
//printf("j: Sub Outlet\n");
				double sL, dsLdW;

				sL  = sqrt(pL/pow(rhoL,GAMMA));
				if (d == 3) {
					for (var = 0; var < Nvar; var++) {
						dsLdW = 0.5*sqrt(pow(rhoL,GAMMA)/pL)/pow(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow(rhoL,GAMMA)-GAMMA*pL*pow(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow(GAMMA,-1.0/GM1)*2.0/GM1*pow(c,-GM3/GM1)*
								 pow(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						vt = vL - VnL*n2;
						wt = wL - VnL*n3;

						u = un + ut;
						v = vn + vt;
						w = wn + wt;
						V2 = u*u+v*v+w*w;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dvdW = dvLdW[var]+n2*cnst1;
						dwdW = dwLdW[var]+n3*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = drhodW*w + rho*dwdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
					}
				} else if (d == 2) {
					for (var = 0; var < Nvar; var++) {
						dsLdW = 0.5*sqrt(pow(rhoL,GAMMA)/pL)/pow(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow(rhoL,GAMMA)-GAMMA*pL*pow(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow(GAMMA,-1.0/GM1)*2.0/GM1*pow(c,-GM3/GM1)*
								 pow(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						vt = vL - VnL*n2;

						u = un + ut;
						v = vn + vt;
						V2 = u*u+v*v;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dvdW = dvLdW[var]+n2*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW));
					}
				} else if (d == 1) {
					for (var = 0; var < Nvar; var++) {
						dsLdW = 0.5*sqrt(pow(rhoL,GAMMA)/pL)/pow(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow(rhoL,GAMMA)-GAMMA*pL*pow(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow(GAMMA,-1.0/GM1)*2.0/GM1*pow(c,-GM3/GM1)*
								 pow(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						u = un + ut;
						V2 = u*u;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW));
					}
				}
			}
		}
		for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWdW_ptr[i]++;
	}
}

static void jacobian_boundary_SlipWall(struct S_BC *const BCdata)
{
	/*
	 *	Jacobian Matrix [var * eq]
	 *
	 *	dWBdWL = [  1  0          0         0         0
	 *	            0  1-2*n1*n1   -2*n2*n1  -2*n3*n1 0
	 *	            0   -2*n1*n2  1-2*n2*n2  -2*n3*n2 0
	 *	            0   -2*n1*n3   -2*n2*n3 1-2*n3*n3 0
	 *	            0  0          0         0         1 ]
	 */

	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL  = BCdata->nL;

	double *const dWdW = BCdata->dWBdWL;

	// Standard datatypes
	unsigned int i, iMax, n, eq, var, InddWdW;
	double       n1, n2, n3, *dWdW_ptr[Neq*Neq];
	const double *n_ptr;

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
		dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
	}}
	n_ptr = nL;

	if (d == 3) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;
			n2 = *n_ptr++;
			n3 = *n_ptr++;

			InddWdW = 0;

			// *** var 1 ***
			*dWdW_ptr[InddWdW++] = 1.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 2 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n1*n1;
			*dWdW_ptr[InddWdW++] =     - 2.0*n2*n1;
			*dWdW_ptr[InddWdW++] =     - 2.0*n3*n1;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 3 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] =     - 2.0*n1*n2;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n2*n2;
			*dWdW_ptr[InddWdW++] =     - 2.0*n3*n2;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 4 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] =     - 2.0*n1*n3;
			*dWdW_ptr[InddWdW++] =     - 2.0*n2*n3;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n3*n3;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 5 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0;

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dWdW_ptr[i]++;
		}
	} else if (d == 2) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;
			n2 = *n_ptr++;

			InddWdW = 0;

			// *** var 1 ***
			*dWdW_ptr[InddWdW++] = 1.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 2 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n1*n1;
			*dWdW_ptr[InddWdW++] =     - 2.0*n2*n1;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 3 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] =     - 2.0*n1*n2;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n2*n2;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 4 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0;

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dWdW_ptr[i]++;
		}
	} else if (d == 1) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;

			InddWdW = 0;

			// *** var 1 ***
			*dWdW_ptr[InddWdW++] = 1.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 2 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n1*n1;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 3 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0;

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dWdW_ptr[i]++;
		}
	}
}

static void jacobian_boundary_BackPressure(struct S_BC *const BCdata)
{
	/*
	 *	Jacobian Matrices [var * eq]
	 *
	 *	Supersonic (Outlet)
	 *	dWBdWL = [  1  0  0  0  0
	 *	            0  1  0  0  0
	 *	            0  0  1  0  0
	 *	            0  0  0  1  0
	 *	            0  0  0  0  1 ]
	 *
	 *	Subsonic (Outlet)
	 *		See below.
	 */

	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL  = BCdata->nL,
	             *const WL  = BCdata->WL;

	double *const dWdW = BCdata->dWBdWL;

	// Standard datatypes
	unsigned int n, eq, var;
	double       *dWdW_ptr[Neq*Neq];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *n_ptr, *WL_ptr[Neq];

	// silence
	rhovL_ptr = rhowL_ptr = n_ptr = NULL;

	for (eq  = 0; eq  < Neq;  eq++)  {
		WL_ptr[eq] = &WL[eq*NnTotal];
		for (var = 0; var < Nvar; var++) {
			dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
		}
	}

	double zeros[NnTotal];

	for (n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	rhoL_ptr  = WL_ptr[0];
	rhouL_ptr = WL_ptr[1];
	EL_ptr    = WL_ptr[d+1];

	n_ptr = nL;

	if (d == 3) {
		rhovL_ptr = WL_ptr[2];
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhovL_ptr = WL_ptr[2];
		rhowL_ptr = zeros;
	} else if (d == 1) {
		rhovL_ptr = zeros;
		rhowL_ptr = zeros;
	}

	for (n = 0; n < NnTotal; n++) {
		unsigned int InddWdW = 0;
		double rhoL, rhoL_inv, uL, vL, wL, EL, VL, V2L, pL, cL, c2L;

		// Inner VOLUME
		rhoL     = *rhoL_ptr++;
		rhoL_inv = 1.0/rhoL;

		uL   = (*rhouL_ptr++)*rhoL_inv;
		vL   = (*rhovL_ptr++)*rhoL_inv;
		wL   = (*rhowL_ptr++)*rhoL_inv;
		EL   = *EL_ptr++;

		V2L = uL*uL+vL*vL+wL*wL;
		VL  = sqrt(V2L);

		pL  = GM1*(EL-0.5*rhoL*V2L);


		c2L = GAMMA*pL/rhoL;
		cL  = sqrt(c2L);

/*
		double n1, n2, n3, VnL;

		// silence
		n2 = n3 = 0.0;

		n1 = *n_ptr++;
		if      (d == 3) { n2 = *n_ptr++; n3 = *n_ptr++; }
		else if (d == 2) { n2 = *n_ptr++; n3 = 0.0;      }
		else if (d == 1) { n2 = 0.0;      n3 = 0.0;      }

		VnL = uL*n1+vL*n2+wL*n3;
		if (VnL < 0.0) // Inlet
			printf("\nWarning: Velocity Inflow in jacobian_boundary_BackPressure.\n");
*/

		if (fabs(VL) >= cL) { // Supersonic
			for (var = 0; var < Nvar; var++) {
			for (eq = 0; eq < Neq; eq++) {
				if (var != eq)
					*dWdW_ptr[InddWdW++] = 0.0;
				else
					*dWdW_ptr[InddWdW++] = 1.0;
			}}
		} else {
			double pBack, rho, u, v, w, V2, drhoLdW[Nvar], duLdW[Nvar], dvLdW[Nvar], dwLdW[Nvar], dpLdW[Nvar];

			if (d == 3) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = -wL; dpLdW[4]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
				dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;
			} else if (d == 2) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;
				dwLdW[0] = 0.0;          dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = 0.0;
			} else if (d == 1) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;
				dvLdW[0] = 0.0;          dvLdW[1] = 0.0;      dvLdW[2] = 0.0;      dvLdW[3] = 0.0;
				dwLdW[0] = 0.0;          dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = 0.0;
			}
			for (var = 0; var < Nvar; var++)
				dpLdW[var] *= GM1;

			pBack = DB.pBack;
			rho  = GAMMA*pBack/c2L;
			u    = uL;
			v    = vL;
			w    = wL;
			V2   = V2L;

			for (var = 0; var < Nvar; var++) {
				double drhodW, dudW, dvdW, dwdW, dc2LdW;

				dc2LdW = GAMMA/(rhoL*rhoL)*(dpLdW[var]*rhoL-pL*drhoLdW[var]);
				drhodW = -GAMMA*pBack/(c2L*c2L)*dc2LdW;

				// Note: Using VL for the boundary
				dudW   = duLdW[var];
				dvdW   = dvLdW[var];
				dwdW   = dwLdW[var];

				*dWdW_ptr[InddWdW++] = drhodW;

				*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
				*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
				if (d == 3)
					*dWdW_ptr[InddWdW++] = drhodW*w + rho*dwdW;
				*dWdW_ptr[InddWdW++] = 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
			}
		}

		for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWdW_ptr[i]++;
	}
}

static void jacobian_boundary_Total_TP(struct S_BC *const BCdata)
{

	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL  = BCdata->nL,
	             *const WL  = BCdata->WL;

	double *const dWdW = BCdata->dWBdWL;

	// Initialize DB Parameters
	double Rg      = DB.Rg,
	       p_Total = DB.p_Total,
	       T_Total = DB.T_Total;

	// Standard datatypes
	double       *dWdW_ptr[Neq*Neq];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *n_ptr, *WL_ptr[Neq];

	// silence
	rhowL_ptr = NULL;

	for (size_t eq = 0; eq  < Neq;  eq++)  {
		WL_ptr[eq] = &WL[eq*NnTotal];
		for (size_t var = 0; var < Nvar; var++) {
			dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
		}
	}

	double zeros[NnTotal];

	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	rhoL_ptr  = WL_ptr[0];
	rhouL_ptr = WL_ptr[1];
	rhovL_ptr = WL_ptr[2];
	EL_ptr    = WL_ptr[d+1];

	if (d == 3) {
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhowL_ptr = zeros;
	}

	n_ptr = nL;

	for (size_t n = 0; n < NnTotal; n++) {
		unsigned int InddWdW = 0;
		double       rhoL, rhoL_inv, uL, vL, wL, EL, V2L, pL, cL, HL, VnL, RL, n1, n2, n3;

		// silence
		n3 = 0.0;

		// Inner VOLUME
		rhoL     = *rhoL_ptr++;
		rhoL_inv = 1.0/rhoL;

		uL   = (*rhouL_ptr++)*rhoL_inv;
		vL   = (*rhovL_ptr++)*rhoL_inv;
		wL   = (*rhowL_ptr++)*rhoL_inv;
		EL   = *EL_ptr++;

		V2L = uL*uL+vL*vL+wL*wL;

		pL  = GM1*(EL-0.5*rhoL*V2L);
		cL = sqrt(GAMMA*pL/rhoL);

		HL = (EL+pL)*rhoL_inv;

		n1 = *n_ptr++;
		n2 = *n_ptr++;
		if (d == 3)
			n3 = *n_ptr++;

		VnL = uL*n1+vL*n2+wL*n3;

		RL = VnL + 2.0/GM1*cL;

		// Solve for c
		double aQ, bQ, cQ, term1, term2, cM, cP, cMult, c, Vn, M, T, p, rho, u, v, w;

		// silence
		cMult = 0.0;

		aQ =  1.0 + 2.0/GM1;
		bQ = -2.0*RL;
		cQ =  0.5*GM1*(RL*RL - 2.0*HL);

		term1 = -bQ/(2.0*aQ);
		term2 = sqrt(bQ*bQ-4*aQ*cQ)/(2.0*aQ);

		cM = term1-term2;
		cP = term1+term2;

		// c = max(cM,cP)
		if (cM > cP) {
			c = cM;
			cMult = -1.0;
		} else {
			c = cP;
			cMult = 1.0;
		}

		Vn = RL - 2.0/GM1*c;

//		if (Vn > EPS)
//			printf("\nWarning: Velocity Outflow in jacobian_boundary_Total_TP.\n");

		M = Vn/c;

		T = T_Total/(1+0.5*GM1*M*M);
		p = p_Total*pow(T/T_Total,GAMMA/GM1);

		rho = p/(Rg*T);
		u   = Vn*n1;
		v   = Vn*n2;
		w   = Vn*n3;


		double drhoLdW[Nvar], duLdW[Nvar], dvLdW[Nvar], dwLdW[Nvar], dELdW[Nvar], dpLdW[Nvar];

		if (d == 3) {
			drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
			dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = -wL; dpLdW[4]   = 1.0;
			dELdW[0]   = 0.0;     dELdW[1]   = 0.0; dELdW[2]   = 0.0; dELdW[3]   = 0.0; dELdW[4]   = 1.0;

			duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
			dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
			dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;
		} else if (d == 2) {
			drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
			dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = 1.0;
			dELdW[0]   = 0.0;     dELdW[1]   = 0.0; dELdW[2]   = 0.0; dELdW[3]   = 1.0;

			duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
			dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;
			dwLdW[0] = 0.0;          dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = 0.0;
		}
		for (size_t var = 0; var < Nvar; var++)
			dpLdW[var] *= GM1;

		for (size_t var = 0; var < Nvar; var++) {
			double dcLdW, dHLdW, dVnLdW, dRLdW,
			       dbQdW, dcQdW, dterm1dW, dterm2dW,
			       dcdW, dVndW, dMdW, dTdW, dpdW,
			       drhodW, dudW, dvdW, dwdW, dEdW;

			dcLdW  = 0.5/cL*GAMMA*(dpLdW[var]*rhoL-pL*drhoLdW[var])*(rhoL_inv*rhoL_inv);
			dHLdW  = ((dELdW[var]+dpLdW[var])*rhoL-(EL+pL)*drhoLdW[var])*(rhoL_inv*rhoL_inv);
			dVnLdW = duLdW[var]*n1+dvLdW[var]*n2+dwLdW[var]*n3;

			dRLdW = dVnLdW + 2.0/GM1*dcLdW;

//			daQdW =  0.0;
			dbQdW = -2.0*dRLdW;
			dcQdW =  0.5*GM1*(2.0*RL*dRLdW - 2.0*dHLdW);

			dterm1dW = -dbQdW/(2.0*aQ);
			dterm2dW = 0.5/sqrt(bQ*bQ-4.0*aQ*cQ)*(2.0*bQ*dbQdW-4.0*aQ*dcQdW)/(2.0*aQ);

			dcdW  = dterm1dW + cMult*dterm2dW;
			dVndW = dRLdW - 2.0/GM1*dcdW;

			dMdW = dVndW/c-Vn/(c*c)*dcdW;

			dTdW = T_Total*(-1.0)*pow(1+0.5*GM1*M*M,-2.0)*0.5*GM1*2.0*M*dMdW;
			dpdW = p_Total*GAMMA/GM1*pow(T/T_Total,GAMMA/GM1-1.0)*dTdW/T_Total;

			drhodW = dpdW/(Rg*T)-p/(Rg*T*T)*dTdW;
			dudW   = dVndW*n1;
			dvdW   = dVndW*n2;
			dwdW   = dVndW*n3;
			dEdW   = dpdW/GM1+0.5*(drhodW*(u*u+v*v+w*w)+rho*2.0*(u*dudW+v*dvdW+w*dwdW));

			*dWdW_ptr[InddWdW++] = drhodW;

			*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
			*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
			if (d == 3)
				*dWdW_ptr[InddWdW++] = drhodW*w + rho*dwdW;
			*dWdW_ptr[InddWdW++] = dEdW;
		}

		for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWdW_ptr[i]++;
	}
}

static void jacobian_boundary_SupersonicInflow(struct S_BC *const BCdata)
{
	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double *const dWdW = BCdata->dWBdWL;

	// Standard datatypes
	double       *dWdW_ptr[Neq*Neq];

	for (size_t eq  = 0; eq  < Neq;  eq++)  {
	for (size_t var = 0; var < Nvar; var++) {
		dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
	}}

	for (size_t n = 0; n < NnTotal; n++) {
		unsigned int InddWdW = 0;
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t eq = 0; eq < Neq; eq++) {
			*dWdW_ptr[InddWdW++] = 0.0;
		}}

		for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWdW_ptr[i]++;
	}
}

static void jacobian_boundary_SupersonicOutflow(struct S_BC *const BCdata)
{
	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double *const dWdW = BCdata->dWBdWL;

	// Standard datatypes
	double       *dWdW_ptr[Neq*Neq];

	for (size_t eq  = 0; eq  < Neq;  eq++)  {
	for (size_t var = 0; var < Nvar; var++) {
		dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
	}}

	for (size_t n = 0; n < NnTotal; n++) {
		unsigned int InddWdW = 0;
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t eq = 0; eq < Neq; eq++) {
			if (var != eq)
				*dWdW_ptr[InddWdW++] = 0.0;
			else
				*dWdW_ptr[InddWdW++] = 1.0;
		}}

		for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWdW_ptr[i]++;
	}
}

static void jacobian_boundary_NoSlip_Dirichlet(struct S_BC *const BCdata)
{

	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double *const dWBdWL = BCdata->dWBdWL,
	       *      dWBdWL_ptr[Neq*Nvar];

	for (size_t eq  = 0; eq  < Neq;  eq++)  {
	for (size_t var = 0; var < Nvar; var++) {
		dWBdWL_ptr[eq*Nvar+var] = &dWBdWL[(eq*Nvar+var)*NnTotal];
	}}

	for (size_t n = 0; n < NnTotal; n++) {
		size_t InddWdW = 0;
		for (size_t var = 0; var < Nvar; var++) {
			if (0) {
				EXIT_BASIC; // Update this if the entropy variable formulation is used. See boundary_NoSlip_Dirichlet.
			} else {
				// When all boundary conditions are imposed.
				for (size_t eq = 0; eq < Neq; eq++) {
					if (var != eq)
						*dWBdWL_ptr[InddWdW++] = 0.0;
					else
						*dWBdWL_ptr[InddWdW++] = -1.0;
				}
			}
		}

		for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWBdWL_ptr[i]++;
	}
}

static void jacobian_boundary_NoSlip_Adiabatic(struct S_BC *const BCdata)
{

	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Neq     = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double const *const WL = BCdata->WL;

	double *const dWBdWL = BCdata->dWBdWL,
	       *      dWBdWL_ptr[Neq*Nvar];

	double const *WL_ptr[Nvar];

	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
	}

	double zeros[NnTotal];
	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	double const *rhoL_ptr  = WL_ptr[0],
	             *rhouL_ptr = WL_ptr[1],
	             *rhovL_ptr = WL_ptr[2],
	             *rhowL_ptr            ;

	if (d == 3) {
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhowL_ptr = zeros;
	} else {
		EXIT_UNSUPPORTED;
	}

	for (size_t eq  = 0; eq  < Neq;  eq++)  {
	for (size_t var = 0; var < Nvar; var++) {
		dWBdWL_ptr[eq*Nvar+var] = &dWBdWL[(eq*Nvar+var)*NnTotal];
	}}

	for (size_t n = 0; n < NnTotal; n++) {
		double const rhoL     = *rhoL_ptr++,
		             rhoL_inv = 1.0/rhoL,

		             uL   = (*rhouL_ptr++)*rhoL_inv,
		             vL   = (*rhovL_ptr++)*rhoL_inv,
		             wL   = (*rhowL_ptr++)*rhoL_inv,

		             V2L = uL*uL+vL*vL+wL*wL;

		double drhoLdW[Nvar], duLdW[Nvar], dvLdW[Nvar], dwLdW[Nvar], dELdW[Nvar], dpLdW[Nvar];
		if (d == 3) {
			drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
			dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = -wL; dpLdW[4]   = 1.0;
			dELdW[0]   = 0.0;     dELdW[1]   = 0.0; dELdW[2]   = 0.0; dELdW[3]   = 0.0; dELdW[4]   = 1.0;

			duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
			dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
			dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;
		} else if (d == 2) {
			drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
			dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = 1.0;
			dELdW[0]   = 0.0;     dELdW[1]   = 0.0; dELdW[2]   = 0.0; dELdW[3]   = 1.0;

			duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
			dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;
			dwLdW[0] = 0.0;          dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = 0.0;
		} else {
			EXIT_UNSUPPORTED;
		}
		for (size_t var = 0; var < Nvar; var++)
			dpLdW[var] *= GM1;

		double u = 0.0, v = 0.0, w = 0.0;
		if (strstr(DB.TestCase,"PlaneCouette") ||
		    strstr(DB.TestCase,"TaylorCouette")) {
			; // Do nothing
		} else {
			EXIT_UNSUPPORTED;
		}

		size_t InddWdW = 0;
		for (size_t var = 0; var < Nvar; var++) {
			*dWBdWL_ptr[InddWdW++] = drhoLdW[var];
			*dWBdWL_ptr[InddWdW++] = -(drhoLdW[var]*uL+rhoL*duLdW[var]) + 2.0*drhoLdW[var]*u;
			*dWBdWL_ptr[InddWdW++] = -(drhoLdW[var]*vL+rhoL*dvLdW[var]) + 2.0*drhoLdW[var]*v;
			if (d == 3)
				*dWBdWL_ptr[InddWdW++] = -(drhoLdW[var]*wL+rhoL*dwLdW[var]) + 2.0*drhoLdW[var]*w;
			*dWBdWL_ptr[InddWdW++] = dELdW[var];
		}

		for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWBdWL_ptr[i]++;
	}
}
