// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__solver_poisson_c_h__INCLUDED
#define DPG__solver_poisson_c_h__INCLUDED

extern void compute_qhat_VOLUME_c (void);
extern void compute_qhat_FACET_c  (void);
extern void finalize_qhat_c       (void);
extern void compute_uhat_VOLUME_c (void);
extern void compute_uhat_FACET_c  (void);

#endif // DPG__solver_poisson_c_h__INCLUDED