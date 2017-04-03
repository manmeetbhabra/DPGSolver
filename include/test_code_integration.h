// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_h__INCLUDED
#define DPG__test_code_integration_h__INCLUDED

#include <stdbool.h>

extern void code_startup             (int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv);
extern void code_cleanup             (void);
extern void code_startup_mod_prmtrs  (int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv,
                                      const unsigned int phase);
extern void code_startup_mod_ctrl    (int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv,
                                      const unsigned int phase);
extern void evaluate_mesh_regularity (double *mesh_quality);
extern void check_convergence_orders (const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                                      const unsigned int PMax, unsigned int *pass, const bool PrintEnabled);
extern void check_mesh_regularity    (const double *mesh_quality, const unsigned int NML, unsigned int *pass,
                                      const bool PrintEnabled);
extern void set_PrintName            (char *name_type, char *PrintName, bool *TestTRI);

#endif // DPG__test_code_integration_h__INCLUDED
