// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__element_functions_h__INCLUDED
#define DPG__element_functions_h__INCLUDED

extern void             initialize_ELEMENTs (void);
extern void             finalize_ELEMENTs   (void);
extern unsigned int     get_IndFType        (const unsigned int Eclass, const unsigned int f);
extern unsigned int     is_ELEMENT_present  (const unsigned int type);
extern unsigned int     get_Eclass          (const unsigned int type);
extern struct S_ELEMENT *get_ELEMENT_type   (const unsigned int type);
extern struct S_ELEMENT *get_ELEMENT_F_type (const unsigned int type, const unsigned int f);
extern struct S_ELEMENT *get_ELEMENT_Eclass (const unsigned int type, const unsigned int IndEclass);
extern struct S_ELEMENT *get_ELEMENT_FACE  (const unsigned int type, const unsigned int IndEclass);

#endif // DPG__element_functions_h__INCLUDED
