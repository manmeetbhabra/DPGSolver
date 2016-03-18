/*
  Purpose:
    Define parameters
 
  Comments:
 
  Notation:
 
  References:
 
*/

#ifndef DPG__parameters_h__INCLUDED
#define DPG__parameters_h__INCLUDED

// Common variables
#define PI    3.1415926535897932
#define GAMMA 1.4

// Element types (Gmsh convention)
#define POINT 15
#define LINE  1
#define TRI   2
#define QUAD  3
#define TET   4
#define HEX   5
#define WEDGE 6
#define PYR   7

/*  Node Tolerance (for physical coordinate comparison)  */
#define NODETOL 1.e-12

/*  Value close to machine zero  */
#define EPS     1.e-15

/*  Min/Max String Lengths  */
#define STRLEN_MIN 64
#define STRLEN_MAX 512




#endif // DPG__parameters_h_INCLUDED
