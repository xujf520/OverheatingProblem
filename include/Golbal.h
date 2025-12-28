// global.h
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdbool.h>

// Define test case enumeration type in header file
typedef enum {
    TEST_1D_SHOCKTUBE,          // 1D shock tube (Sod problem)
    TEST_2D_SHOCKTUBE_CASE1,    // 2D Riemann problem Case 1
    TEST_2D_SHOCKTUBE_CASE2,    // 2D Riemann problem Case 2
    TEST_2D_SHOCKTUBE_CASE3,    // 2D Riemann problem Case 3
    TEST_2D_SHOCKTUBE_CASE4,    // 2D Riemann problem Case 4
    TEST_2D_SHOCKTUBE_CASE5,    // 2D Riemann problem Case 5
    TEST_TAYLOR_GREEN_VORTEX,   // Taylor-Green vortex
    TEST_GAUSSIAN_PULSE,        // Gaussian pulse
    TEST_KELVIN_HELMHOLTZ,      // Kelvin-Helmholtz instability
    TEST_RAYLEIGH_TAYLOR,       // Rayleigh-Taylor instability
    TEST_DOUBLE_MACH_REFLECTION,// Double Mach reflection
    TEST_BLAST_WAVE,           // Spherical blast wave
    TEST_NOH_PROBLEM,           // Noh problem
    TEST_BACKWARD_STEP         // Backward step flow
} TestCase2D;

// Computational control variables
const int var = 4;    
//static int Control_Compution = 0;
//static int Control_output = 2;
static bool Control_Out = false;

const int Block = 2;

// Boundary condition control:
//static bool Periodicity = false;
//static bool Reflect = false;

// Program computational parameters
static double Time = 0;       // Computational domain parameter 
static int Ite = 0;

// Ghost cell parameter
static int GhostCell = 4;

/*                                            *********                                          */
/*                                           Declare global variables                           */
/*                                            *********                                          */
// Program computational parameters
static double CFL = 0.1;
// Reconstruction method control variable {0 is 0th order, 2 is 2nd order TVD scheme; 3 is 3rd order WENO, 5 is 5th order WENO reconstruction}
// TVD includes Vanleer Limiter, Minmod limiter, etc. Modify in CFD_convection.h: 
static int Recon_Accur = 5;
// Whether to use characteristic reconstruction
static bool Characteriz = true;
// Whether to enable gravity:
static bool Source;
static double Gravity;

// Computational output parameters, etc.
static int Control_Compution;
static int Control_output;

// Material property parameters:
static double M_gamma;
#endif