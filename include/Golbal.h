// global.h
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdbool.h>

// Define test case enumeration type in header file
typedef enum {
    TEST_PRECISION,             //
    TEST_1D_SHOCKTUBE,          // 1D Sod shock tube problem - classical benchmark
    TEST_1D_CONTACTWAVE,        // 1D Contact discontinuity problem - for testing contact resolution
    TEST_1D_SHOCKIMPACT,        // 1D Symmetric shock collision problem - tests overheating behavior
    TEST_1D_IMPACTWALL,         // 1D Shock reflection from rigid wall - wall heating test
    TEST_1D_DOUBLERARE,         // 1D Double rarefaction wave problem - low density test
    TEST_1D_NOHPROBLEM,         // 1D Noh problem - strong shock implosion test
    TEST_2D_SHOCKTUBE_CASE1,    // 2D Riemann problem Case 1 - Configuration 1 of 2D shock tube
    TEST_2D_SHOCKTUBE_CASE2,    // 2D Riemann problem Case 2 - Configuration 2 of 2D shock tube
    TEST_2D_SHOCKTUBE_CASE3,    // 2D Riemann problem Case 3 - Configuration 3 of 2D shock tube
    TEST_2D_SHOCKTUBE_CASE4,    // 2D Riemann problem Case 4 - Configuration 4 of 2D shock tube
    TEST_2D_SHOCKTUBE_CASE5,    // 2D Riemann problem Case 5 - Configuration 5 of 2D shock tube
    TEST_TAYLOR_GREEN_VORTEX,   // Taylor-Green vortex - decaying turbulence benchmark
    TEST_GAUSSIAN_PULSE,        // Gaussian pulse - smooth flow test case
    TEST_KELVIN_HELMHOLTZ,      // Kelvin-Helmholtz instability - shear layer instability
    TEST_RAYLEIGH_TAYLOR,       // Rayleigh-Taylor instability - buoyancy-driven instability
    TEST_DOUBLE_MACH_REFLECTION,// Double Mach reflection - complex shock interaction test
    TEST_BLAST_WAVE,            // Spherical blast wave - multidimensional strong shock
    TEST_NOH_PROBLEM,           // 2D/3D Noh problem - multidimensional implosion test
    TEST_OddEven_Decoupling,    // Odd-Even Decoupling - carbuncle phenomenon test
    TEST_BACKWARD_STEP          // Backward step flow - separation and recirculation test
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
static double CFL = 0.5;
// Reconstruction method control variable {0 is 0th order, 2 is 2nd order TVD scheme; 3 is 3rd order WENO, 5 is 5th order WENO reconstruction}
// TVD includes Vanleer Limiter, Minmod limiter, etc. Modify in CFD_convection.h: 
static int Recon_Accur = 5;
// Whether to use characteristic reconstruction
static bool Characteriz = false;
// Whether to enable gravity:
static bool Source;
static double Gravity;

// Computational output parameters, etc.
static int Control_Compution;
static int Control_output;

// Material property parameters:
static double M_gamma;
#endif