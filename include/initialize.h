#ifndef INITIALIZATION_H
#define INITIALIZATION_H


#include <stdio.h>
#include "Golbal.h"
#include "math.h"
#include "function.h"
#include "Boundary_Condition.h"


// Define The Constant Pi
#ifndef PI
    #define PI 3.14159265358979323846
#endif

static inline void initEuler2D(int var, int rows, int cols, double (*x)[rows][cols]) {
    int i, j, k;

    #pragma omp parallel for collapse(3)
    for (i = 0; i < var; i++) 
        for (j = 0; j < rows; j++) 
           for (k = 0; k < cols; k++) 
                x[i][j][k] = 0.0;

}

//==============================================================================
// 1D Euler Equations with Exact Solution Test
// Initial conditions: ρ(x,0) = 1 + 0.2*sin(x), u(x,0) = 1, p(x,0) = 1
// Exact solution: ρ(x,t) = 1 + 0.2*sin(x - t)
// Domain: [0, 2π]
// Boundary conditions: Periodic
// Final time: t = 2.0
//==============================================================================
static inline void initEulerpri1D_PrecisionTest(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    // Set The Boundary Condition To Periodic
    bc_config.left = BC_PERIODICITY;      
    bc_config.right = BC_PERIODICITY;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 2.0 * M_PI;    // 计算域长度 [0, 2π]
    *Ly = 1.0;           // y方向长度
    
    int nx = rows - 2 * GhostCell;
    double dx = *Lx / nx;
   
    //#pragma omp parallel
    for (i = GhostCell; i < rows-GhostCell; i++) {
        x[0][i][GhostCell] = 1.0 + 0.2 * sin((i - GhostCell + 0.5) * dx);
        x[1][i][GhostCell] = 1.0;
        x[2][i][GhostCell] = 0.0;
        x[3][i][GhostCell] = 1.0;
    }

    *Time = 2.0;
}

//==============================================================================
// 1D Shock Tube Test (Sod Problem)
// Left state:  (ρ, u, v, p) = (1, 0, 0, 1)
// Right state: (ρ, u, v, p) = (0.125, 0, 0, 0.1)
// Domain: [0,1] × [0,1]
// Boundary conditions: All outflow
// Reference time: t = 0.2
//==============================================================================
static inline void initEulerpri1D_Shocktube(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Shock Tube Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;    
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 ){
                x[0][i][j] = 1.0;
                x[1][i][j] = 0.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
            else {
                x[0][i][j] = 0.125;
                x[1][i][j] = 0.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 0.1;
            }
        }
    }

    *Time = 0.2;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

//==============================================================================
// 1D Contact Discontinuity Test
// Left state:  (ρ, u, v, p) = (1, 0, 0, 1)
// Right state: (ρ, u, v, p) = (0.125, 0, 0, 1)
// Domain: [0,1] × [0,1]
// Boundary conditions: All outflow
// Reference time: t = 0.2
// Purpose: Tests numerical diffusion at contact interfaces
//==============================================================================
static inline void initEulerpri1D_ContactWave(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Contact Discontinuity Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;    
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 ){
                x[0][i][j] = 1.0;
                x[1][i][j] = 0.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
            else {
                x[0][i][j] = 0.125;
                x[1][i][j] = 0.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
        }
    }

    *Time = 0.2;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

//==============================================================================
// 1D Shock Impact (Symmetric Shock Collision) Test
// Left state:  (ρ, u, v, p) = (1, 4, 0, 1)
// Right state: (ρ, u, v, p) = (1, -4, 0, 1)
// Domain: [0,1] × [0,1]
// Boundary conditions: All outflow
// Reference time: t = 0.2
// Purpose: Tests overheating behavior in shock collision
//==============================================================================
static inline void initEulerpri1D_ShockImpact(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Symmetric Shock Collision Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;    
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 ){
                x[0][i][j] = 1.0;
                x[1][i][j] = 4.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
            else {
                x[0][i][j] = 1.0;
                x[1][i][j] = -4.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
        }
    }

    *Time = 0.2;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

//==============================================================================
// 1D Double Rarefaction Wave Test
// Left state:  (ρ, u, v, p) = (1, -2, 0, 0.4)
// Right state: (ρ, u, v, p) = (1, 2, 0, 0.4)
// Domain: [0,1] × [0,1]
// Boundary conditions: All outflow
// Reference time: t = 0.15
// Purpose: Tests performance in low-density and vacuum regions
//==============================================================================
static inline void initEulerpri1D_DoubleRare(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Double Rarefaction Wave Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;    
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 ){
                x[0][i][j] = 1.0;
                x[1][i][j] = -2.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
            else {
                x[0][i][j] = 1.0;
                x[1][i][j] = 2.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
        }
    }

    *Time = 0.15;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

//==============================================================================
// 1D Shock-Wall Impact (Shock Reflection) Test
// Initial state: (ρ, u, v, p) = (1, 1, 0, 0.001)
// Domain: [0,1] × [0,1]
// Boundary conditions: Left/outflow, Right/reflecting wall
// Reference time: t = 0.6
// Purpose: Tests wall overheating phenomenon in shock capturing
//==============================================================================
static inline void initEulerpri1D_ShockImpactWall(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Shock Reflection from Rigid Wall Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_REFLECTION;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;    
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            x[0][i][j] = 1.0;
            x[1][i][j] = 1.0;
            x[2][i][j] = 0.0;
            x[3][i][j] = 1.0;
        }
    }

    *Time = 0.6;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  ReFelection\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

//==============================================================================
// 1D Noh Problem (Strong Implosion Test)
// Initial state: (ρ, u, v, p) = (1, -1, 0, 1e-6)
// Domain: [0,1] × [0,1]
// Boundary conditions: Left/outflow, Right/reflecting wall
// Reference time: t = 0.6
// Purpose: Tests strong converging shocks and analytical solution
//==============================================================================
static inline void initEulerpri1D_NohProblem(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Noh Implosion Problem Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_REFLECTION;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;
    bc_config.top = BC_OUTFLOW;

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;    
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            x[0][i][j] = 1.0;
            x[1][i][j] = -1.0;
            x[2][i][j] = 0.0;
            x[3][i][j] = 1e-6;
        }
    }

    *Time = 0.6;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Reflection\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}


// Case 1
static inline void initEulerpri2D_Shocktube1(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 1 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;    
    bc_config.top = BC_OUTFLOW;       

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0; 

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 2;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 3;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
        }
    }
    *Time = 0.23;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 2
static inline void initEulerpri2D_Shocktube2(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    
    printf("=== 2D Shock Tube Case 2 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;    
    bc_config.top = BC_OUTFLOW;       

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0; 

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 2;       // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = -0.3;    // v
                x[3][i][j] = 1;       // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1;       // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = -0.4;    // v
                x[3][i][j] = 1;       // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1.0625;  // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = 0.2145;  // v
                x[3][i][j] = 0.4;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 0.5197;  // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = -1.125;  // v
                x[3][i][j] = 0.4;     // p
            }
        }
    }

    *Time = 0.30;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 3
static inline void initEulerpri2D_Shocktube3(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 3 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;    
    bc_config.top = BC_OUTFLOW;       

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0; 

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 0.5065;  // ρ
                x[1][i][j] = 0.8939;  // u
                x[2][i][j] = 0;       // v
                x[3][i][j] = 0.35;    // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1.5;     // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = 0;       // v
                x[3][i][j] = 1.5;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1.1;     // ρ
                x[1][i][j] = 0.8939;  // u
                x[2][i][j] = 0.8939;  // v
                x[3][i][j] = 1.1;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 0.5065;  // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = 0.8939;  // v
                x[3][i][j] = 0.35;    // p
            }
        }
    }

    *Time = 0.25;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 4
static inline void initEulerpri2D_Shocktube4(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 4 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;    
    bc_config.top = BC_OUTFLOW;       

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0; 

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 2;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 3;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
        }
    }

    *Time = 0.30;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 5
static inline void initEulerpri2D_Shocktube5(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 5 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;    
    bc_config.top = BC_OUTFLOW;       

    int i, j;
    *Lx = 1.0;    
    *Ly = 1.0;

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 1;       // ρ
                x[1][i][j] = -0.6259; // u
                x[2][i][j] = 0.1;     // v
                x[3][i][j] = 1;       // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 0.5197;  // ρ
                x[1][i][j] = 0.1;     // u
                x[2][i][j] = 0.1;     // v
                x[3][i][j] = 0.4;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 0.8;     // ρ
                x[1][i][j] = 0.1;     // u
                x[2][i][j] = 0.1;     // v
                x[3][i][j] = 1;       // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 1;       // ρ
                x[1][i][j] = 0.1;     // u
                x[2][i][j] = -0.6259; // v
                x[3][i][j] = 1;       // p
            }
        }
    }

    *Time = 0.25;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}


// Taylor Green Vortex Problem
static inline void initEulerpri2D_TaylorGreenVortex(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Taylor-Green Vortex Test Case ===\n");

    bc_config.left = BC_PERIODICITY;   
    bc_config.right = BC_PERIODICITY;
    bc_config.bottom = BC_PERIODICITY;
    bc_config.top = BC_PERIODICITY;

    *Lx = 2.0 * M_PI;  
    *Ly = 2.0 * M_PI;
    double rho0 = 1.0;       // Reference Density
    double p0 = 1.0;         // Reference Pressure
    double U0 = 1.0;         // Reference Speed
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    
    // Precomputed Constant
    double inv_nx = 1.0 / nx;
    double inv_ny = 1.0 / ny;
    double prefactor = rho0 * U0 * U0 / 4.0;
    
    //Initialize The Entire Grid In Parallel
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Calculate Physical Coordinate Position
            double x_pos = (double)i * inv_nx * (*Lx);
            double y_pos = (double)j * inv_ny * (*Ly);
            
            // Taylor Green Vortex Velocity Field
            double sin_x = sin(x_pos);
            double cos_x = cos(x_pos);
            double sin_y = sin(y_pos);
            double cos_y = cos(y_pos);
            
            double u = U0 * sin_x * cos_y;      
            double v = -U0 * cos_x * sin_y;    
            double pressure = p0 + prefactor * (cos(2.0 * x_pos) + cos(2.0 * y_pos));
            
            x[0][i][j] = rho0;   
            x[1][i][j] = u;       
            x[2][i][j] = v;       
            x[3][i][j] = pressure;
        }
    }

    *Time = 5.0;

    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Periodic\n");
    printf("  Right:  Periodic\n");
    printf("  Bottom: Periodic\n");
    printf("  Top:    Periodic\n");
}

// Two Dimensional Gaussian Density Pulse Problem All Boundaries Are Outlet Boundaries
static inline void initEulerpri2D_GaussianPulse(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== 2D Gaussian Density Pulse Test Case ===\n");

    // Set boundary conditions for open domain
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_OUTFLOW;    
    bc_config.top = BC_OUTFLOW;       

    int i, j;
    *Lx = 10.0;    
    *Ly = 10.0; 
    double x0 = (*Lx)/2.0;  // Pulse Center X Coordinate
    double y0 = (*Ly)/2.0;  // Pulse Center Y Coordinate
    double A = 1.0;      // Pulse Amplitude
    double sigma = 1.0;  // Pulse Width
    double rho0 = 1.0;   // Background Density
    double p0 = 1.0;     // Background Pressure
    double u0 = 0.5;     
    double v0 = 0.3;     
    double x_pos, y_pos, r2;
    double nx = rows - 2* GhostCell;
    double ny = cols - 2* GhostCell;


    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            x_pos = (double)i / nx * (*Lx);
            y_pos = (double)j / ny * (*Ly);
            // Square Of The Distance To The Pulse Center
            r2 = (x_pos - x0)*(x_pos - x0) + (y_pos - y0)*(y_pos - y0);
            // Set Initial Conditions Gaussian Density Pulse With Uniform Background Flow
            x[0][i][j] = rho0 + A * exp(-r2 / (sigma * sigma));  
            x[1][i][j] = u0;      
            x[2][i][j] = v0;      
            x[3][i][j] = p0;  
        }
    }

    *Time = 4.0;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Gaussian Pulse: Center at (%f, %f), Amplitude = %f\n", x0, y0, A);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}


// KH Instability with Sharp Interface (ICs A - Based on Springel 2009)
static inline void initEulerpri2D_KelvinHelmholtz(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== KH Instability - Sharp Interface (ICs A - Springel 2009) ===\n");
    
    // Boundary conditions: Periodic in all directions for standard KH test
    bc_config.left = BC_PERIODICITY;
    bc_config.right = BC_PERIODICITY;
    bc_config.bottom = BC_PERIODICITY;  // Changed from reflection to periodic
    bc_config.top = BC_PERIODICITY;     // Changed from reflection to periodic
    
    int i, j;
    *Lx = 1.0;    // Domain length in x-direction
    *Ly = 1.0;    // Domain length in y-direction (back to 1.0)
    
    // Fluid properties (from Springel 2009)
    double rho1 = 2.0;     // Central slab density (|y-0.5| < 0.25)
    double rho2 = 1.0;     // Outer regions density (|y-0.5| > 0.25)
    double u1 = 0.5;       // Central slab x-velocity
    double u2 = -0.5;      // Outer regions x-velocity (opposite direction)
    double p0 = 2.5;       // Uniform pressure
    double gamma = 5.0/3.0; // Adiabatic index
    
    // Perturbation parameters (from equation 4 in the paper)
    double w0 = 0.1;       // Perturbation amplitude
    int n = 4;             // Mode number (n=4 for four wavelengths in domain)
    double sigma = 0.05;   // Gaussian width for perturbation
    
    // Physical grid dimensions (excluding ghost cells)
    int nx_phys = rows - 2 * GhostCell;
    int ny_phys = cols - 2 * GhostCell;
    
    // Initialize flow field - only physical cells
    for (i = GhostCell; i < rows - GhostCell; i++) {
        for (j = GhostCell; j < cols - GhostCell; j++) {
            // Calculate physical coordinates (0 to Lx, 0 to Ly)
            // Using cell-centered coordinates
            double x_pos = (double)(i - GhostCell + 0.5) / nx_phys * (*Lx);
            double y_pos = (double)(j - GhostCell + 0.5) / ny_phys * (*Ly);
            
            // SHARP INTERFACE: Central slab at |y-0.5| < 0.25
            // This is the key feature of ICs A - discontinuous interface
            if (fabs(y_pos - 0.5) < 0.25) {
                // Central slab (denser fluid)
                x[0][i][j] = rho1;  // Density = 2.0
                x[1][i][j] = u1;    // x-velocity = 0.5 (to the right)
            } else {
                // Outer regions (lighter fluid)
                x[0][i][j] = rho2;  // Density = 1.0
                x[1][i][j] = u2;    // x-velocity = -0.5 (to the left)
            }
            
            // VERTICAL VELOCITY PERTURBATION (equation 4 from paper)
            // Original form: v_y(x,y) = w0 * sin(nπx) * 
            // [exp(-(y-0.25)^2/(2σ^2)) + exp(-(y-0.75)^2/(2σ^2))]
            double vy_pert = w0 * sin(n * M_PI * x_pos);
            double gaussian1 = exp(-pow(y_pos - 0.25, 2) / (2.0 * sigma * sigma));
            double gaussian2 = exp(-pow(y_pos - 0.75, 2) / (2.0 * sigma * sigma));
            x[2][i][j] = vy_pert * (gaussian1 + gaussian2);  // y-velocity perturbation
            
            x[3][i][j] = p0;  // Pressure (uniform throughout domain)
        }
    }
    
    *Time = 2.0;  // Simulation time as in Springel (2009)
    
    printf("Computational domain: Lx = %.2f, Ly = %.2f\n", *Lx, *Ly);
    printf("Grid configuration:\n");
    printf("  Total cells: %d x %d (including ghost cells)\n", rows, cols);
    printf("  Physical cells: %d x %d\n", nx_phys, ny_phys);
    printf("  Ghost cells: %d per side\n", GhostCell);
    printf("Fluid Properties (Springel 2009 setup):\n");
    printf("  Central slab (|y-0.5| < 0.25): ρ = %.1f, u = %.1f\n", rho1, u1);
    printf("  Outer regions (|y-0.5| > 0.25): ρ = %.1f, u = %.1f\n", rho2, u2);
    printf("  Velocity difference: Δu = u1 - u2 = %.1f - (%.1f) = %.1f\n", u1, u2, u1 - u2);
    printf("  Density ratio: ρ1/ρ2 = %.1f/%.1f = %.1f\n", rho1, rho2, rho1/rho2);
    printf("  Pressure: p = %.1f (uniform)\n", p0);
    printf("  Adiabatic index: γ = %.3f\n", gamma);
    printf("Mach numbers:\n");
    double a1 = sqrt(gamma * p0 / rho1);  // Sound speed in dense fluid
    double a2 = sqrt(gamma * p0 / rho2);  // Sound speed in light fluid
    printf("  Central slab: M1 = |u1|/a1 = %.2f/%.2f = %.2f\n", fabs(u1), a1, fabs(u1)/a1);
    printf("  Outer regions: M2 = |u2|/a2 = %.2f/%.2f = %.2f\n", fabs(u2), a2, fabs(u2)/a2);
    printf("Perturbation (equation 4):\n");
    printf("  Amplitude: w0 = %.2f\n", w0);
    printf("  Mode number: n = %d (four wavelengths in domain)\n", n);
    printf("  Gaussian width: σ = %.3f\n", sigma);
    printf("  Wavenumber: kx = nπ = %.3f\n", n * M_PI);
    printf("  Wavelength: λ = 2/n = %.2f\n", 2.0/n);
    printf("Final simulation time: t = %.1f (as in Springel 2009)\n", *Time);
    printf("Boundary Conditions: Periodic in both x and y directions\n");
    printf("\nKey features of ICs A (Springel 2009):\n");
    printf("1. SHARP interface at |y-0.5| = 0.25 (discontinuous density and velocity)\n");
    printf("2. Central slab: ρ=2.0, u=0.5; Outer regions: ρ=1.0, u=-0.5\n");
    printf("3. Perturbation localized at both interfaces (y=0.25 and y=0.75)\n");
    printf("4. Uniform pressure: p=2.5, γ=5/3\n");
    printf("5. Shown in Robertson et al. (2010) to exhibit Galilean non-invariance at low resolution\n");
    printf("   due to numerical diffusion affecting small-scale modes\n");
}


static inline void initEulerpri2D_RayleighTaylor(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Rayleigh-Taylor Instability Test Case ===\n");

     // Set Boundary Conditions
    bc_config.left = BC_REFLECTION;     
    bc_config.right = BC_REFLECTION;
    bc_config.bottom = BC_FIXED_VALUE;  
    bc_config.top = BC_FIXED_VALUE;

    // Set Computational Domain
    *Lx = 0.25;
    *Ly = 1.0;
    *Time = 1.95;

    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;

    //#pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // Computational Physical Coordinates
            double x_pos = (double)(i + 0.5 - GhostCell) / nx * (*Lx);
            double y_pos = (double)(j + 0.5 - GhostCell) / ny * (*Ly);
            
            // Set Different Initial States Based On The Y Coordinate Position
            if (y_pos < 0.5) {
                // Lower Fluid 0 < y < 1/2
                x[0][i][j] = 2.0;                           
                x[1][i][j] = 0.0;                          
                x[2][i][j] = -0.025 * sqrt((2.0 * y_pos + 1.0)*M_gamma/2.0) \
                                    * cos(8.0 * M_PI * x_pos);
                x[3][i][j] = 2.0 * y_pos + 1.0;            
            } else {
                // Upper Fluid 1/2 < y ≤ 1
                x[0][i][j] = 1.0;                            
                x[1][i][j] = 0.0;                          
                x[2][i][j] = -0.025 * sqrt((y_pos + 1.5)*M_gamma/1.0) \
                                    * cos(8.0 * M_PI * x_pos); 
                x[3][i][j] = y_pos + 1.5;                  
            }
        }
    }

    printf("Computational domain: [0, %f] x [0, %f]\n", *Lx, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Density ratio: 2:1 (heavy fluid below, light fluid above)\n");
    printf("Initial perturbation: v = -0.025 cos(8πx)\n");
    printf("Gravity direction: -y (implied by pressure gradient)\n");

}


// Two Dimensional Double Mach Reflection Problem
// Double Mach Reflection Initialization Function
static inline void initEulerpri2D_DoubleMachReflection(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    double post_shock_state[4] = {8.0, 7.145, -4.125, 116.83333};
    double pre_shock_state[4] = {1.4, 0.0, 0.0, 1.0};
    double shock_slope = sqrt(3.0);  // tan(60°)
    double shock_start_x = 1.0/6.0;
    double current_time = 0.0;                                  // The Current Time Needs To Be Updated In The Main Loop
    printf("=== Double Mach Reflection Test Case ===\n");

    // Set Boundary Conditions
    bc_config.left = BC_INFLOW;                
    bc_config.right = BC_OUTFLOW;               
    bc_config.bottom = BC_REFLECTION;   
    bc_config.top = BC_REFLECTION;         

    // Set Computational Domain
    *Lx = 4.0;
    *Ly = 1.0;
    *Time = 0.2;

    // Set Shockwave State 
    post_shock_state[0] = 8.0;        
    post_shock_state[1] = 7.1447;        
    post_shock_state[2] = -4.125;       
    post_shock_state[3] = 116.5;    

    pre_shock_state[0] = 1.4;         
    pre_shock_state[1] = 0.0;           
    pre_shock_state[2] = 0.0;           
    pre_shock_state[3] = 1.0;           

    //shock_slope = 1.732;          
    shock_start_x = 1.0/6.0;

    // Initialize The Flow Field
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows-GhostCell; i++) {
        for (int j = GhostCell; j < cols-GhostCell; j++) {
            // Calculate Physical Coordinates
            double x_pos = (double)(i + 0.5 - GhostCell) / nx * (*Lx);
            double y_pos = (double)(j + 0.5 - GhostCell) / ny * (*Ly);
            
            // Initialize The Flow Field Based On The Shock Wave Position
            // y = shock_slope * (x - shock_start_x)
            if (y_pos < shock_slope * (x_pos - shock_start_x)) {
                // Post Shock Region
                x[0][i][j] = pre_shock_state[0];   
                x[1][i][j] = pre_shock_state[1];   
                x[2][i][j] = pre_shock_state[2];
                x[3][i][j] = pre_shock_state[3];   
                
            } else {
                x[0][i][j] = post_shock_state[0];  
                x[1][i][j] = post_shock_state[1];  
                x[2][i][j] = post_shock_state[2];
                x[3][i][j] = post_shock_state[3];  
               
            }
        }
    }

    printf("Computational domain: [0, %f] x [0, %f]\n", *Lx, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Shock Mach number: 10\n");
    printf("Shock angle: 60 degrees\n");
    printf("Pre-shock state: (ρ, u, v, p) = (%f, %f, %f, %f)\n", 
           pre_shock_state[0], pre_shock_state[1], pre_shock_state[2], pre_shock_state[3]);
    printf("Post-shock state: (ρ, u, v, p) = (%f, %f, %f, %f)\n", 
           post_shock_state[0], post_shock_state[1], post_shock_state[2], post_shock_state[3]);
    printf("Reflecting wall from x = %f to x = %f\n", 1.0/6.0, 4.0);
}



static inline void initEulerpri2D_BlastWave(int var, int rows, int cols, 
                                           double (*x)[rows][cols], 
                                           double *Lx, double *Ly, double *Time) {
    printf("=== 2D Blast Wave Test Case ===\n");
    printf("Strong spherical explosion with initial pressure discontinuity\n");
    
    // Set Boundary Conditions
    bc_config.left = BC_REFLECTION;    
    bc_config.right = BC_REFLECTION;   
    bc_config.bottom = BC_REFLECTION;  
    bc_config.top = BC_REFLECTION;     
    
    // Set Computational Domain
    *Lx = 1.0;                                  // [-0.5, 0.5]
    *Ly = 1.0;                                  // [-0.5, 0.5]
    *Time = 0.01;                              // Typical Time In The Literature
    
    // High Pressure Zone Central Explosion
    double high_pressure = 10.0;   
    double high_density = 1.0;     
    
    // Low Pressure Area External Environment
    double low_pressure = 0.1;     
    double low_density = 1.0; 
    
    // Initial Speed Is Zero
    double init_velocity_x = 0.0;
    double init_velocity_y = 0.0;
    
    // Blast Radius
    double blast_radius = 0.1;
    
    //Grid Parameters
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    double dx = *Lx / nx;
    double dy = *Ly / ny;
    
    // Calculate The Center Coordinates Of The Domain For Spherical Symmetry
    double center_x = 0.0;
    double center_y = 0.0;
    
    int high_pressure_cells = 0;
    int low_pressure_cells = 0;
    
    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // Compute Physical Coordinates From 0 5 To 0 5
            double x_pos = (i + 0.5 - GhostCell) * dx - 0.5;
            double y_pos = (j + 0.5 - GhostCell) * dy - 0.5;
            
            // Distance To The Center
            double dx_center = x_pos - center_x;
            double dy_center = y_pos - center_y;
            double radius = sqrt(dx_center * dx_center + dy_center * dy_center);
            
            if (radius <= blast_radius) {
                // High Pressure Blast Area
                x[0][i][j] = high_density;
                x[1][i][j] = init_velocity_x;
                x[2][i][j] = init_velocity_y;
                x[3][i][j] = high_pressure;
                high_pressure_cells++;
            } else {
                // 低压环境区域
                x[0][i][j] = low_density;
                x[1][i][j] = init_velocity_x;
                x[2][i][j] = init_velocity_y;
                x[3][i][j] = low_pressure;
                low_pressure_cells++;
            }
        }
    }
    
    printf("Computational domain: [%f, %f] x [%f, %f]\n", -0.5, 0.5, -0.5, 0.5);
    printf("Final time: t = %f\n", *Time);
    printf("High pressure region (r <= %f):\n", blast_radius);
    printf("  Density: ρ = %f, Pressure: p = %f\n", high_density, high_pressure);
    printf("Low pressure region (r > %f):\n", blast_radius);
    printf("  Density: ρ = %f, Pressure: p = %f\n", low_density, low_pressure);
    printf("Initial velocity: (u, v) = (%f, %f)\n", init_velocity_x, init_velocity_y);
    printf("Gamma: γ = %f\n", M_gamma);
    printf("Grid cells in high pressure region: %d\n", high_pressure_cells);
    printf("Grid cells in low pressure region: %d\n", low_pressure_cells);
    printf("Pressure ratio: %f\n", high_pressure / low_pressure);
    printf("Expected features:\n");
    printf("  1. Strong outward propagating shock wave\n");
    printf("  2. Contact discontinuity separating fluids\n");
    printf("  3. Rarefaction wave moving inward\n");
    printf("  4. Complex shock interactions at late times\n");
}


static inline void initEulerpri2D_NohProblem(int var, int rows, int cols,
                                            double (*x)[rows][cols],
                                            double *Lx, double *Ly, double *Time) {
    printf("=== 2D Noh Problem Test Case ===\n");
    printf("Strong converging shock wave with infinite strength\n");
    
    // Set Boundary Conditions Approximation For Spherically Symmetric Problems
    bc_config.left = BC_REFLECTION;  
    bc_config.right = BC_REFLECTION; 
    bc_config.bottom = BC_REFLECTION;
    bc_config.top = BC_REFLECTION;   
    
    // Calculate The Domain Size First Quadrant Using Symmetry
    *Lx = 1.0;                                                          // [0, 1]
    *Ly = 1.0;                                                          // [0, 1]
    *Time = 0.6;                                                        // Typical Computation Time
    
    // Physical Parameters The Noh Problem Commonly Uses 5/3
    // Initial Conditions Uniform Medium Collapsing Inward
    double init_density = 1.0;
    double init_pressure = 1.0e-6; 
    
    // Radial Velocity Inward Magnitude 1 Pointing Toward The Origin
    double init_speed = 1.0; 
    
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    double dx = *Lx / nx;
    double dy = *Ly / ny;
    
    //Origin Coordinates Bottom Left
    double origin_x = 0.0;
    double origin_y = 0.0;
    
    // Add Small Perturbations To Avoid Complete Symmetry Optional Helps With Numerical Stability
    double perturbation_amplitude = 0.001;
    

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // Computational Physical Coordinates
            double x_pos = (i + 0.5 - GhostCell) * dx;
            double y_pos = (j + 0.5 - GhostCell) * dy;
            
            // Calculate The Distance To The Origin
            double distance = sqrt(x_pos * x_pos + y_pos * y_pos);
            
            // Calculate The Radial Unit Vector Pointing Towards The Origin
            double u_comp, v_comp;
            if (distance > 1.0e-10) {
                u_comp = -x_pos / distance;  
                v_comp = -y_pos / distance; 
            } else {
                u_comp = 0.0;
                v_comp = 0.0;
            }
            
            // Add Tiny Random Perturbations To Avoid Complete Symmetry
            double perturbation = perturbation_amplitude * 
                                 (2.0 * ((double)rand()/RAND_MAX) - 1.0);
            u_comp += perturbation;
            v_comp += perturbation;
            
            // Normalization Maintains Speed Magnitude
            double norm = sqrt(u_comp * u_comp + v_comp * v_comp);
            if (norm > 1.0e-10) {
                u_comp = u_comp * init_speed / norm;
                v_comp = v_comp * init_speed / norm;
            }
            
            // Set Conserved Variables
            x[0][i][j] = init_density;
            x[1][i][j] = init_density * u_comp;
            x[2][i][j] = init_density * v_comp;
            x[3][i][j] = init_pressure/(M_gamma-1.0) + 
                        0.5 * init_density * (u_comp * u_comp + v_comp * v_comp);
            
            // Add Small Density Perturbations Near The Center Region Optional
            if (distance < 0.05) {
                x[0][i][j] = init_density * (1.0 + perturbation_amplitude * 
                                            ((double)rand()/RAND_MAX - 0.5));
            }
        }
    }
    
    // Theoretical Solution Parameters Of Noh Problem
    double shock_speed = (M_gamma - 1.0) / 2.0;  
    double density_jump = (M_gamma + 1.0) / (M_gamma - 1.0); 
    double post_shock_density = init_density * density_jump;
    double post_shock_pressure = 0.5 * init_density * init_speed * init_speed * 
                                (M_gamma + 1.0);
    
    printf("Computational domain: [%f, %f] x [%f, %f]\n", 0.0, *Lx, 0.0, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Initial conditions:\n");
    printf("  Density: ρ₀ = %f\n", init_density);
    printf("  Pressure: p₀ = %e (approximately zero)\n", init_pressure);
    printf("  Radial velocity: v_r = -%f (inward)\n", init_speed);
    printf("  Gamma: γ = %f (5/3)\n", M_gamma);
    printf("\nTheoretical post-shock state (self-similar solution):\n");
    printf("  Shock speed: v_s = %f\n", shock_speed);
    printf("  Shock position at t=%f: r_s = %f\n", *Time, shock_speed * (*Time));
    printf("  Density jump ratio: %f\n", density_jump);
    printf("  Post-shock density: ρ₁ = %f\n", post_shock_density);
    printf("  Post-shock pressure: p₁ = %f\n", post_shock_pressure);
    printf("  Post-shock velocity: v₁ = 0 (stagnation)\n");
    printf("\nNumerical challenges:\n");
    printf("  1. Infinite strength shock (p₀ ≈ 0)\n");
    printf("  2. Strong convergence to center singularity\n");
    printf("  3. Preservation of spherical symmetry\n");
    printf("  4. Shock overheating problem verification\n");
    
    // 计算激波最终位置
    double shock_final_radius = shock_speed * (*Time);
    printf("\nExpected shock front at t=%f: radius = %f\n", *Time, shock_final_radius);
    
    if (shock_final_radius > 0.5) {
        printf("Warning: Shock may reach boundary before final time!\n");
        printf("Consider reducing final time or increasing domain size.\n");
    }
}


// Case: Odd-Even Decoupling Test (Based on Quirk 1994)
static inline void initEulerpri2D_OddEvenDecoupling(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Odd-Even Decoupling Test (Quirk 1994) ===\n");

    // Boundary conditions: Outflow (simulating a long duct)
    bc_config.left = BC_OUTFLOW;      
    bc_config.right = BC_OUTFLOW;     
    bc_config.bottom = BC_REFLECTION;    
    bc_config.top = BC_REFLECTION;       

    int i, j;
    *Lx = 20.0;    // Pipe Length Long Enough To Observe Shock Wave Propagation
    *Ly = 1.0;     // Pipe Width

    // Shock parameters: Strong shock M=6, γ=1.4
    double gamma = 1.4;
    double Ms = 6.0;  // Shock Mach number
    
    // Pre-shock state (right side of shock)
    double rho_pre = 1.4;
    double u_pre = 0.0;    // Fluid at rest before shock
    double v_pre = 0.0;
    double p_pre = 1.0;
    double a_pre = sqrt(gamma * p_pre / rho_pre);  // Sound speed
    
    // Calculate post-shock state (left side, using Rankine-Hugoniot relations)
    double rho_post = rho_pre * ((gamma + 1.0) * Ms * Ms) / ((gamma - 1.0) * Ms * Ms + 2.0);
    double p_post = p_pre * (2.0 * gamma * Ms * Ms - (gamma - 1.0)) / (gamma + 1.0);
    double u_post = Ms * a_pre * (1.0 - rho_pre / rho_post);  // Post-shock velocity (positive to right)
    
    printf("Shock parameters: M=%.2f, γ=%.2f\n", Ms, gamma);
    printf("Pre-shock state: ρ=%.4f, u=%.4f, p=%.4f\n", rho_pre, u_pre, p_pre);
    printf("Post-shock state: ρ=%.4f, u=%.4f, p=%.4f\n", rho_post, u_post, p_post);
    
    // Odd-even perturbation amplitude (very small, simulating numerical error or grid imperfection)
    double perturbation = 1e-6;
    
    // Physical grid dimensions (excluding ghost cells)
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;

    #pragma omp parallel for collapse(2)
    // Initialize only physical cells (ghost cells handled by boundary conditions)
    for (i = GhostCell; i < rows - GhostCell; i++) {
        for (j = GhostCell; j < cols - GhostCell; j++) {
            
            // Calculate physical coordinates (cell-centered coordinates)
            // Using i+0.5 to get cell center in x-direction
            double x_pos = (double)(i + 0.5 - GhostCell) / nx * (*Lx);
            double y_pos = (double)(j + 0.5 - GhostCell) / ny * (*Ly);
            
            // Shock location: initially positioned at x = Lx/4
            if (x_pos < (*Lx) * 0.25) {
                // Post-shock region (left side of shock)
                x[0][i][j] = rho_post;  // Density
                x[1][i][j] = u_post;    // x-velocity
                x[2][i][j] = v_pre;     // y-velocity (zero)
                x[3][i][j] = p_post;    // Pressure
            } else {
                // Pre-shock region (right side of shock)
                x[0][i][j] = rho_pre;   // Density
                x[1][i][j] = u_pre;     // x-velocity
                x[2][i][j] = v_pre;     // y-velocity (zero)
                x[3][i][j] = p_pre;     // Pressure
            }
            
            // KEY STEP: Apply odd-even perturbation to density and pressure in post-shock region
            // This simulates the grid centerline perturbation described in Quirk's paper
            /*if (x_pos < (*Lx) * 1) {
                if (j % 2 == 0) {  // Even grid points in y-direction
                    x[0][i][j] += perturbation;  // Increase density
                    x[3][i][j] += perturbation;  // Increase pressure
                } else {           // Odd grid points in y-direction
                    x[0][i][j] -= perturbation;  // Decrease density
                    x[3][i][j] -= perturbation;  // Decrease pressure
                }
            }*/
            
            // Optional: Apply small perturbation in pre-shock region to simulate grid imperfection
            // Uncomment if you want to be closer to Quirk's original setup
            /*
            if (x_pos >= (*Lx) * 0.25) {
                if ((i + j) % 2 == 0) {  // Different parity pattern (diagonal)
                    x[0][i][j] += perturbation * 0.1;
                    x[3][i][j] += perturbation * 0.1;
                } else {
                    x[0][i][j] -= perturbation * 0.1;
                    x[3][i][j] -= perturbation * 0.1;
                }
            }
            */
        }
    }

    *Time = 2.0;  // Sufficient time for shock to propagate and decoupling to develop
    
    printf("Computational domain: Lx = %.2f, Ly = %.2f\n", *Lx, *Ly);
    printf("Final simulation time: t = %.2f\n", *Time);
    printf("Odd-Even perturbation amplitude: %e\n", perturbation);
    printf("Boundary Conditions: All outflow\n");
    printf("Initial shock location: x = %.2f\n", *Lx * 0.25);
    printf("Physical cells initialized: %d x %d\n", (int)nx, (int)ny);
    printf("Ghost cells: %d per boundary\n", GhostCell);
    printf("\nTest Configuration Notes:\n");
    printf("1. Designed to trigger Odd-Even Decoupling phenomenon with Roe-type solvers\n");
    printf("2. Strong shock (M=6) aligned with x-direction\n");
    printf("3. Small perturbations (±%e) applied to density and pressure in post-shock region\n", perturbation);
    printf("4. Perturbations applied with opposite signs on even/odd y-grid lines\n");
    printf("5. Expected: With Roe solver, decoupling should develop along shock front\n");
}


double L_fixed_value_state[4] = {1.4, (1.4)*3.0, 0.0, (1.0)/0.4 + 0.5 * (1.4) *3.0 *3.0}; 
double R_fixed_value_state[4] = {1.4, 1.4*3.0, 0.0, 1.0/0.4 + 0.5 * 1.4 *3.0 *3.0}; 
double B_fixed_value_state[4] = {1.4, 1.4*3.0, 0.0, 1.0/0.4 + 0.5 * 1.4 *3.0 *3.0}; 
double T_fixed_value_state[4] = {1.4, 1.4*3.0, 0.0, 1.0/0.4 + 0.5 * 1.4 *3.0 *3.0};
static inline void initEulerpri2D_BackwardStep(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Woodward & Colella Mach 3.0 flow over a forward facing step (1984) ===\n");

    // Complete computational domain boundary conditions
    bc_config.left = BC_FIXED_VALUE;
    bc_config.right = BC_OUTFLOW;
    bc_config.bottom = BC_REFLECTION;
    bc_config.top = BC_REFLECTION;

    // Geometric parameters of the domain [0.0, 3.0] × [0.0, 1.0]
    *Lx = 3.0;
    *Ly = 1.0;

    // Final simulation time: t = 4.0
    *Time = 4.0;

    // Step geometric parameters
    double step_height = 0.2;
    double step_position = 0.6;

    // Inflow conditions with Mach number 3
    double inflow_mach = 3.0;
    double inflow_density = 1.4;
    double inflow_pressure = 1.0;

    // Calculate derived inflow quantities
    double sound_speed = sqrt(M_gamma * inflow_pressure / inflow_density);
    double inflow_velocity_x = inflow_mach * sound_speed;
    double inflow_velocity_y = 0.0;
    double inflow_energy = inflow_pressure / (M_gamma - 1.0) +
                          0.5 * inflow_density * (inflow_velocity_x * inflow_velocity_x);

    // double inflow_state[4] = {inflow_density, inflow_density * inflow_velocity_x, inflow_velocity_y, inflow_energy}; // Default values

    // Print setup summary
    printf("Woodward & Colella (1984) setup:\n");
    printf("Computational domain: [0, %f] x [0, %f]\n", *Lx, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Step configuration:\n");
    printf("  - Step height: %f (20%% of channel height)\n", step_height);
    printf("  - Step position: x = %f\n", step_position);
    printf("Inflow conditions (Mach 3):\n");
    printf("  - Mach number: %f\n", inflow_mach);
    printf("  - Density: %f\n", inflow_density);
    printf("  - Pressure: %f\n", inflow_pressure);
    printf("  - Velocity: (%f, %f)\n", inflow_velocity_x, inflow_velocity_y);
    printf("  - Sound speed: %f\n", sound_speed);
    printf("  - Gamma: %f\n", M_gamma);

    // Compute grid parameters
    int nx = rows - 2 * GhostCell;      // Number of interior cells in x-direction
    int ny = cols - 2 * GhostCell;      // Number of interior cells in y-direction
    double dx = *Lx / nx;
    double dy = *Ly / ny;

    // Initialize entire domain with inflow state
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            x[0][i][j] = inflow_density;      // Density
            x[1][i][j] = inflow_velocity_x;   // x-velocity
            x[2][i][j] = inflow_velocity_y;   // y-velocity
            x[3][i][j] = inflow_pressure;     // Pressure
        }
    }

    // Verify physical plausibility of initial conditions
    printf("\nInitial condition verification:\n");
    printf("Total cells in solid region (step): approximately %d\n",
           (int)(step_position / dx * step_height / dy));
    printf("Total cells in fluid region: approximately %d\n",
           nx * ny - (int)(step_position / dx * step_height / dy));

    // Compute and display key parameters
    double reynolds_number = 0.0;  // Inviscid Euler equations → infinite Reynolds number
    printf("Reynolds number: infinite (inviscid Euler equations)\n");
    printf("Flow features expected (from Woodward & Colella 1984):\n");
    printf("  1. Strong shock from step corner\n");
    printf("  2. Expansion fan at step corner\n");
    printf("  3. Recompression shock downstream\n");
    printf("  4. Separation region behind step\n");
    printf("  5. Shear layer development\n");
}




// Map Primitive Variables To Conserved Variables
static inline void initEulerconser2D(int var, int rows, int cols, double (*x)[rows][cols], double (*y)[rows][cols]) {
    int j,k;
    #pragma omp parallel for collapse(2)
    for (j = 0; j < rows; j++) {
        for ( k = 0; k < cols; k++){
            double rho, u, v, p;
            rho = x[0][j][k];
            u = x[1][j][k];
            v = x[2][j][k];
            p = x[3][j][k];
            y[0][j][k] = rho;
            y[1][j][k] = u * rho;
            y[2][j][k] = v * rho;
            y[3][j][k] = p /(M_gamma-1) + 0.5 * rho * (pow(u,2) + pow(v,2));
        }
    }
}


// Map Conserved Variables To Flux Functions
static inline void initEulerflux2D(int var, int rows, int cols, double (*y)[rows][cols], double (*f)[rows][cols],double (*g)[rows][cols]) {
    int j,k;

    #pragma omp parallel for collapse(2)

    for (j = 0; j < rows; j++) {
        for ( k = 0; k < cols; k++){
            double rho, u, v, p;
            rho = y[0][j][k];
            u = y[1][j][k]/y[0][j][k];
            v = y[2][j][k]/y[0][j][k];
            p = (y[3][j][k] - 0.5 * rho * (pow(u,2) + pow(v,2)))*(M_gamma-1);
            //x dir initilize
            f[0][j][k] = rho * u;
            f[1][j][k] = rho * u * u + p;
            f[2][j][k] = rho * u * v;
            f[3][j][k] = (y[3][j][k] + p) * u;
            //y dir initilize
            g[0][j][k] = rho * v;
            g[1][j][k] = rho * u * v;
            g[2][j][k] = rho * v * v + p;
            g[3][j][k] = (y[3][j][k] + p) * v;

        }
    }
}

// Map Conserved Variables To Primitive Variables
static inline void Con_to_Pri_2D(int var, int rows, int cols, double (*x)[rows][cols], double (*y)[rows][cols]){

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell-1; i <= rows-GhostCell; i++) {
        for (int j = GhostCell-1; j <= cols-GhostCell; j++){
            double rho, u, v, p;
            rho = y[0][i][j];
            u = y[1][i][j]/y[0][i][j];
            v = y[2][i][j]/y[0][i][j];
            p = (y[3][i][j] - 0.5 * rho * (pow(u,2) + pow(v,2)))*(M_gamma-1);
            //initilize
            x[0][i][j] = rho;
            x[1][i][j] = u;
            x[2][i][j] = v;
            x[3][i][j] = p;
        }
    }
}


// Initialization Function For Unified 2 D Euler Equation Test Cases
static inline void initEulerTestCase(TestCase2D test_case, 
                                    int var, int rows, int cols,
                                    double (*x)[rows][cols],
                                    double *Lx, double *Ly, double *Time) {
    
    // Select The Corresponding Initialization Function According To The Test Cases
    switch(test_case) {
        case TEST_PRECISION:
            initEulerpri1D_PrecisionTest(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_1D_SHOCKTUBE:
            initEulerpri1D_Shocktube(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_1D_CONTACTWAVE:
            initEulerpri1D_ContactWave(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_1D_SHOCKIMPACT:
            initEulerpri1D_ShockImpact(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_1D_IMPACTWALL:
            initEulerpri1D_ShockImpactWall(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_1D_DOUBLERARE:
            initEulerpri1D_DoubleRare(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_1D_NOHPROBLEM:
            initEulerpri1D_NohProblem(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE1:
            initEulerpri2D_Shocktube1(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE2:
            initEulerpri2D_Shocktube2(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE3:
            initEulerpri2D_Shocktube3(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE4:
            initEulerpri2D_Shocktube4(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE5:
            initEulerpri2D_Shocktube5(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_TAYLOR_GREEN_VORTEX:
            initEulerpri2D_TaylorGreenVortex(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_GAUSSIAN_PULSE:
            initEulerpri2D_GaussianPulse(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_KELVIN_HELMHOLTZ:
            initEulerpri2D_KelvinHelmholtz(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_RAYLEIGH_TAYLOR:
            initEulerpri2D_RayleighTaylor(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_DOUBLE_MACH_REFLECTION:
            initEulerpri2D_DoubleMachReflection(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_BLAST_WAVE:
            initEulerpri2D_BlastWave(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_NOH_PROBLEM:
            initEulerpri2D_NohProblem(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_BACKWARD_STEP:
            initEulerpri2D_BackwardStep(var, rows, cols, x, Lx, Ly, Time);
            break;

        case TEST_OddEven_Decoupling :
            initEulerpri2D_OddEvenDecoupling(var, rows, cols, x, Lx, Ly, Time);
            break;

        default:
            printf("Error: Unknown test case selected!\n");
            // Default Use Of 1 D Shock Tube
            initEulerpri1D_Shocktube(var, rows, cols, x, Lx, Ly, Time);
            break;
    }
}


static inline void Init_Euler_2D(TestCase2D test_case,
                                int var, int rows, int cols, int GC, 
                                double (*x)[rows][cols], double (*y)[rows][cols], 
                                double (*f)[rows][cols], double (*g)[rows][cols], 
                                double *Lx, double *Ly, double *Time) {
    
    printf("========== Initializing 2D Euler Equations ==========\n");
    
    // Initialize Array Set To Zero
    initEuler2D(var, rows, cols, x);
    initEuler2D(var, rows, cols, y);
    initEuler2D(var, rows, cols, f);
    initEuler2D(var, rows, cols, g);
    
    // Initialize Using A Unified Test Case Function
    printf("Selected test case: ");
    switch(test_case) {
        case TEST_PRECISION: printf("Scheme Precision Test\n"); break;
        case TEST_1D_SHOCKTUBE: printf("1D Sod Shock Tube\n"); break;
       case TEST_1D_CONTACTWAVE: printf("1D Contact Discontinuity Test\n"); break;
        case TEST_1D_SHOCKIMPACT: printf("1D Shock Impact Problem\n"); break;
        case TEST_1D_IMPACTWALL: printf("1D Shock-Wall Impact Problem\n"); break;
        case TEST_1D_DOUBLERARE: printf("1D Double Rarefaction Wave Test\n"); break;
        case TEST_1D_NOHPROBLEM: printf("1D Noh Problem Test\n"); break;
        case TEST_2D_SHOCKTUBE_CASE1: printf("2D Riemann Problem Case 1\n"); break;
        case TEST_2D_SHOCKTUBE_CASE2: printf("2D Riemann Problem Case 2\n"); break;
        case TEST_2D_SHOCKTUBE_CASE3: printf("2D Riemann Problem Case 3\n"); break;
        case TEST_2D_SHOCKTUBE_CASE4: printf("2D Riemann Problem Case 4\n"); break;
        case TEST_2D_SHOCKTUBE_CASE5: printf("2D Riemann Problem Case 5\n"); break;
        case TEST_TAYLOR_GREEN_VORTEX: printf("Taylor-Green Vortex\n"); break;
        case TEST_GAUSSIAN_PULSE: printf("Gaussian Density Pulse\n"); break;
        case TEST_KELVIN_HELMHOLTZ: printf("Kelvin-Helmholtz Instability (Smooth)\n"); break;
        case TEST_RAYLEIGH_TAYLOR: printf("Rayleigh-Taylor Instability\n"); break;
        case TEST_DOUBLE_MACH_REFLECTION: printf("Double Mach Reflection\n"); break;
        case TEST_BLAST_WAVE: printf("Blast Wave (Spherical Explosion)\n"); break;
        case TEST_NOH_PROBLEM: printf("Noh Problem (Converging Shock)\n"); break;
        case TEST_BACKWARD_STEP: printf("Backward Step Flow\n"); break;
        case TEST_OddEven_Decoupling: printf("Odd Even Decoupling\n"); break;
    }
    
    // Call A Unified Initialization Function
    initEulerTestCase(test_case, var, rows, cols, x, Lx, Ly, Time);
    
    initEulerconser2D(var, rows, cols, x, y);
    
    initEulerflux2D(var, rows, cols, y, f, g);
    
    printf("========== Initialization Complete ==========\n\n");
}


static inline void Init_Precision(TestCase2D test_case,
                                int var, int rows, int cols, int GC, 
                                double (*x)[rows][cols], double (*y)[rows][cols], 
                                double (*f)[rows][cols], double (*g)[rows][cols], 
                                double *Lx, double *Ly, double *Time) {
    

    // Initialize Array Set To Zero
    initEuler2D(var, rows, cols, x);
    initEuler2D(var, rows, cols, y);
    initEuler2D(var, rows, cols, f);
    initEuler2D(var, rows, cols, g);
    
    //
    initEulerpri1D_PrecisionTest(var, rows, cols, x, Lx, Ly, Time);
    
    initEulerconser2D(var, rows, cols, x, y);
    
    initEulerflux2D(var, rows, cols, y, f, g);
    
}

/*                                            *********                                          */
/*                                        Test Case Name Functions                               */
/*                                            *********                                          */

// Get full test case name (for display purposes)
const char* getTestCaseName(TestCase2D test_case) {
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: return "1D Sod Shock Tube Problem";
        case TEST_1D_CONTACTWAVE: return "1D Contact Discontinuity Test";
        case TEST_1D_SHOCKIMPACT: return "1D Symmetric Shock Collision";
        case TEST_1D_IMPACTWALL: return "1D Shock Reflection from Rigid Wall";
        case TEST_1D_DOUBLERARE: return "1D Double Rarefaction Wave Problem";
        case TEST_1D_NOHPROBLEM: return "1D Noh Implosion Problem";
        case TEST_2D_SHOCKTUBE_CASE1: return "2D Riemann Problem - Case 1";
        case TEST_2D_SHOCKTUBE_CASE2: return "2D Riemann Problem - Case 2";
        case TEST_2D_SHOCKTUBE_CASE3: return "2D Riemann Problem - Case 3";
        case TEST_2D_SHOCKTUBE_CASE4: return "2D Riemann Problem - Case 4";
        case TEST_2D_SHOCKTUBE_CASE5: return "2D Riemann Problem - Case 5";
        case TEST_TAYLOR_GREEN_VORTEX: return "Taylor-Green Vortex";
        case TEST_GAUSSIAN_PULSE: return "Gaussian Pulse";
        case TEST_KELVIN_HELMHOLTZ: return "Kelvin-Helmholtz Instability";
        case TEST_RAYLEIGH_TAYLOR: return "Rayleigh-Taylor Instability";
        case TEST_DOUBLE_MACH_REFLECTION: return "Double Mach Reflection";
        case TEST_BACKWARD_STEP: return "Backward-Facing Step Flow";
        case TEST_OddEven_Decoupling: return "Odd-Even Decoupling Test";
        case TEST_BLAST_WAVE: return "Spherical Blast Wave";
        case TEST_NOH_PROBLEM: return "2D Noh Problem";
        default: return "Unknown Test Case";
    }
}

// Get short test case name (for file naming)
const char* getTestCaseShortName(TestCase2D test_case) {
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: return "SOD";
        case TEST_1D_CONTACTWAVE: return "CONTACT";
        case TEST_1D_SHOCKIMPACT: return "SHOCKIMPACT";
        case TEST_1D_IMPACTWALL: return "IMPACTWALL";
        case TEST_1D_DOUBLERARE: return "DOUBLERARE";
        case TEST_1D_NOHPROBLEM: return "NOH1D";
        case TEST_2D_SHOCKTUBE_CASE1: return "RIEMANN1";
        case TEST_2D_SHOCKTUBE_CASE2: return "RIEMANN2";
        case TEST_2D_SHOCKTUBE_CASE3: return "RIEMANN3";
        case TEST_2D_SHOCKTUBE_CASE4: return "RIEMANN4";
        case TEST_2D_SHOCKTUBE_CASE5: return "RIEMANN5";
        case TEST_TAYLOR_GREEN_VORTEX: return "TGV";
        case TEST_GAUSSIAN_PULSE: return "GAUSSIAN";
        case TEST_KELVIN_HELMHOLTZ: return "KH";
        case TEST_RAYLEIGH_TAYLOR: return "RT";
        case TEST_DOUBLE_MACH_REFLECTION: return "DMR";
        case TEST_BACKWARD_STEP: return "BACKSTEP";
        case TEST_OddEven_Decoupling: return "OEDC";
        case TEST_BLAST_WAVE: return "BLAST";
        case TEST_NOH_PROBLEM: return "NOH2D";
        default: return "UNKNOWN";
    }
}


// Get The Description Information Of The Test Case
static inline const char* getTestCaseDescription(TestCase2D test_case) {
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: return "1D Sod shock tube problem - classical benchmark";
        case TEST_1D_CONTACTWAVE: return "1D Contact discontinuity problem - for testing contact resolution";
        case TEST_1D_SHOCKIMPACT: return "1D Symmetric shock collision problem - tests overheating behavior";
        case TEST_1D_IMPACTWALL: return "1D Shock reflection from rigid wall - wall heating test";
        case TEST_1D_DOUBLERARE: return "1D Double rarefaction wave problem - low density test";
        case TEST_1D_NOHPROBLEM: return "1D Noh problem - strong shock implosion test";  
        case TEST_2D_SHOCKTUBE_CASE1: return "2D Riemann problem case 1 - four interacting states";
        case TEST_2D_SHOCKTUBE_CASE2: return "2D Riemann problem case 2 - complex wave interactions";
        case TEST_2D_SHOCKTUBE_CASE3: return "2D Riemann problem case 3 - shock interactions";
        case TEST_2D_SHOCKTUBE_CASE4: return "2D Riemann problem case 4 - symmetric 4-quadrant";
        case TEST_2D_SHOCKTUBE_CASE5: return "2D Riemann problem case 5 - shock-rarefaction";
        case TEST_TAYLOR_GREEN_VORTEX: return "Taylor-Green vortex - periodic decaying turbulence";
        case TEST_GAUSSIAN_PULSE: return "Gaussian density pulse in uniform flow";
        case TEST_KELVIN_HELMHOLTZ: return "Kelvin-Helmholtz instability - shear layer with smooth interface";
        case TEST_RAYLEIGH_TAYLOR: return "Rayleigh-Taylor instability - heavy fluid over light fluid";
        case TEST_DOUBLE_MACH_REFLECTION: return "Double Mach reflection - strong shock interaction";
        case TEST_BLAST_WAVE: return "Blast wave - strong spherical explosion";
        case TEST_NOH_PROBLEM: return "Noh problem - infinite strength converging shock";
        case TEST_BACKWARD_STEP: return "Backward step flow - supersonic flow with separation";
        case TEST_OddEven_Decoupling: return "Odd Even Decoupling  - Test Riemann Problem";
        default: return "Unknown test case";
    }
}


/*static inline void Init_BACKWARD_STEP(int var, int rows, int cols, int GC, 
                                double (*x)[rows][cols], double (*y)[rows][cols], 
                                double (*f)[rows][cols], double (*g)[rows][cols], 
                                double *Lx, double *Ly, double *Time) {
    
    printf("========== Initializing 2D Euler Equations ==========\n");
    
    // Initialize Array Set To Zero
    initEuler2D(var, rows, cols, x);
    initEuler2D(var, rows, cols, y);
    initEuler2D(var, rows, cols, f);
    initEuler2D(var, rows, cols, g);


    for (int i = 0; i < rows; i++) {
        for (int  j = 0; j < cols; j++){
            x[0][i][j] = 1.4;
            x[1][i][j] = 3.0;
            x[2][i][j] = 0.0;
            x[3][i][j] = 1.0;
        }
    }
    
    // 将原始变量转换为守恒变量
    initEulerconser2D(var, rows, cols, x, y);
    
    // 计算通量
    initEulerflux2D(var, rows, cols, y, f, g);
    
    printf("========== Initialization Complete ==========\n\n");
}*/

#endif