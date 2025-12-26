#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mesh.h"
#include "initialize.h"
#include "scheme.h"
#include "function.h"

// Define Boundary Type Macros In The Header File
#define L "left"
#define R "right"  
#define T "top"
#define B "bottom"

//
typedef enum {
    BC_OUTFLOW = 0,     
    BC_REFLECTION = 1,  
    BC_PERIODICITY = 2,
    BC_INFLOW = 3,     
    BC_FIXED_VALUE = 4  
} BC_Type;

// Boundary Configuration Structure
typedef struct {
    BC_Type left;
    BC_Type right;
    BC_Type bottom;
    BC_Type top;
} BoundaryConfig;

// Global Border Configuration Declaration
extern BoundaryConfig bc_config;

// Global Entry Status Declaration In A Certain Header File
extern double inflow_state[4]; // [密度, x动量, y动量, 压力]

// Global Constant State Declaration
extern double L_fixed_value_state[4]; 
extern double R_fixed_value_state[4]; 
extern double B_fixed_value_state[4]; 
extern double T_fixed_value_state[4]; 

double inflow_state[4] = {1.4, 1.4*3.0, 0.0, 1.0/0.4 + 0.5 * 1.4 *3.0 *3.0}; // 默认值


/**
 * Apply fixed-value boundary conditions for 2D Euler equations
 * This function sets ghost cell values to specified fixed states for each boundary
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param x Solution array [var][rows][cols]
 * @param boundary_type String identifier for boundary: "L" (left), "R" (right), "B" (bottom), "T" (top)
 * @param Ghost_cell Number of ghost cells on each side
 */
static inline void BC_FixedValue_2D(int var, int rows, int cols, double (*x)[rows][cols], 
                                    const char* boundary_type, int Ghost_cell) {
    double *fixed_state = NULL;
    
    // Select fixed state array based on boundary type
    if (strcmp(boundary_type, L) == 0) {
        fixed_state = L_fixed_value_state;      // Left boundary fixed state
    }
    else if (strcmp(boundary_type, R) == 0) {
        fixed_state = R_fixed_value_state;      // Right boundary fixed state
    }
    else if (strcmp(boundary_type, B) == 0) {
        fixed_state = B_fixed_value_state;      // Bottom boundary fixed state
    }
    else if (strcmp(boundary_type, T) == 0) {
        fixed_state = T_fixed_value_state;      // Top boundary fixed state
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
        return;
    }
    
    // Apply fixed-value boundary conditions to ghost cells
    if (strcmp(boundary_type, L) == 0) {
        // Left boundary: set ghost cells to fixed state
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                // Set all ghost cells at left boundary to the same fixed value
                x[i][Ghost_cell-1][k] = fixed_state[i];  // First ghost cell
                x[i][Ghost_cell-2][k] = fixed_state[i];  // Second ghost cell
                x[i][Ghost_cell-3][k] = fixed_state[i];  // Third ghost cell
                x[i][Ghost_cell-4][k] = fixed_state[i];  // Fourth ghost cell
            }
        }
    }
    else if (strcmp(boundary_type, R) == 0) {
        // Right boundary: set ghost cells to fixed state
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][rows - Ghost_cell][k] = fixed_state[i];      // First ghost cell
                x[i][rows - Ghost_cell + 1][k] = fixed_state[i];  // Second ghost cell
                x[i][rows - Ghost_cell + 2][k] = fixed_state[i];  // Third ghost cell
                x[i][rows - Ghost_cell + 3][k] = fixed_state[i];  // Fourth ghost cell
            }
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // Bottom boundary: set ghost cells to fixed state
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][Ghost_cell-1] = fixed_state[i];      // First ghost cell
                x[i][j][Ghost_cell-2] = fixed_state[i];      // Second ghost cell
                x[i][j][Ghost_cell-3] = fixed_state[i];      // Third ghost cell
                x[i][j][Ghost_cell-4] = fixed_state[i];      // Fourth ghost cell
            }
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // Top boundary: set ghost cells to fixed state
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][cols - Ghost_cell] = fixed_state[i];      // First ghost cell
                x[i][j][cols - Ghost_cell + 1] = fixed_state[i];  // Second ghost cell
                x[i][j][cols - Ghost_cell + 2] = fixed_state[i];  // Third ghost cell
                x[i][j][cols - Ghost_cell + 3] = fixed_state[i];  // Fourth ghost cell
            }
        }
    }
}

/**
 * Apply inflow boundary conditions for 2D Euler equations
 * This function sets ghost cell values to specified inflow state for each boundary
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param x Solution array [var][rows][cols]
 * @param boundary_type String identifier for boundary: "L" (left), "R" (right), "B" (bottom), "T" (top)
 * @param Ghost_cell Number of ghost cells on each side
 */
static inline void BC_Inflow_2D(int var, int rows, int cols, double (*x)[rows][cols], 
                                const char* boundary_type, int Ghost_cell) {
    // Use global inflow_state array as inflow condition
    if (strcmp(boundary_type, L) == 0) {
        // Left boundary inflow (e.g., supersonic inflow for backward step problem)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][Ghost_cell-1][k] = inflow_state[i];  // First ghost cell from left
                x[i][Ghost_cell-2][k] = inflow_state[i];  // Second ghost cell
                x[i][Ghost_cell-3][k] = inflow_state[i];  // Third ghost cell
                x[i][Ghost_cell-4][k] = inflow_state[i];  // Fourth ghost cell
            }
        }
    }
    else if (strcmp(boundary_type, R) == 0) {
        // Right boundary inflow (less common but included for completeness)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][rows - Ghost_cell][k] = inflow_state[i];      // First ghost cell from right
                x[i][rows - Ghost_cell + 1][k] = inflow_state[i];  // Second ghost cell
                x[i][rows - Ghost_cell + 2][k] = inflow_state[i];  // Third ghost cell
                x[i][rows - Ghost_cell + 3][k] = inflow_state[i];  // Fourth ghost cell
            }
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // Bottom boundary inflow (e.g., injection flow)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][Ghost_cell-1] = inflow_state[i];  // First ghost cell from bottom
                x[i][j][Ghost_cell-2] = inflow_state[i];  // Second ghost cell
                x[i][j][Ghost_cell-3] = inflow_state[i];  // Third ghost cell
                x[i][j][Ghost_cell-4] = inflow_state[i];  // Fourth ghost cell
            }
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // Top boundary inflow (e.g., free stream)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][cols - Ghost_cell] = inflow_state[i];      // First ghost cell from top
                x[i][j][cols - Ghost_cell + 1] = inflow_state[i];  // Second ghost cell
                x[i][j][cols - Ghost_cell + 2] = inflow_state[i];  // Third ghost cell
                x[i][j][cols - Ghost_cell + 3] = inflow_state[i];  // Fourth ghost cell
            }
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}


/**
 * Apply outflow boundary conditions for 2D Euler equations
 * Simple zero-gradient (Neumann) condition for outflow boundaries
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param x Solution array [var][rows][cols]
 * @param boundary_type String identifier for boundary: "L" (left), "R" (right), "B" (bottom), "T" (top)
 * @param Ghost_cell Number of ghost cells on each side
 */
static inline void BC_OutFlow_2D(int var, int rows, int cols, double (*x)[rows][cols], 
                                 const char* boundary_type, int Ghost_cell){
    // Apply boundary conditions based on boundary type string
    if (strcmp(boundary_type, L) == 0) {
        // Left boundary outflow condition (zero-gradient extrapolation)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k <= cols-Ghost_cell; k++){
                x[i][Ghost_cell-1][k] = x[i][Ghost_cell][k];  // Extrapolate from first interior cell
                x[i][Ghost_cell-2][k] = x[i][Ghost_cell][k];
                x[i][Ghost_cell-3][k] = x[i][Ghost_cell][k];
                x[i][Ghost_cell-4][k] = x[i][Ghost_cell][k];
            }
        }   
    }
    else if (strcmp(boundary_type, R) == 0) {
        // Right boundary non-reflective outflow condition
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k <= cols-Ghost_cell; k++){
                x[i][rows-Ghost_cell][k] = x[i][rows-Ghost_cell-1][k];  // Extrapolate from last interior cell
                x[i][rows-Ghost_cell+1][k] = x[i][rows-Ghost_cell-1][k];
                x[i][rows-Ghost_cell+2][k] = x[i][rows-Ghost_cell-1][k];
                x[i][rows-Ghost_cell+3][k] = x[i][rows-Ghost_cell-1][k];
            }                   
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // Bottom boundary outflow condition
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j <= rows-Ghost_cell; j++){
                x[i][j][Ghost_cell-1] = x[i][j][Ghost_cell];
                x[i][j][Ghost_cell-2] = x[i][j][Ghost_cell];
                x[i][j][Ghost_cell-3] = x[i][j][Ghost_cell];
                x[i][j][Ghost_cell-4] = x[i][j][Ghost_cell];
            }                   
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // Top boundary outflow condition
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j <= rows-Ghost_cell; j++){
                x[i][j][cols-Ghost_cell] = x[i][j][cols-Ghost_cell-1];
                x[i][j][cols-Ghost_cell+1] = x[i][j][cols-Ghost_cell-1];
                x[i][j][cols-Ghost_cell+2] = x[i][j][cols-Ghost_cell-1];
                x[i][j][cols-Ghost_cell+3] = x[i][j][cols-Ghost_cell-1];
            }                   
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}

/**
 * Apply periodic boundary conditions for 2D Euler equations
 * Ghost cells are filled from the opposite side of the domain
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param x Solution array [var][rows][cols]
 * @param boundary_type String identifier for boundary: "L" (left), "R" (right), "B" (bottom), "T" (top)
 * @param Ghost_cell Number of ghost cells on each side
 */
static inline void BC_Periodicity_2D(int var, int rows, int cols, double (*x)[rows][cols], 
                                     const char* boundary_type, int Ghost_cell){
    if (strcmp(boundary_type, L) == 0) {
        // Left boundary <- Right interior region
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){  
                x[i][Ghost_cell-1][k] = x[i][rows - Ghost_cell - 1][k];  // Left ghost <- Right interior
                x[i][Ghost_cell-2][k] = x[i][rows - Ghost_cell - 2][k];
                x[i][Ghost_cell-3][k] = x[i][rows - Ghost_cell - 3][k];
                x[i][Ghost_cell-4][k] = x[i][rows - Ghost_cell - 4][k];
            }
        }   
    }
    else if (strcmp(boundary_type, R) == 0) {
        // Right boundary <- Left interior region
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){
                x[i][rows - Ghost_cell][k] = x[i][Ghost_cell][k];         // Right ghost <- Left interior
                x[i][rows - Ghost_cell + 1][k] = x[i][Ghost_cell + 1][k];
                x[i][rows - Ghost_cell + 2][k] = x[i][Ghost_cell + 2][k];
                x[i][rows - Ghost_cell + 3][k] = x[i][Ghost_cell + 3][k];
            }                   
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // Bottom boundary <- Top interior region
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                x[i][j][Ghost_cell-1] = x[i][j][cols - Ghost_cell - 1];   // Bottom ghost <- Top interior
                x[i][j][Ghost_cell-2] = x[i][j][cols - Ghost_cell - 2];
                x[i][j][Ghost_cell-3] = x[i][j][cols - Ghost_cell - 3];
                x[i][j][Ghost_cell-4] = x[i][j][cols - Ghost_cell - 4];
            }                   
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // Top boundary <- Bottom interior region
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                x[i][j][cols - Ghost_cell] = x[i][j][Ghost_cell];         // Top ghost <- Bottom interior
                x[i][j][cols - Ghost_cell + 1] = x[i][j][Ghost_cell + 1];
                x[i][j][cols - Ghost_cell + 2] = x[i][j][Ghost_cell + 2];
                x[i][j][cols - Ghost_cell + 3] = x[i][j][Ghost_cell + 3];
            }                   
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}

/**
 * Apply reflection (slip-wall) boundary conditions for 2D Euler equations
 * For inviscid flows: normal velocity reverses sign, tangential velocity remains unchanged
 * Assumes variable order: [density, x-momentum, y-momentum, energy, ...]
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param x Solution array [var][rows][cols]
 * @param boundary_type String identifier for boundary: "L" (left), "R" (right), "B" (bottom), "T" (top)
 * @param Ghost_cell Number of ghost cells on each side
 */
static inline void BC_Reflection_2D(int var, int rows, int cols, double (*x)[rows][cols], 
                                    const char* boundary_type, int Ghost_cell){
    // Assumes variable order: [density, x-momentum, y-momentum, energy, ...]
    // For reflection boundary: normal velocity reverses sign, tangential velocity remains unchanged, 
    // other variables (density, energy) remain unchanged
    
    if (strcmp(boundary_type, L) == 0) {
        // Left boundary wall reflection: x-velocity reverses sign, y-velocity unchanged
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){
                if (i == 1) { // x-momentum (normal velocity component for left boundary)
                    x[i][Ghost_cell-1][k] = -x[i][Ghost_cell][k];      // Reverse normal velocity
                    x[i][Ghost_cell-2][k] = -x[i][Ghost_cell+1][k];
                    x[i][Ghost_cell-3][k] = -x[i][Ghost_cell+2][k];
                    x[i][Ghost_cell-4][k] = -x[i][Ghost_cell+3][k];
                } else { // Density, y-momentum, energy remain unchanged
                    x[i][Ghost_cell-1][k] = x[i][Ghost_cell][k];
                    x[i][Ghost_cell-2][k] = x[i][Ghost_cell+1][k];
                    x[i][Ghost_cell-3][k] = x[i][Ghost_cell+2][k];
                    x[i][Ghost_cell-4][k] = x[i][Ghost_cell+3][k];
                }
            }
        }   
    }
    else if (strcmp(boundary_type, R) == 0) {
        // Right boundary wall reflection: x-velocity reverses sign, y-velocity unchanged
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){
                if (i == 1) { // x-momentum (normal velocity component for right boundary)
                    x[i][rows - Ghost_cell][k] = -x[i][rows - Ghost_cell - 1][k];
                    x[i][rows - Ghost_cell + 1][k] = -x[i][rows - Ghost_cell - 2][k];
                    x[i][rows - Ghost_cell + 2][k] = -x[i][rows - Ghost_cell - 3][k];
                    x[i][rows - Ghost_cell + 3][k] = -x[i][rows - Ghost_cell - 4][k];
                } else { // Density, y-momentum, energy remain unchanged
                    x[i][rows - Ghost_cell][k] = x[i][rows - Ghost_cell - 1][k];
                    x[i][rows - Ghost_cell + 1][k] = x[i][rows - Ghost_cell - 2][k];
                    x[i][rows - Ghost_cell + 2][k] = x[i][rows - Ghost_cell - 3][k];
                    x[i][rows - Ghost_cell + 3][k] = x[i][rows - Ghost_cell - 4][k];
                }
            }                   
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // Bottom boundary wall reflection: y-velocity reverses sign, x-velocity unchanged
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                if (i == 2) { // y-momentum (normal velocity component for bottom boundary)
                    x[i][j][Ghost_cell-1] = -x[i][j][Ghost_cell];
                    x[i][j][Ghost_cell-2] = -x[i][j][Ghost_cell+1];
                    x[i][j][Ghost_cell-3] = -x[i][j][Ghost_cell+2];
                    x[i][j][Ghost_cell-4] = -x[i][j][Ghost_cell+3];
                } else { // Density, x-momentum, energy remain unchanged
                    x[i][j][Ghost_cell-1] = x[i][j][Ghost_cell];
                    x[i][j][Ghost_cell-2] = x[i][j][Ghost_cell+1];
                    x[i][j][Ghost_cell-3] = x[i][j][Ghost_cell+2];
                    x[i][j][Ghost_cell-4] = x[i][j][Ghost_cell+3];
                }
            }                   
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // Top boundary wall reflection: y-velocity reverses sign, x-velocity unchanged
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                if (i == 2) { // y-momentum (normal velocity component for top boundary)
                    x[i][j][cols - Ghost_cell] = -x[i][j][cols - Ghost_cell - 1];
                    x[i][j][cols - Ghost_cell + 1] = -x[i][j][cols - Ghost_cell - 2];
                    x[i][j][cols - Ghost_cell + 2] = -x[i][j][cols - Ghost_cell - 3];
                    x[i][j][cols - Ghost_cell + 3] = -x[i][j][cols - Ghost_cell - 4];
                } else { // Density, x-momentum, energy remain unchanged
                    x[i][j][cols - Ghost_cell] = x[i][j][cols - Ghost_cell - 1];
                    x[i][j][cols - Ghost_cell + 1] = x[i][j][cols - Ghost_cell - 2];
                    x[i][j][cols - Ghost_cell + 2] = x[i][j][cols - Ghost_cell - 3];
                    x[i][j][cols - Ghost_cell + 3] = x[i][j][cols - Ghost_cell - 4];
                }
            }                   
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}


/**
 * Boundary condition function for Double Mach Reflection problem
 * Special boundary conditions based on Woodward & Colella (1984) setup
 * This function handles all boundaries in a single call
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param x Solution array [var][rows][cols]
 * @param Ghost_cell Number of ghost cells on each side
 * @param current_time Current simulation time (used for moving shock position)
 */
static inline void BC_DoubleMach_2D(int var, int rows, int cols, double (*x)[rows][cols], 
                                    int Ghost_cell, double current_time) {
    // Physical parameters for Double Mach Reflection problem (Woodward & Colella, 1984)
    double Lx = 4.0;          // Domain length in x-direction
    double Ly = 1.0;          // Domain height in y-direction
    double shock_slope = 1.732;  // tan(60°) for 60-degree incident shock
    double shock_start_x = 1.0/6.0;  // Initial shock position at bottom wall
    double wall_start_x = 1.0/6.0;   // Starting point of reflecting wall
    double gamma = 1.4;       // Specific heat ratio
    
    // Pre-shock state (ahead of shock) - primitive variables
    double pre_shock_rho = 1.4;  // Density
    double pre_shock_u = 0.0;    // x-velocity
    double pre_shock_v = 0.0;    // y-velocity
    double pre_shock_p = 1.0;    // Pressure
    
    // Post-shock state (behind shock) - primitive variables
    double post_shock_rho = 8.0;      // Density
    double post_shock_u = 7.145;      // x-velocity
    double post_shock_v = -4.125;     // y-velocity
    double post_shock_p = 116.83333;  // Pressure
    
    // Convert primitive variables to conservative variables
    // Pre-shock conservative variables
    double pre_rho = pre_shock_rho;
    double pre_rhou = pre_shock_rho * pre_shock_u;
    double pre_rhov = pre_shock_rho * pre_shock_v;
    double pre_E = pre_shock_p/(gamma-1.0) + 0.5*pre_shock_rho*(pre_shock_u*pre_shock_u + pre_shock_v*pre_shock_v);
    
    // Post-shock conservative variables
    double post_rho = post_shock_rho;
    double post_rhou = post_shock_rho * post_shock_u;
    double post_rhov = post_shock_rho * post_shock_v;
    double post_E = post_shock_p/(gamma-1.0) + 0.5*post_shock_rho*(post_shock_u*post_shock_u + post_shock_v*post_shock_v);
    
    // Number of interior cells in x-direction
    double nx = rows - 2 * Ghost_cell;
    
    // ========== 1. Left Boundary: Post-shock inflow ==========
    // The left boundary is completely in the post-shock region
    #pragma omp parallel for collapse(1)
    for (int k = Ghost_cell; k < cols-Ghost_cell; k++){
        // Density (conservative variable)
        x[0][Ghost_cell-1][k] = post_rho;
        x[0][Ghost_cell-2][k] = post_rho;
        x[0][Ghost_cell-3][k] = post_rho;
        x[0][Ghost_cell-4][k] = post_rho;
        
        // x-momentum
        x[1][Ghost_cell-1][k] = post_rhou;
        x[1][Ghost_cell-2][k] = post_rhou;
        x[1][Ghost_cell-3][k] = post_rhou;
        x[1][Ghost_cell-4][k] = post_rhou;
        
        // y-momentum
        x[2][Ghost_cell-1][k] = post_rhov;
        x[2][Ghost_cell-2][k] = post_rhov;
        x[2][Ghost_cell-3][k] = post_rhov;
        x[2][Ghost_cell-4][k] = post_rhov;
        
        // Total energy
        x[3][Ghost_cell-1][k] = post_E;
        x[3][Ghost_cell-2][k] = post_E;
        x[3][Ghost_cell-3][k] = post_E;
        x[3][Ghost_cell-4][k] = post_E;
    }
    
    // ========== 2. Right Boundary: Outflow ==========
    // Zero-gradient (Neumann) condition for outflow
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < var; i++){
        for (int k = Ghost_cell; k < cols-Ghost_cell; k++){
            x[i][rows-Ghost_cell][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+1][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+2][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+3][k] = x[i][rows-Ghost_cell-1][k];
        }                   
    }
    
    // ========== 3. Bottom Boundary: Mixed boundary ==========
    // Part reflecting wall (x > wall_start_x), part post-shock inflow (x < wall_start_x)
    #pragma omp parallel for collapse(1)
    for (int j = Ghost_cell; j < rows-Ghost_cell; j++){
        // Calculate physical coordinate (cell center)
        double x_pos = (double)(j+0.5-Ghost_cell)/ nx * Lx;
        
        if (x_pos >= wall_start_x) {
            // Reflecting wall region (right portion of bottom boundary)
            
            // Density - extrapolate from interior
            x[0][j][Ghost_cell-1] = x[0][j][Ghost_cell];
            x[0][j][Ghost_cell-2] = x[0][j][Ghost_cell+1];
            x[0][j][Ghost_cell-3] = x[0][j][Ghost_cell+2];
            x[0][j][Ghost_cell-4] = x[0][j][Ghost_cell+3];
            
            // x-momentum - unchanged (tangential component)
            x[1][j][Ghost_cell-1] = x[1][j][Ghost_cell];
            x[1][j][Ghost_cell-2] = x[1][j][Ghost_cell+1];
            x[1][j][Ghost_cell-3] = x[1][j][Ghost_cell+2];
            x[1][j][Ghost_cell-4] = x[1][j][Ghost_cell+3];
            
            // y-momentum - reverse sign (normal component: ρv → -ρv)
            x[2][j][Ghost_cell-1] = -x[2][j][Ghost_cell];
            x[2][j][Ghost_cell-2] = -x[2][j][Ghost_cell+1];
            x[2][j][Ghost_cell-3] = -x[2][j][Ghost_cell+2];
            x[2][j][Ghost_cell-4] = -x[2][j][Ghost_cell+3];
            
            // Total energy - unchanged
            x[3][j][Ghost_cell-1] = x[3][j][Ghost_cell];
            x[3][j][Ghost_cell-2] = x[3][j][Ghost_cell+1];
            x[3][j][Ghost_cell-3] = x[3][j][Ghost_cell+2];
            x[3][j][Ghost_cell-4] = x[3][j][Ghost_cell+3];
            
        } else {
            // Non-wall region (left portion of bottom boundary) - post-shock inflow
           
            // Density
            x[0][j][Ghost_cell-1] = post_rho;
            x[0][j][Ghost_cell-2] = post_rho;
            x[0][j][Ghost_cell-3] = post_rho;
            x[0][j][Ghost_cell-4] = post_rho;
            
            // x-momentum
            x[1][j][Ghost_cell-1] = post_rhou;
            x[1][j][Ghost_cell-2] = post_rhou;
            x[1][j][Ghost_cell-3] = post_rhou;
            x[1][j][Ghost_cell-4] = post_rhou;
            
            // y-momentum
            x[2][j][Ghost_cell-1] = post_rhov;
            x[2][j][Ghost_cell-2] = post_rhov;
            x[2][j][Ghost_cell-3] = post_rhov;
            x[2][j][Ghost_cell-4] = post_rhov;
            
            // Total energy
            x[3][j][Ghost_cell-1] = post_E;
            x[3][j][Ghost_cell-2] = post_E;
            x[3][j][Ghost_cell-3] = post_E;
            x[3][j][Ghost_cell-4] = post_E;
        }
    }
    
    // ========== 4. Top Boundary: Shock tracking boundary ==========
    // Calculate shock intersection with top boundary (shock moves with time)
    // shock_x_at_top = initial_x + (Ly/tanθ) + (shock_speed * time / tanθ)
    // Shock speed = 20.0 for Mach 10 shock (from Woodward & Colella, 1984)
    double shock_x_at_top = shock_start_x + (Ly / shock_slope) + (20.0 * current_time / shock_slope);
    
    #pragma omp parallel for collapse(1)
    for (int j = Ghost_cell; j < rows-Ghost_cell; j++){
        // Calculate physical coordinate (cell center)
        double x_pos = (double)(j + 0.5 - Ghost_cell) / nx * Lx;
        
        if (x_pos < shock_x_at_top) {
            // Post-shock region (left of shock intersection)
            
            // Density
            x[0][j][cols-Ghost_cell] = post_rho;
            x[0][j][cols-Ghost_cell+1] = post_rho;
            x[0][j][cols-Ghost_cell+2] = post_rho;
            x[0][j][cols-Ghost_cell+3] = post_rho;
            
            // x-momentum
            x[1][j][cols-Ghost_cell] = post_rhou;
            x[1][j][cols-Ghost_cell+1] = post_rhou;
            x[1][j][cols-Ghost_cell+2] = post_rhou;
            x[1][j][cols-Ghost_cell+3] = post_rhou;
            
            // y-momentum
            x[2][j][cols-Ghost_cell] = post_rhov;
            x[2][j][cols-Ghost_cell+1] = post_rhov;
            x[2][j][cols-Ghost_cell+2] = post_rhov;
            x[2][j][cols-Ghost_cell+3] = post_rhov;
            
            // Total energy
            x[3][j][cols-Ghost_cell] = post_E;
            x[3][j][cols-Ghost_cell+1] = post_E;
            x[3][j][cols-Ghost_cell+2] = post_E;
            x[3][j][cols-Ghost_cell+3] = post_E;
            
        } else {
            // Pre-shock region (right of shock intersection)
            
            // Density
            x[0][j][cols-Ghost_cell] = pre_rho;
            x[0][j][cols-Ghost_cell+1] = pre_rho;
            x[0][j][cols-Ghost_cell+2] = pre_rho;
            x[0][j][cols-Ghost_cell+3] = pre_rho;
            
            // x-momentum
            x[1][j][cols-Ghost_cell] = pre_rhou;
            x[1][j][cols-Ghost_cell+1] = pre_rhou;
            x[1][j][cols-Ghost_cell+2] = pre_rhou;
            x[1][j][cols-Ghost_cell+3] = pre_rhou;
            
            // y-momentum
            x[2][j][cols-Ghost_cell] = pre_rhov;
            x[2][j][cols-Ghost_cell+1] = pre_rhov;
            x[2][j][cols-Ghost_cell+2] = pre_rhov;
            x[2][j][cols-Ghost_cell+3] = pre_rhov;
            
            // Total energy
            x[3][j][cols-Ghost_cell] = pre_E;
            x[3][j][cols-Ghost_cell+1] = pre_E;
            x[3][j][cols-Ghost_cell+2] = pre_E;
            x[3][j][cols-Ghost_cell+3] = pre_E;
         
        }
    }
}


// Boundary Condition Function For Two Dimensional Backward Step Flow According To Woodward Colella 1984
static inline void BC_BackwardStep_2D(int var, int rows1, int cols1, double (*x)[rows1][cols1],int rows2, int cols2, double (*y)[rows2][cols2]) {


    // ========== Block1左边界：超声速入流 ==========
    // 整个左边界都是固定来流条件
    for (int k = GhostCell; k < cols1-GhostCell; k++) {
        // 密度
        x[0][GhostCell-1][k] = 1.4;
        x[0][GhostCell-2][k] = 1.4;
        x[0][GhostCell-3][k] = 1.4;
        x[0][GhostCell-4][k] = 1.4;
        
        // x动量
        x[1][GhostCell-1][k] = 1.4 * 3.0;
        x[1][GhostCell-2][k] = 1.4 * 3.0;
        x[1][GhostCell-3][k] = 1.4 * 3.0;
        x[1][GhostCell-4][k] = 1.4 * 3.0;
        
        // y动量
        x[2][GhostCell-1][k] = 0.0;
        x[2][GhostCell-2][k] = 0.0;
        x[2][GhostCell-3][k] = 0.0;
        x[2][GhostCell-4][k] = 0.0;
        
        // 总能量
        x[3][GhostCell-1][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
        x[3][GhostCell-2][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
        x[3][GhostCell-3][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
        x[3][GhostCell-4][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
    }
    
    // ========== block1右边界：超声速出流 ==========
    // 零梯度外推
    for (int i = 0; i < var; i++) {
        for (int k = GhostCell; k < cols1-GhostCell; k++) {
            x[i][rows1-GhostCell][k] = x[i][rows1-GhostCell-1][k];
            x[i][rows1-GhostCell+1][k] = x[i][rows1-GhostCell-1][k];
            x[i][rows1-GhostCell+2][k] = x[i][rows1-GhostCell-1][k];
            x[i][rows1-GhostCell+3][k] = x[i][rows1-GhostCell-1][k];
        }
    }
    
    // ========== block1下边界：反射壁面 ==========
    // 整个边界都是反射壁面
    for (int i = 0; i < var; i++){
        for (int j = GhostCell; j < rows1 - GhostCell; j++){
            if (i == 2) { // y方向动量（速度反向）
                x[i][j][GhostCell-1] = -x[i][j][GhostCell];
                x[i][j][GhostCell-2] = -x[i][j][GhostCell+1];
                x[i][j][GhostCell-3] = -x[i][j][GhostCell+2];
                x[i][j][GhostCell-4] = -x[i][j][GhostCell+3];
            } else { // 密度、x方向动量、能量等保持不变
                x[i][j][GhostCell-1] = x[i][j][GhostCell];
                x[i][j][GhostCell-2] = x[i][j][GhostCell+1];
                x[i][j][GhostCell-3] = x[i][j][GhostCell+2];
                x[i][j][GhostCell-4] = x[i][j][GhostCell+3];
            }
        }                   
    }
    
    // ========== block1上边界：反射壁面 ==========
    // 整个上边界都是反射壁面
    for (int i = 0; i < var; i++){
        for (int j = GhostCell; j < rows1 - GhostCell; j++){
            if (i == 2) { // y方向动量（速度反向）
                x[i][j][cols1 - GhostCell] = -x[i][j][cols1 - GhostCell - 1];
                x[i][j][cols1 - GhostCell + 1] = -x[i][j][cols1 - GhostCell - 2];
                x[i][j][cols1 - GhostCell + 2] = -x[i][j][cols1 - GhostCell - 3];
                x[i][j][cols1 - GhostCell + 3] = -x[i][j][cols1 - GhostCell - 4];
            } else { // 密度、x方向动量、能量等保持不变
                x[i][j][cols1 - GhostCell] = x[i][j][cols1 - GhostCell - 1];
                x[i][j][cols1 - GhostCell + 1] = x[i][j][cols1 - GhostCell - 2];
                x[i][j][cols1 - GhostCell + 2] = x[i][j][cols1 - GhostCell - 3];
                x[i][j][cols1 - GhostCell + 3] = x[i][j][cols1 - GhostCell - 4];
            }
        }                   
    }

    // ========== Block2左边界：超声速入流 ==========
    // 整个左边界都是固定来流条件
    for (int k = GhostCell; k < cols2 - GhostCell; k++) {
        // 密度
        y[0][GhostCell-1][k] = 1.4;
        y[0][GhostCell-2][k] = 1.4;
        y[0][GhostCell-3][k] = 1.4;
        y[0][GhostCell-4][k] = 1.4;
        
        // x动量
        y[1][GhostCell-1][k] = 1.4 * 3.0;
        y[1][GhostCell-2][k] = 1.4 * 3.0;
        y[1][GhostCell-3][k] = 1.4 * 3.0;
        y[1][GhostCell-4][k] = 1.4 * 3.0;
        
        // y动量
        y[2][GhostCell-1][k] = 0.0;
        y[2][GhostCell-2][k] = 0.0;
        y[2][GhostCell-3][k] = 0.0;
        y[2][GhostCell-4][k] = 0.0;
        
        // 总能量
        y[3][GhostCell-1][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
        y[3][GhostCell-2][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
        y[3][GhostCell-3][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
        y[3][GhostCell-4][k] = 1.0/(M_gamma-1)+0.5*1.4*3*3;
    }
    
    // ========== block2右边界：反射边界 ==========
    // 零梯度外推
    for (int i = 0; i < var; i++){
            for (int k = GhostCell; k < cols2 - GhostCell; k++){
                if (i == 1) { // x方向动量（速度反向）
                    y[i][rows2 - GhostCell][k] = -y[i][rows2 -GhostCell - 1][k];
                    y[i][rows2 - GhostCell + 1][k] = -y[i][rows2 -GhostCell - 2][k];
                    y[i][rows2 - GhostCell + 2][k] = -y[i][rows2 -GhostCell - 3][k];
                    y[i][rows2 - GhostCell + 3][k] = -y[i][rows2 -GhostCell - 4][k];
                } else { // 密度、y方向动量、能量等保持不变
                    y[i][rows2 - GhostCell][k] = y[i][rows2 -GhostCell - 1][k];
                    y[i][rows2 - GhostCell + 1][k] = y[i][rows2 -GhostCell - 2][k];
                    y[i][rows2 - GhostCell + 2][k] = y[i][rows2 -GhostCell - 3][k];
                    y[i][rows2 - GhostCell + 3][k] = y[i][rows2 -GhostCell - 4][k];
                }
            }                   
        }
    
    // ========== block2下边界：反射壁面 ==========
    // 整个下边界都是反射壁面
    for (int i = 0; i < var; i++){
        for (int j = GhostCell; j < rows2 - GhostCell; j++){
            if (i == 2) { // y方向动量（速度反向）
                y[i][j][GhostCell-1] = -y[i][j][GhostCell];
                y[i][j][GhostCell-2] = -y[i][j][GhostCell+1];
                y[i][j][GhostCell-3] = -y[i][j][GhostCell+2];
                y[i][j][GhostCell-4] = -y[i][j][GhostCell+3];
            } else { // 密度、x方向动量、能量等保持不变
                y[i][j][GhostCell-1] = y[i][j][GhostCell];
                y[i][j][GhostCell-2] = y[i][j][GhostCell+1];
                y[i][j][GhostCell-3] = y[i][j][GhostCell+2];
                y[i][j][GhostCell-4] = y[i][j][GhostCell+3];
            }
        }                   
    }
    
    // ========== block2上边界,处理相互作用==========
    // block2上边界的相互作用边界
    for (int i = 0; i < var; i++){
        for (int j = GhostCell; j < rows2 - GhostCell; j++){
            //处理Block2
            y[i][j][cols2 - GhostCell] = x[i][j][GhostCell];
            y[i][j][cols2 - GhostCell + 1] = x[i][j][GhostCell + 1];
            y[i][j][cols2 - GhostCell + 2] = x[i][j][GhostCell + 2];
            y[i][j][cols2 - GhostCell + 3] = x[i][j][GhostCell + 3];

            //处理Block1
            x[i][j][0] = y[i][j][cols2 - 2 * GhostCell];
            x[i][j][1] = y[i][j][cols2 - 2 * GhostCell + 1];
            x[i][j][2] = y[i][j][cols2 - 2 * GhostCell + 2];
            x[i][j][3] = y[i][j][cols2 - 2 * GhostCell + 3];

        }                   
    }

}


/**
 * Determine boundary condition configuration for each block in multi-block decomposition
 * This function sets the appropriate boundary condition types for each block face
 * 
 * @param test_case Test case identifier from TestCase2D enumeration
 * @param Block_name Block identifier (0 for block 0, 1 for block 1, etc.)
 */
static inline void Get_block_BC(int test_case, int Block_name) {
    
    switch (test_case) {
        // Test cases that require block decomposition
        case TEST_BACKWARD_STEP:
            // Backward step flow: two blocks (main flow and step region)
            switch (Block_name) {
                case 0:  // Main flow region (above the step)
                    bc_config.left = BC_FIXED_VALUE;      // Inflow boundary (Mach 3 inflow)
                    bc_config.right = BC_OUTFLOW;         // Outflow boundary
                    bc_config.bottom = BC_REFLECTION;     // Top of step wall (reflection)
                    bc_config.top = BC_REFLECTION;        // Channel top wall (reflection)
                    break;
                
                case 1:  // Step region (below the step)
                    bc_config.left = BC_FIXED_VALUE;      // Inflow boundary (should match main region)
                    bc_config.right = BC_REFLECTION;      // Step face (vertical wall, reflection)
                    bc_config.bottom = BC_REFLECTION;     // Channel bottom wall (reflection)
                    bc_config.top = BC_OUTFLOW;           // Interface with main region (outflow)
                    break;
                    
                default:
                    printf("Warning: Block_name value %d out of range (should be 0 or 1)\n", Block_name);
                    printf("Location: Get_block_BC function, test_case=%d\n", test_case);
                    break;
            }
            break;
            
        default:
            // For test cases without block decomposition, no action needed
            break; 
    }
}

/**
 * Exchange data between overlapping regions of neighboring blocks
 * This function transfers solution data from overlap region to block ghost cells
 * 
 * @param var Number of variables (typically 4: density, x-momentum, y-momentum, energy)
 * @param Block_rows Number of rows in the block (including ghost cells)
 * @param Block_cols Number of columns in the block (including ghost cells)
 * @param BU Block solution array [var][Block_rows][Block_cols]
 * @param Overlap_X Number of cells in x-direction in overlap region
 * @param Overlap_Y Number of cells in y-direction in overlap region
 * @param OverlapU Overlap region solution array [var][Overlap_X][Overlap_Y]
 * @param test_case Test case identifier from TestCase2D enumeration
 * @param Block_name Block identifier (0 for block 0, 1 for block 1, etc.)
 */
static inline void Block_DataExchange(int var, int Block_rows, int Block_cols, double (*BU)[Block_rows][Block_cols],
                                        int Overlap_X, int Overlap_Y, double (*OverlapU)[Overlap_X][Overlap_Y],
                                        int test_case, int Block_name) {
    
    switch (test_case) {
        case TEST_BACKWARD_STEP:
            // Data exchange for backward step problem
            switch (Block_name) {
                case 0: 
                    // Block 0 (main region): receives data from block 1 at bottom boundary
                    for (int k = 0; k < var; k++) {  // Loop over variables
                        for (int i = GhostCell; i < Overlap_X + GhostCell; i++) {  // Loop over x-direction (overlap region)
                            for (int j = 0; j < GhostCell; j++) {  // Loop over ghost cells in y-direction
                                // Transfer data from overlap region to block's bottom ghost cells
                                BU[k][i][j] = OverlapU[k][i - GhostCell][j];
                            }
                        }                   
                    }
                    break;
                
                case 1: 
                    // Block 1 (step region): receives data from block 0 at top boundary
                    for (int k = 0; k < var; k++) {  // Loop over variables
                        for (int i = GhostCell; i < Overlap_X + GhostCell; i++) {  // Loop over x-direction (overlap region)
                            for (int j = 0; j < GhostCell; j++) {  // Loop over ghost cells in y-direction
                                // Transfer data from overlap region to block's top ghost cells
                                // Note: OverlapU contains data starting from j=GhostCell for proper alignment
                                BU[k][i][Block_cols - GhostCell + j] = OverlapU[k][i - GhostCell][j + GhostCell];
                            }
                        }                   
                    }
                    break;
                    
                default:
                    printf("Warning: Block_name value %d out of range (should be 0 or 1)\n", Block_name);
                    printf("Location: Get_block_BC function, test_case=%d\n", test_case);
                    break;
            }
            break;
            
        default:
            // For test cases without block decomposition, no data exchange needed
            break; 
    }
}


// Update The Boundary Condition Function Pointer Array
static inline void Boundary_Conditions(int var, int rows, int cols, double (*y)[rows][cols], int GC) {

    void (*bc_funcs[])(int, int, int, double (*)[*][*], const char*, int) = {
        BC_OutFlow_2D,      // 0: BC_OUTFLOW
        BC_Reflection_2D,   // 1: BC_REFLECTION
        BC_Periodicity_2D,  // 2: BC_PERIODICITY
        BC_Inflow_2D,       // 3: BC_INFLOW
        BC_FixedValue_2D    // 4: BC_FIXED_VALUE
    };
    
    const char* sides[] = {L, R, B, T};
    BC_Type bc_types[] = {bc_config.left, bc_config.right, bc_config.bottom, bc_config.top};
    
    for (int i = 0; i < 4; i++) {
        bc_funcs[bc_types[i]](var, rows, cols, y, sides[i], GC);
    }
}


#endif // BOUNDARY_CONDITION_H