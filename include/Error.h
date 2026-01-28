#ifndef ERROR_H
#define ERROR_H

#include <stdio.h>
#include <math.h>
#include "Reconstruction.h"
#include "function.h"
#include "Golbal.h"


// Structure to store error data for convergence analysis
typedef struct {
    int grid_points;        // Number of grid points (Nx)
    double delta_x;         // Grid spacing (Δx)
    double L1_error;        // L1 error norm
    double L2_error;        // L2 error norm
    double Linf_error;      // L-infinity (L∞) error norm
    double L1_order;        // Convergence order calculated from L1 error
    double L2_order;        // Convergence order calculated from L2 error
    double Linf_order;      // Convergence order calculated from L∞ error
} ErrorData;

// Global array to store error data for different grid resolutions
ErrorData error_table[100];  // Storage for error data from different mesh refinements (max 100 levels)
int error_count = 0;         // Current number of processed grid resolutions


//计算有限体积格式误差L1误差
static inline void Scheme_Error_L1(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    

    for (int i = 0; i < rows; i++)
    {
        Error[i] = 0.0;
    }
    

    for (int i = 0; i < rows; i++){
        double E1 = 0;
        for (int j = GC ; j <= cols-GC-2; j++){
            double u[10];
            for (int h = 0; h < 10; h++)
                u[h] = x[i][j+h-2];

            double Inter = delta_x * numerical_integration(Smooth_function,(j-GC)*delta_x, (j-GC+1.0)*delta_x, 1000, 2, &u[2]);
            E1 = E1 + Inter;
        } 
        Error[i] = E1;
    }

    printf("Mass_Error = %f \n", Error[0]);
}


//计算有限体积格式误差L2误差
static inline void Scheme_Error_L2(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    
    int l =  cols- 2*GC;

    for (int i = 0; i < rows; i++)
    {
        Error[i] = 0.0;
    }
    

//    WENO3_Reconstruction(rows,cols,GC,x,Error_L,Error_R);

    for (int i = 0; i < rows; i++){
        double E1 = 0;
        double Error_max = 0.0;
        for (int j = GC ; j <= cols-GC-2; j++){
            Error_max = max_of_two(Error_max, fabs(x[i][j] + 1.0/( 4.0 * PI )*(cos(4.0 * PI * (j+1-GC) * delta_x)-cos(4.0 * PI * (j-GC) * delta_x))/delta_x - 2.0));
 //          E1 = E1 + pow(x[i][j] * delta_x +  1.0/( 4.0 * PI )*(cos(4.0 * PI * (j+1-GC) * delta_x)-cos(4.0 * PI * (j-GC) * delta_x)) - 2.0 *delta_x,2) ;
            E1 = E1 + delta_x * pow(sin(4.0 * PI * (j-GC) * delta_x) + 2.0 - Error_R[i][j-1], 2);
        } 
        Error[i] = Error_max;
    }

    printf("Mass_Error = %f \n", Error[0]);
}

//计算有限体积格式误差L∞误差
static inline void Scheme_Error_LinFin(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    
    int l =  cols- 2*GC;

    for (int i = 0; i < rows; i++)
        Error[i] = 0.0;


    for (int i = 0; i < rows; i++){
        double E1 = 0;
        double Error_max = 0.0;
        for (int j = GC ; j <= cols-GC-2; j++){
            Error_max = max_of_two(Error_max, fabs(x[i][j] + 1.0/( 4.0 * PI )*(cos(4.0 * PI * (j+1-GC) * delta_x)-cos(4.0 * PI * (j-GC) * delta_x))/delta_x - 2.0));
        } 
        Error[i] = Error_max;
    }

    printf("Mass_Error = %f \n", Error[0]);
}



/**
 * Calculate and store L1 error for density in 1D Euler equations test case
 * L1 error = Δx * Σ|ρ_exact - ρ_numerical|, where ρ_numerical is cell average
 * Exact solution: ρ(x,t) = 1 + 0.2*sin(x - t)
 * 
 * @param var Number of variables (unused in this function but kept for consistency)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells at each boundary
 * @param Conser Conservative variables array [var][rows][cols], Conser[0][i][GC] stores density ρ
 * @param Delta Grid spacing in x-direction (Δx)
 * @param Time Current simulation time (t)
 * @param grid_points Number of grid points in x-direction (without ghost cells)
 * @param store_index Index in error_table array for storing error data
 */
static inline void TEST_Scheme_Error_L1(int var, int rows, int cols, int GC, 
                               double (*Conser)[rows][cols], double Delta, double Time,
                               int grid_points, int store_index) {
    
    double Mass_Error = 0.0;
    
    // Calculate exact solution at cell centers
    double Exact[rows];
    for (int i = GC; i < rows-GC; i++) {
        // x coordinate at cell center
        double x_center = (i - GC + 0.5) * Delta;
        // Exact solution at cell center: ρ(x,t) = 1 + 0.2*sin(x - t)
        Exact[i] = 1.0 + 0.2 * sin(x_center - Time);
    }
    
    // Calculate L1 error: Δx * Σ|ρ_exact - ρ_numerical|
    for (int i = GC; i < rows-GC; i++) {
        Mass_Error += fabs(Exact[i] - Conser[0][i][GC]) * Delta;
    }
    
    // Store error data for later convergence analysis
    error_table[store_index].grid_points = grid_points;
    error_table[store_index].delta_x = Delta;
    error_table[store_index].L1_error = Mass_Error;
    
    printf("L1 Error = %.10e \n", Mass_Error);
}

/**
 * Calculate and store L2 error for density in 1D Euler equations test case
 * L2 error = sqrt(Δx * Σ(ρ_exact - ρ_numerical)^2)
 * Exact solution: ρ(x,t) = 1 + 0.2*sin(x - t)
 * 
 * @param var Number of variables (unused in this function but kept for consistency)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells at each boundary
 * @param Conser Conservative variables array [var][rows][cols], Conser[0][i][GC] stores density ρ
 * @param Delta Grid spacing in x-direction (Δx)
 * @param Time Current simulation time (t)
 * @param store_index Index in error_table array for storing error data
 */
static inline void TEST_Scheme_Error_L2(int var, int rows, int cols, int GC, 
                               double (*Conser)[rows][cols], double Delta, double Time,
                               int store_index) {
    
    double Mass_Error = 0.0;
    
    // Calculate L2 error: sqrt(Δx * Σ(ρ_exact - ρ_numerical)^2)
    for (int i = GC; i < rows-GC; i++) {
        // x coordinate at cell center
        double x_center = (i - GC + 0.5) * Delta;
        // Exact solution at cell center
        double rho_exact = 1.0 + 0.2 * sin(x_center - Time);
        // Numerical solution (cell average of density)
        double rho_num = Conser[0][i][GC];
        
        // Accumulate squared error
        double error = rho_exact - rho_num;
        Mass_Error += error * error * Delta;
    }
    
    // Take square root to get L2 norm
    Mass_Error = sqrt(Mass_Error);
    
    // Store L2 error
    error_table[store_index].L2_error = Mass_Error;
    
    printf("L2 Error = %.10e \n", Mass_Error);
}

/**
 * Calculate and store L-infinity error for density in 1D Euler equations test case
 * L∞ error = max|ρ_exact - ρ_numerical| across all cells
 * Exact solution: ρ(x,t) = 1 + 0.2*sin(x - t)
 * 
 * @param var Number of variables (unused in this function but kept for consistency)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells at each boundary
 * @param Conser Conservative variables array [var][rows][cols], Conser[0][i][GC] stores density ρ
 * @param Delta Grid spacing in x-direction (Δx)
 * @param Time Current simulation time (t)
 * @param store_index Index in error_table array for storing error data
 */
static inline void TEST_Scheme_Error_Linf(int var, int rows, int cols, int GC, 
                                  double (*Conser)[rows][cols], double Delta, double Time,
                                  int store_index) {
    
    double max_error = 0.0;
    
    // Find maximum absolute error across all cells
    for (int i = GC; i < rows-GC; i++) {
        // x coordinate at cell center
        double x_center = (i - GC + 0.5) * Delta;
        // Exact solution at cell center
        double rho_exact = 1.0 + 0.2 * sin(x_center - Time);
        // Numerical solution (cell average of density)
        double rho_num = Conser[0][i][GC];
        
        // Calculate absolute error
        double error = fabs(rho_exact - rho_num);
        
        // Update maximum error if current error is larger
        if (error > max_error) {
            max_error = error;
        }
    }
    
    // Store L-infinity error
    error_table[store_index].Linf_error = max_error;
    
    printf("L∞ Error = %.10e \n", max_error);
}

/**
 * Calculate convergence orders based on errors from successive grid refinements
 * Convergence order = log10(error_coarse/error_fine) / log10(2) when grid is doubled
 * This function populates the order fields in error_table array
 * 
 * Note: First grid (index 0) has no order since there's no coarser grid for comparison
 */
/**
 * Calculate convergence orders based on errors from successive grid refinements
 * Convergence order = ln(error_coarse/error_fine) / ln(h_coarse/h_fine)
 * where h = Δx is the grid spacing
 * For fixed domain length L, h = L/N, so h_coarse/h_fine = N_fine/N_coarse
 * 
 * Note: First grid (index 0) has no order since there's no coarser grid for comparison
 */
static inline void calculate_orders() {
    
    for (int i = 1; i < error_count; i++) {
        // Calculate the grid spacing ratio: h_coarse / h_fine = N_fine / N_coarse
        double h_ratio = (double)error_table[i].grid_points / (double)error_table[i-1].grid_points;
        
        // Calculate L1 convergence order: ln(error_coarse/error_fine) / ln(h_coarse/h_fine)
        if (error_table[i-1].L1_error > 0 && error_table[i].L1_error > 0) {
            error_table[i].L1_order = log(error_table[i-1].L1_error / error_table[i].L1_error) / 
                                      log(h_ratio);
        } else {
            error_table[i].L1_order = 0.0;
        }
        
        // Calculate L2 convergence order
        if (error_table[i-1].L2_error > 0 && error_table[i].L2_error > 0) {
            error_table[i].L2_order = log(error_table[i-1].L2_error / error_table[i].L2_error) / 
                                      log(h_ratio);
        } else {
            error_table[i].L2_order = 0.0;
        }
        
        // Calculate L-infinity convergence order
        if (error_table[i-1].Linf_error > 0 && error_table[i].Linf_error > 0) {
            error_table[i].Linf_order = log(error_table[i-1].Linf_error / error_table[i].Linf_error) / 
                                        log(h_ratio);
        } else {
            error_table[i].Linf_order = 0.0;
        }
        
        // Optional: Print debugging information
        printf("Order calculation for grid %d (N=%d) to %d (N=%d): h_ratio = %.4f, L1_order = %.4f\n",
               i-1, error_table[i-1].grid_points, i, error_table[i].grid_points,
               h_ratio, error_table[i].L1_order);
    }
}

/**
 * Print error table in a formatted way for publication
 * Format follows Table 3.3 from the reference paper:
 * Grid points, Δx, L1 error, Order, L2 error, Order, L∞ error, Order
 * 
 * Note: First row shows "---" for orders since no coarser grid is available
 */
static inline void print_error_table() {
    printf("\n\n");
    printf("===============================================================================\n");
    printf("Table 3.3: Errors and numerical orders of accuracy of the density\n");
    printf("===============================================================================\n");
    printf("Grid points  Δx        L1 error    Order     L2 error    Order     L∞ error    Order\n");
    printf("------------------------------------------------------------------------------\n");
    
    // First row: no convergence order available (coarsest grid)
    printf("%-11d  %-9.6f  %-10.2e  %-8s  %-10.2e  %-8s  %-10.2e  %-8s\n",
           error_table[0].grid_points,
           error_table[0].delta_x,
           error_table[0].L1_error, "---",
           error_table[0].L2_error, "---",
           error_table[0].Linf_error, "---");
    
    // Subsequent rows: include convergence orders
    for (int i = 1; i < error_count; i++) {
        printf("%-11d  %-9.6f  %-10.2e  %-8.2f  %-10.2e  %-8.2f  %-10.2e  %-8.2f\n",
               error_table[i].grid_points,
               error_table[i].delta_x,
               error_table[i].L1_error, error_table[i].L1_order,
               error_table[i].L2_error, error_table[i].L2_order,
               error_table[i].Linf_error, error_table[i].Linf_order);
    }
    printf("===============================================================================\n");
}

#endif
