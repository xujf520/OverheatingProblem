#ifndef CFD_CONVECTION_H
#define CFD_CONVECTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Golbal.h"
#include "Reconstruction.h"
#include "scheme.h"
#include "function.h"

/*                                      ******************                                          */
/*                          Reconstruction Step Based on Conservative Variables                     */
/*                                      ******************                                          */

                                    /*……………………………………………………*/
                                /*Approximate Riemann Solver*/
                                    /*……………………………………………………*/

/**
 * Flux reconstruction using approximate Riemann solvers
 * This function performs spatial reconstruction and flux calculation in both x and y directions
 * 
 * @param AR_scheme Approximate Riemann solver scheme identifier:
 *                  1: HLL, 2: HLLC, 3: Roe, 11: HLLHC, 22: HLLCHC, 33: RoeHC
 * @param var Number of variables (typically 4 for Euler equations)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells
 * @param y Solution array in conservative variables [var][rows][cols]
 * @param f Flux array in x-direction [var][rows][cols]
 * @param g Flux array in y-direction [var][rows][cols]
 * @param dt Time step (currently not used but kept for compatibility)
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 */
static inline void Flux_Reconstruction_RP(int AR_scheme, int var, int rows, int cols, int GC, double (*y)[rows][cols],
                                          double (*f)[rows][cols], double (*g)[rows][cols], double dt, double dx, double dy) {
    int i, j, k;

    // Allocate memory for temporary arrays
    double (*Flux)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conserl)[rows][cols] = malloc(var * sizeof(double[rows][cols]));  // Left reconstructed states
    double (*Conserr)[rows][cols] = malloc(var * sizeof(double[rows][cols]));  // Right reconstructed states
    
    // Check if memory allocation was successful
    if (Flux == NULL || Conserl == NULL || Conserr == NULL) {
        fprintf(stderr, "Memory allocation failed in Flux_Reconstruction_RP\n");
        // Free any allocated memory before returning
        free(Flux);
        free(Conserl);
        free(Conserr);
        return;
    }

    // Dimensional splitting: process x-direction first, then y-direction
    
    // ========== X-Direction Reconstruction ==========
    int space_dir = 1;  // 1 for x-direction, 2 for y-direction
    
    // Select reconstruction method based on accuracy order
    switch (Recon_Accur) {
        case 1:  // First-order Godunov
            Reconstruction_Godunov(space_dir, var, rows, cols, GC, y, Conserl, Conserr, dx);
            break;
        case 2:  // Second-order TVD
            TVD_Reconstruction(space_dir, var, rows, cols, GC, y, Conserl, Conserr, dx, dy);
            break;
        case 3:  // Third-order WENO
            WENO3_Reconstruction(space_dir, var, rows, cols, GC, y, Conserl, Conserr);
            break;
        case 5:  // Fifth-order WENO
            WENO5_Reconstruction(space_dir, var, rows, cols, GC, y, Conserl, Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n", Recon_Accur);
            exit(1);
    }

    // Evolution process: compute fluxes using approximate Riemann solver
    // AR_scheme is Approximate Riemann Solver identifier
    switch (AR_scheme) {
        case 1:   // HLL Riemann solver
            HLL_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 2:   // HLLC Riemann solver
            HLLC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 3:   // Roe Riemann solver
            Roe_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 11:  // HLLHC (HLL with Heat Conduction)
            HLLHC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 22:  // HLLCHC (HLLC with Heat Conduction)
            HLLCHC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 33:  // RoeHC (Roe with Heat Conduction)
            RoeHC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 5:   // RS_Marquina (not implemented in current code)
            // RS_Marquina(3, rows, y, Flux, dt, dx);
            break;
        default:
            // ER_Flux(var, rows, GC, prileft, priright, Flux);
            // Add appropriate handling logic based on actual requirements
            break;
    }

    // Store x-direction fluxes in f array (for interior cells including one ghost cell)
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < var; i++) {
        for (int j = GC - 1; j <= rows - GC; j++) {
            for (int k = GC - 1; k <= cols - GC; k++) {
                f[i][j][k] = Flux[i][j][k];
            }
        }
    }

    // ========== Y-Direction Reconstruction ==========
    space_dir = 2;  // Switch to y-direction
    
    // Select reconstruction method based on accuracy order
    switch (Recon_Accur) {
        case 1:  // First-order Godunov
            Reconstruction_Godunov(space_dir, var, rows, cols, GC, y, Conserl, Conserr, dx);
            break;
        case 2:  // Second-order TVD
            TVD_Reconstruction(space_dir, var, rows, cols, GC, y, Conserl, Conserr, dx, dy);
            break;
        case 3:  // Third-order WENO
            WENO3_Reconstruction(space_dir, var, rows, cols, GC, y, Conserl, Conserr);
            break;
        case 5:  // Fifth-order WENO
            WENO5_Reconstruction(space_dir, var, rows, cols, GC, y, Conserl, Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n", Recon_Accur);
            exit(1);
    }

    
    // Evolution process: compute fluxes using approximate Riemann solver
    // AR_scheme is Approximate Riemann Solver identifier
    switch (AR_scheme) {
        case 1:   // HLL Riemann solver
            HLL_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 2:   // HLLC Riemann solver
            HLLC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 3:   // Roe Riemann solver
            Roe_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 11:  // HLLHC (HLL with Heat Conduction)
            HLLHC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 22:  // HLLCHC (HLLC with Heat Conduction)
            HLLCHC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 33:  // RoeHC (Roe with Heat Conduction)
            RoeHC_Flux(space_dir, var, rows, cols, GC, Conserl, Conserr, Flux);
            break;
        case 5:   // RS_Marquina (not implemented in current code)
            // RS_Marquina(3, rows, y, Flux, dt, dx);
            break;
        default:
            // ER_Flux(var, rows, GC, prileft, priright, Flux);
            // Add appropriate handling logic based on actual requirements
            break;
    }

    // Store y-direction fluxes in g array (for interior cells including one ghost cell)
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < var; k++) {
        for (int i = GC - 1; i <= rows - GC; i++) {
            for (int j = GC - 1; j <= cols - GC; j++) {
                g[k][i][j] = Flux[k][i][j];
            }
        }
    }

    // Free allocated memory
    free(Conserl);
    free(Conserr);
    free(Flux);
}


/*                                      ******************                                          */
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */

                                    /*……………………………………………………*/
                                /*具有人工热传导的近似黎曼求解数值方法*/
                                    /*……………………………………………………*/
/*……………………………………………………………………………………………………*/


static inline void Flux_Reconstruction_RP_Heat(int AR_scheme, int var, int rows, int GC, double (*y)[rows],double (*z)[rows]\
                                ,double dt,double dx, double u_Refer) {
    int i,j;
    double Conserl[var][rows],Conserr[var][rows];
    double prileft[var][rows], priright[var][rows];
    double Flux[var][rows];
    
     switch (Recon_Accur){
        case 1:
//            Reconstruction_Godunov(var,rows,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
//            TVD_Reconstruction(var,rows,GC,y,Conserl,Conserr,dx);
            break;
        case 3:
//            WENO3_Reconstruction(var,rows,GC,y,Conserl,Conserr);
            break;
        case 5:
//            WENO5_Reconstruction(var,rows,GC,y,Conserl,Conserr);
//            WENO5_Reconstruction_C(var,rows,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }
    
//    Con_to_Pri_1D(3,rows,prileft,Conserl);
//   Con_to_Pri_1D(3,rows,priright,Conserr);


    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 2:
            HLL_Flux_HeatConduction(var, rows, GC, prileft,priright,Flux,u_Refer);
            break;
        case 3:
            HLLC_Flux_HeatConduction(var, rows, GC, prileft,priright,Flux,u_Refer);
            break;
        case 4:
            Roe_Flux_HeatConduction(var, rows, GC, prileft,priright,Flux,u_Refer);
            break;
        case 5:
            RS_Marquina(3,rows,y,z,dt,dx);
            break;
        default:
            ER_Flux_PlusHeat(var, rows, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }


    for ( i = 0; i < var; i++){
        for ( j = 1; j < rows-1; j++){
            z[i][j] = Flux[i][j];
        }
    }

}




#endif 
