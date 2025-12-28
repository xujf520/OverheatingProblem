#ifndef TIME_ADVANCE_H
#define TIME_ADVANCE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Golbal.h"
#include "mesh.h"
#include "initialize.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "scheme.h"



int Time_Step = 0;
//static double mass_error = 0.0;
//static double momentum_error = 0.0;
//static double energy_error = 0.0;

void OutputFluxData_file();
void OutputConservationErrors_file(); 


// Unified Rk1 Time Advancement Function
/**
 * First-order Runge-Kutta (Forward Euler) time advancement with multi-block decomposition
 * Unified function that handles both single-block and multi-block cases
 * 
 * @param AR_scheme Approximate Riemann solver scheme identifier
 * @param var Number of variables (typically 4 for Euler equations)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells
 * @param y Solution array in conservative variables [var][rows][cols]
 * @param dt Time step
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @param test_case Test case identifier from TestCase2D enumeration
 */

static inline void RK1_TimeAd_Unified(int AR_scheme, int var, int rows, int cols, int GC, 
                                      double (*y)[rows][cols], double dt, double dx, double dy, 
                                      int test_case) {
    int i, j, k;
    
    // Calculate number of interior cells
    int LNX = rows - 2 * GC;
    int LNY = cols - 2 * GC;

    // Set up block decomposition based on test case
    BlockGeometricParameters BGP = setup_block_decomposition(test_case, LNX, LNY, GC);
    
    // Get overlap region information for data exchange between blocks
    DataExchangeArea DEA = Get_OverLap(test_case, LNX, LNY, GC);

    // Allocate memory for overlap region data exchange
    double (*Lap_U)[DEA.OverLap_X[0]][DEA.OverLap_Y[0]] = 
        malloc(var * sizeof(double[DEA.OverLap_X[0]][DEA.OverLap_Y[0]]));

    // Check memory allocation for overlap region
    if (Lap_U == NULL) {
        fprintf(stderr, "Memory allocation failed in RK1_TimeAd_Unified Point1\n");
        free(Lap_U); 
        return;
    }

    // Initialize overlap region array to zero
    for (k = 0; k < var; k++) {
        for (i = 0; i < DEA.OverLap_X[0]; i++) {
            for (j = 0; j < DEA.OverLap_Y[0]; j++) {
                Lap_U[k][i][j] = 0.0;
            }
        }
    }

    // Copy data from global array to overlap region
    // This extracts the region where blocks will exchange data
    for (k = 0; k < var; k++) {
        for (i = 0; i < DEA.OverLap_X[0]; i++) {
            for (j = 0; j < DEA.OverLap_Y[0]; j++) {
                Lap_U[k][i][j] = y[k][i + DEA.start_OLX[0]][j + DEA.start_OLY[0]];
            }
        }
    }

    // Process each block in the domain decomposition
    for (int Blockn = 0; Blockn < BGP.num_blocks; Blockn++) {
        // Allocate memory for each block's arrays
        double (*BlockU)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockSI)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockF)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockG)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockS)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));

        // Check memory allocation for block arrays
        if (BlockU == NULL || BlockSI == NULL || BlockF == NULL || BlockG == NULL || BlockS == NULL) {
            fprintf(stderr, "Memory allocation failed in RK1_TimeAd_Unified (backward step path)\n");
            free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
            return;
        }

        // Initialize block arrays to zero
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                    BlockU[k][i][j] = 0.0;
                    BlockSI[k][i][j] = 0.0;
                    BlockF[k][i][j] = 0.0;
                    BlockG[k][i][j] = 0.0;
                    BlockS[k][i][j] = 0.0;
                }
            }
        }

        // Copy relevant portion from global array to block array
        // Only copy interior cells (excluding ghost cells in the copy operation)
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    BlockU[k][i][j] = y[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]];
                }
            }
        }

        // Calculate source terms (e.g., gravity)
        if (Source) {
            Source_Gravity(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, BlockU, BlockS); 
        }
        
        // Apply boundary conditions
        if (BGP.num_blocks == 1) {
            // Single block: use standard boundary conditions
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
        }
        else {
            // Multiple blocks: get block-specific boundary conditions
            Get_block_BC(test_case, Blockn);
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
            
            // Exchange data between blocks through overlap region
            Block_DataExchange(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, 
                              DEA.OverLap_X[0], DEA.OverLap_Y[0], Lap_U, test_case, Blockn);
        }
    
//        double (*Pri)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
//            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
//        Con_to_Pri_2D(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], Pri, BlockU);
//        free(Pri);

        // Calculate fluxes and spatial discretization terms for the block
        Flux_Reconstruction_RP(AR_scheme, var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, 
                               BlockU, BlockF, BlockG, dt, dx, dy);
        Space_Discrete_Item(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, dx, dy, 
                           BlockF, BlockG, BlockS, BlockSI);

        // Update block solution using forward Euler method
        #pragma omp parallel for collapse(3)
        for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
            for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                for (k = 0; k < var; k++) {
                    BlockU[k][i][j] = BlockU[k][i][j] + dt * BlockSI[k][i][j];
                }
            }
        }

        // Copy updated block data back to global array
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    y[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] = BlockU[k][i][j];
                }
            }
        }
        
        // Free memory allocated for this block
        free(BlockU); 
        free(BlockSI); 
        free(BlockF); 
        free(BlockG); 
        free(BlockS);
    }

    // Clean up: free decomposition structures
    free_block_decomposition(&BGP);
    free(Lap_U);  
}



/*static inline void RK3_TimeAd(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt,double dx, double dy) {
    int i,j,k;

    double (*Conser_U1)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U2)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_F)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Source_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    // 检查内存分配是否成功
    if (Flux_F == NULL || Flux_G == NULL) {
        fprintf(stderr, "Memory allocation failed in RK3_TimeAd\n");
        // 释放已分配的内存
        free(Conser_U1);
        free(Conser_U2);
        free(Flux_F);
        free(Flux_G);
        free(Source_G);
        return;
    }    
    //初始化参数
    for ( k = 0; k < var; k++){
        for ( i = 0; i < rows; i++){
            for ( j = 0; j < cols; j++){
                Conser_U1[k][i][j] = Conser_U2[k][i][j] = 0.0;
                Flux_F[k][i][j] = Flux_G[k][i][j] = 0.0;
                Source_G[k][i][j] = 0.0;
            }
        }
    }

    
    // 第一步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,y,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, y, GC,Time); 
    Boundary_Conditions(var, rows, cols, y, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,y,Flux_F,Flux_G,dt,dx,dy);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                Conser_U1[k][i][j] = y[k][i][j] - dt * (Flux_F[i][j][k]-Flux_F[i][j-1][k])/dx \
                                            - dt * (Flux_G[i][j][k]-Flux_G[i][j][k-1])/dy + dt * Source_G[k][i][j]; 


    //第二步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U1,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, Conser_U1, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U1, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U1,Flux_F,Flux_G,dt,dx,dy);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                Conser_U2[k][i][j] = (3.0/4.0) * y[k][i][j] +  (1.0/4.0) * (Conser_U1[i][j][k] \
                                            - dt * (Flux_F[i][j][k]-Flux_F[i][j-1][k])/dx \
                                            - dt * (Flux_G[i][j][k]-Flux_G[i][j][k-1])/dy + dt * Source_G[k][i][j]);
    
    
    //第三步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U2,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, Conser_U2, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U2, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U2,Flux_F,Flux_G,dt,dx,dy);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                y[k][i][j] = (1.0/2.0) * y[k][i][j] +  (2.0/3.0) * (Conser_U2[i][j][k] \
                                                    - dt * (Flux_F[i][j][k]-Flux_F[i][j-1][k])/dx \
                                                    - dt * (Flux_G[i][j][k]-Flux_G[i][j][k-1])/dy + dt * Source_G[k][i][j]);



    free(Conser_U1);
    free(Conser_U2);
    free(Flux_F);
    free(Flux_G);
    free(Source_G); 

}*/


static inline void RK3_TimeAd_Unified(int AR_scheme, int var, int rows, int cols, int GC,
                                     double (*y)[rows][cols], double dt, double dx, double dy,
                                     int test_case) {
    int i, j, k;
    int LNX = rows - 2 * GC;
    int LNY = cols - 2 * GC;

    // Allocate global arrays for RK3 intermediate stages
    double (*Conser_U1)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U2)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    
    if (Conser_U1 == NULL || Conser_U2 == NULL) {
        fprintf(stderr, "Memory allocation failed for RK3 stage arrays\n");
        free(Conser_U1);
        free(Conser_U2);
        return;
    }
    
    // Initialize RK stage arrays
    #pragma omp parallel for collapse(3)
    for (k = 0; k < var; k++) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                Conser_U1[k][i][j] = 0.0;
                Conser_U2[k][i][j] = 0.0;
            }
        }
    }

    // Setup domain decomposition and overlap information
    BlockGeometricParameters BGP = setup_block_decomposition(test_case, LNX, LNY, GC);
    DataExchangeArea DEA = Get_OverLap(test_case, LNX, LNY, GC);

    // Allocate single overlap buffer (reused for all stages)
    double (*Lap_U)[DEA.OverLap_X[0]][DEA.OverLap_Y[0]] = malloc(var * sizeof(double[DEA.OverLap_X[0]][DEA.OverLap_Y[0]]));
    
    if (Lap_U == NULL) {
        fprintf(stderr, "Memory allocation failed for overlap region buffer\n");
        free(Lap_U);
        free_block_decomposition(&BGP);
        free(Conser_U1);
        free(Conser_U2);
        return;
    }

    // ==================== RK3 Stage 1 ====================
    // Initialize overlap buffer from global solution array
    #pragma omp parallel for collapse(3)
    for (k = 0; k < var; k++) {
        for (i = 0; i < DEA.OverLap_X[0]; i++) {
            for (j = 0; j < DEA.OverLap_Y[0]; j++) {
               Lap_U[k][i][j] = y[k][i + DEA.start_OLX[0]][j + DEA.start_OLY[0]];
            }
        }
    }

    // Process each block for stage 1
    for (int Blockn = 0; Blockn < BGP.num_blocks; Blockn++) {
        // Allocate block arrays (reused across RK stages to reduce allocation overhead)
        double (*BlockU)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockSI)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockF)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockG)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockS)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));

        if (BlockU == NULL || BlockSI == NULL || BlockF == NULL || BlockG == NULL || BlockS == NULL) {
            fprintf(stderr, "Memory allocation failed for block arrays in stage 1\n");
            free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
            free(Lap_U); free_block_decomposition(&BGP); free(Conser_U1); free(Conser_U2);
            return;
        }

        // Initialize block arrays (zero out only necessary arrays)
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                    BlockU[k][i][j] = 0.0;
                    BlockSI[k][i][j] = 0.0;
                    BlockF[k][i][j] = 0.0;
                    BlockG[k][i][j] = 0.0;
                }
            }
        }

        // Copy interior region from global array to block array
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    BlockU[k][i][j] = y[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]];
                }
            }
        }

        // Calculate source terms if enabled
        if (Source) {
            Source_Gravity(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, BlockU, BlockS);
        } else {
            // Zero out source array if not used
            #pragma omp parallel for collapse(3)
            for (k = 0; k < var; k++) {
                for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                    for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                        BlockS[k][i][j] = 0.0;
                    }
                }
            }
        }
        
        // Apply boundary conditions
        if (BGP.num_blocks == 1) {
            // Single block: standard boundary conditions
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
        } else {
            // Multiple blocks: block-specific boundary conditions
            Get_block_BC(test_case, Blockn);
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
            Block_DataExchange(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU,
                             DEA.OverLap_X[0], DEA.OverLap_Y[0], Lap_U, test_case, Blockn);
        }
    
        // Calculate fluxes and spatial discretization terms
        Flux_Reconstruction_RP(AR_scheme, var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC,
                               BlockU, BlockF, BlockG, dt, dx, dy);
        Space_Discrete_Item(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, dx, dy,
                           BlockF, BlockG, BlockS, BlockSI);

        // Stage 1 update: U1 = y + dt * SI
        #pragma omp parallel for collapse(3)
        for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
            for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                for (k = 0; k < var; k++) {
                    BlockU[k][i][j] = BlockU[k][i][j] + dt * BlockSI[k][i][j];
                }
            }
        }

        // Copy updated block data to stage 1 global array
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    Conser_U1[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] = BlockU[k][i][j];
                }
            }
        }

        // Free block arrays after stage 1
        free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
    }

    // ==================== RK3 Stage 2 ====================
    // Update overlap buffer from stage 1 solution
    #pragma omp parallel for collapse(3)
    for (k = 0; k < var; k++) {
        for (i = 0; i < DEA.OverLap_X[0]; i++) {
            for (j = 0; j < DEA.OverLap_Y[0]; j++) {
               Lap_U[k][i][j] = Conser_U1[k][i + DEA.start_OLX[0]][j + DEA.start_OLY[0]];
            }
        }
    }

    // Process each block for stage 2
    for (int Blockn = 0; Blockn < BGP.num_blocks; Blockn++) {
        // Reallocate block arrays for stage 2
        double (*BlockU)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockSI)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockF)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockG)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockS)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));

        if (BlockU == NULL || BlockSI == NULL || BlockF == NULL || BlockG == NULL || BlockS == NULL) {
            fprintf(stderr, "Memory allocation failed for block arrays in stage 2\n");
            free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
            free(Lap_U); free_block_decomposition(&BGP); free(Conser_U1); free(Conser_U2);
            return;
        }

        // Initialize block arrays
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                    BlockU[k][i][j] = 0.0;
                    BlockSI[k][i][j] = 0.0;
                    BlockF[k][i][j] = 0.0;
                    BlockG[k][i][j] = 0.0;
                }
            }
        }

        // Copy interior region from stage 1 solution
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    BlockU[k][i][j] = Conser_U1[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]];
                }
            }
        }

        // Calculate source terms if enabled
        if (Source) {
            Source_Gravity(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, BlockU, BlockS);
        } else {
            #pragma omp parallel for collapse(3)
            for (k = 0; k < var; k++) {
                for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                    for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                        BlockS[k][i][j] = 0.0;
                    }
                }
            }
        }
        
        // Apply boundary conditions
        if (BGP.num_blocks == 1) {
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
        } else {
            Get_block_BC(test_case, Blockn);
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
            Block_DataExchange(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU,
                              DEA.OverLap_X[0], DEA.OverLap_Y[0], Lap_U, test_case, Blockn);
        }
    
        // Calculate fluxes and spatial discretization terms
        Flux_Reconstruction_RP(AR_scheme, var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC,
                               BlockU, BlockF, BlockG, dt, dx, dy);
        Space_Discrete_Item(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, dx, dy,
                           BlockF, BlockG, BlockS, BlockSI);

        // Stage 2 update: U2 = (3/4)*y + (1/4)*(U1 + dt*SI)
        #pragma omp parallel for collapse(3)
        for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
            for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                for (k = 0; k < var; k++) {
                    BlockU[k][i][j] = (3.0/4.0) * y[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] + 
                                      (1.0/4.0) * (Conser_U1[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] + dt * BlockSI[k][i][j]);
                }
            }
        }

        // Copy updated block data to stage 2 global array
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    Conser_U2[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] = BlockU[k][i][j];
                }
            }
        }

        // Free block arrays after stage 2
        free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
    }

    // ==================== RK3 Stage 3 ====================
    // Update overlap buffer from stage 2 solution
    #pragma omp parallel for collapse(3)
    for (k = 0; k < var; k++) {
        for (i = 0; i < DEA.OverLap_X[0]; i++) {
            for (j = 0; j < DEA.OverLap_Y[0]; j++) {
               Lap_U[k][i][j] = Conser_U2[k][i + DEA.start_OLX[0]][j + DEA.start_OLY[0]];
            }
        }
    }

    // Process each block for stage 3
    for (int Blockn = 0; Blockn < BGP.num_blocks; Blockn++) {
        // Reallocate block arrays for stage 3
        double (*BlockU)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockSI)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockF)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockG)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));
        double (*BlockS)[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]] = 
            malloc(var * sizeof(double[BGP.LNX_NGC[Blockn]][BGP.LNY_NGC[Blockn]]));

        if (BlockU == NULL || BlockSI == NULL || BlockF == NULL || BlockG == NULL || BlockS == NULL) {
            fprintf(stderr, "Memory allocation failed for block arrays in stage 3\n");
            free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
            //free(Lap_U); free_block_decomposition(&BGP); free(Conser_U1); free(Conser_U2);
            return;
        }

        // Initialize block arrays
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                    BlockU[k][i][j] = 0.0;
                    BlockSI[k][i][j] = 0.0;
                    BlockF[k][i][j] = 0.0;
                    BlockG[k][i][j] = 0.0;
                }
            }
        }

        // Copy interior region from stage 2 solution
        #pragma omp parallel for collapse(3)
        for (k = 0; k < var; k++) {
            for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
                for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                    BlockU[k][i][j] = Conser_U2[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]];
                }
            }
        }

        // Calculate source terms if enabled
        if (Source) {
            Source_Gravity(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, BlockU, BlockS);
        } else {
            #pragma omp parallel for collapse(3)
            for (k = 0; k < var; k++) {
                for (i = 0; i < BGP.LNX_NGC[Blockn]; i++) {
                    for (j = 0; j < BGP.LNY_NGC[Blockn]; j++) {
                        BlockS[k][i][j] = 0.0;
                    }
                }
            }
        }
        
        // Apply boundary conditions
        if (BGP.num_blocks == 1) {
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
        } else {
            Get_block_BC(test_case, Blockn);
            Boundary_Conditions(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU, GC);
            Block_DataExchange(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], BlockU,
                              DEA.OverLap_X[0], DEA.OverLap_Y[0], Lap_U, test_case, Blockn);
        }
    
        // Calculate fluxes and spatial discretization terms
        Flux_Reconstruction_RP(AR_scheme, var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC,
                               BlockU, BlockF, BlockG, dt, dx, dy);
        Space_Discrete_Item(var, BGP.LNX_NGC[Blockn], BGP.LNY_NGC[Blockn], GC, dx, dy,
                           BlockF, BlockG, BlockS, BlockSI);

        // Stage 3 update: y = (1/3)*y + (2/3)*(U2 + dt*SI)
        #pragma omp parallel for collapse(3)
        for (i = GC; i < BGP.LNX[Blockn] + GC; i++) {
            for (j = GC; j < BGP.LNY[Blockn] + GC; j++) {
                for (k = 0; k < var; k++) {
                    y[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] = 
                        (1.0/3.0) * y[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] + 
                        (2.0/3.0) * (Conser_U2[k][i + BGP.start_i[Blockn]][j + BGP.start_j[Blockn]] + dt * BlockSI[k][i][j]);
                }
            }
        }

        // Free block arrays after stage 3
        free(BlockU); free(BlockSI); free(BlockF); free(BlockG); free(BlockS);
    }

    // Clean up allocated memory
    free_block_decomposition(&BGP);
    free(Lap_U);
    free(Conser_U1);
    free(Conser_U2);
}



static inline void RK3_TimeAd(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt,double dx, double dy) {
    int i,j,k;

    double (*Space_Item)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U1)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U2)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_F)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Source_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    // 检查内存分配是否成功
    if (Flux_F == NULL || Flux_G == NULL) {
        fprintf(stderr, "Memory allocation failed in RK3_TimeAd\n");
        // 释放已分配的内存
        free(Space_Item);
        free(Conser_U1);
        free(Conser_U2);
        free(Flux_F);
        free(Flux_G);
        free(Source_G);
        return;
    }    
    //初始化参数

    #pragma omp parallel for collapse(3)
    for ( k = 0; k < var; k++){
        for ( i = 0; i < rows; i++){
            for ( j = 0; j < cols; j++){
                Conser_U1[k][i][j] = Conser_U2[k][i][j] = 0.0;
                Flux_F[k][i][j] = Flux_G[k][i][j] = 0.0;
                Space_Item[k][i][j] = 0.0;
                Source_G[k][i][j] = 0.0;
            }
        }
    }

    // 第一步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,y,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, y, GC,Time); 
    Boundary_Conditions(var, rows, cols, y, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,y,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    #pragma omp parallel for collapse(3)
    for ( k = 0; k < var; k++){
        for ( i = GC; i <= rows-GC-1; i++){
            for ( j = GC; j <= cols-GC-1; j++){
                Conser_U1[k][i][j] = y[k][i][j] + dt * Space_Item[k][i][j];
            }
        }
    }

    //第二步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U1,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, Conser_U1, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U1, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U1,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    #pragma omp parallel for collapse(3)
    for ( k = 0; k < var; k++){
        for ( i = GC; i <= rows-GC-1; i++){
            for ( j = GC; j <= cols-GC-1; j++){
                Conser_U2[k][i][j] = (3.0/4.0) * y[k][i][j] +  (1.0/4.0) * (Conser_U1[k][i][j] + dt * Space_Item[k][i][j]);
            }
        }
    }
    
    //第三步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U2,Source_G); 
    // 施加边界条件
//    BC_DoubleMach_2D(var, rows, cols, Conser_U2, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U2, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U2,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    for ( k = 0; k < var; k++){
        for ( i = GC; i <= rows-GC-1; i++){
            for ( j = GC; j <= cols-GC-1; j++){
                y[k][i][j] = (1.0/3.0) * y[k][i][j] +  (2.0/3.0) * (Conser_U2[k][i][j] + dt * Space_Item[k][i][j]);
            }
        }
    }


    free(Space_Item);
    free(Conser_U1);
    free(Conser_U2);
    free(Flux_F);
    free(Flux_G);
    free(Source_G); 

}


                           /*First Order Euler Time Stepping Method*/
static inline void RK1_TimeAd(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt, double dx, double dy) {
    int i,j,k;
    double Conser_1[3]={0.0}, Conser_2[3]={0.0};
    double (*Space_Item)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_F)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Source_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));

    // 检查内存分配是否成功
    if (Flux_F == NULL || Flux_G == NULL) {
        fprintf(stderr, "Memory allocation failed in RK1_TimeAd\n");
        // 释放已分配的内存
        free(Space_Item);
        free(Flux_F);
        free(Flux_G);
        free(Source_G);
        return;
    }

  
    //初始化参数
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < var; k++) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                Space_Item[k][i][j] = 0.0;
                Flux_F[k][i][j] = 0.0;
                Flux_G[k][i][j] = 0.0;
                Source_G[k][i][j] = 0.0;
            }
        }
    }

    if (Source)
        Source_Gravity(var,rows,cols,GC,y,Source_G); 

   // 第一步计算
    // 施加边界条件
    Boundary_Conditions(var, rows, cols, y, GC); 
//    BC_BackwardStep_2D(var, rows, cols, y, GC,Time); 
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,y,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);


    #pragma omp parallel for collapse(3)
    for (int i = GC; i <= rows - GC - 1; i++) {
        for (int j = GC; j <= cols - GC - 1; j++) {
            for (int k = 0; k < var; k++) {
                y[k][i][j] = y[k][i][j] + dt * Space_Item[k][i][j];
            }
        }
    }
            
        
    free(Space_Item);
    free(Flux_F);
    free(Flux_G);
    free(Source_G);
}


static inline void RK1_TimeAd_BS(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt, double dx, double dy) {
    int i,j,k;
    double Conser_1[3]={0.0}, Conser_2[3]={0.0};


    //划分Block
    
    int LNX = rows - 2 * GhostCell;
    int LNY = cols - 2 * GhostCell;
    int B1_LNX = LNX;
    int B1_LNY = LNY * 0.8;
    int B2_LNX = LNX * 0.2;
    int B2_LNY = LNY * 0.2; 
    int B1_LNX_ngc = LNX + 2 * GhostCell;
    int B1_LNY_ngc = LNY * 0.8 + 2 * GhostCell;
    int B2_LNX_ngc = LNX * 0.2 + 2 * GhostCell;
    int B2_LNY_ngc = LNY * 0.2 + 2 * GhostCell; 

    double (*B1U)[B1_LNX_ngc][B1_LNY_ngc] = malloc(4 * sizeof(double[B1_LNX_ngc][B1_LNY_ngc]));
    double (*B1SI)[B1_LNX_ngc][B1_LNY_ngc] = malloc(4 * sizeof(double[B1_LNX_ngc][B1_LNY_ngc]));
    double (*B1FU)[B1_LNX_ngc][B1_LNY_ngc] = malloc(4 * sizeof(double[B1_LNX_ngc][B1_LNY_ngc]));
    double (*B1GU)[B1_LNX_ngc][B1_LNY_ngc] = malloc(4 * sizeof(double[B1_LNX_ngc][B1_LNY_ngc]));
    double (*B1S)[B1_LNX_ngc][B1_LNY_ngc] = malloc(4 * sizeof(double[B1_LNX_ngc][B1_LNY_ngc]));

    double (*B2SI)[B2_LNX_ngc][B2_LNY_ngc] = malloc(4 * sizeof(double[B2_LNX_ngc][B2_LNY_ngc]));
    double (*B2U)[B2_LNX_ngc][B2_LNY_ngc] = malloc(4 * sizeof(double[B2_LNX_ngc][B2_LNY_ngc]));
    double (*B2FU)[B2_LNX_ngc][B2_LNY_ngc] = malloc(4 * sizeof(double[B2_LNX_ngc][B2_LNY_ngc]));
    double (*B2GU)[B2_LNX_ngc][B2_LNY_ngc] = malloc(4 * sizeof(double[B2_LNX_ngc][B2_LNY_ngc]));
    double (*B2S)[B2_LNX_ngc][B2_LNY_ngc] = malloc(4 * sizeof(double[B2_LNX_ngc][B2_LNY_ngc]));


    // 检查内存分配是否成功
    if (B1U == NULL || B1SI == NULL|| B1FU == NULL || B1GU == NULL || B1S == NULL || B2U == NULL || B2SI == NULL|| B2FU == NULL || B2GU == NULL || B2S == NULL) {
        fprintf(stderr, "Memory allocation failed in RK1_TimeAd\n");
        // 释放已分配的内存
        free(B1U), free(B1SI), free(B1FU), free(B1GU), free(B1S);
        free(B2U), free(B2SI), free(B2FU), free(B2GU), free(B2S);
        return;
    }

    //初始化Block参数
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < var; k++) {
        for (int i = 0; i < B1_LNX_ngc; i++) {
            for (int j = 0; j < B1_LNY_ngc; j++) {
                B1U[k][i][j] = 0.0;
                B1SI[k][i][j] = 0.0;
                B1FU[k][i][j] = 0.0;
                B1GU[k][i][j] = 0.0;
                B1S[k][i][j] = 0.0;
            }
        }
    }
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < var; k++) {
        for (int i = 0; i < B2_LNX_ngc; i++) {
            for (int j = 0; j < B2_LNY_ngc; j++) {
                B2U[k][i][j] = 0.0;
                B2SI[k][i][j] = 0.0;
                B2FU[k][i][j] = 0.0;
                B2GU[k][i][j] = 0.0;
                B2S[k][i][j] = 0.0;
            }
        }
    }

    //初始化每个Blcok的值
    for (int k = 0; k < var; k++) {
        for (int i = GC; i < B1_LNX_ngc-GC; i++) {
            for (int j = GC; j < B1_LNY_ngc-GC; j++) {
                B1U[k][i][j] = y[k][i][j + B2_LNY];
            }
        }
    }

    for (int k = 0; k < var; k++) {
        for (int i = GC; i < B2_LNX_ngc-GC; i++) {
            for (int j = GC; j < B2_LNY_ngc-GC; j++) {
                B2U[k][i][j] = y[k][i][j];
            }
        }
    }

    if (Source){
        Source_Gravity(var,B1_LNX_ngc,B1_LNY_ngc,GC,B1U,B1S); 
        Source_Gravity(var,B2_LNX_ngc,B2_LNY_ngc,GC,B2U,B2S); 
    }
        
   // 第一步计算
    // 施加边界条件

    BC_BackwardStep_2D(var,B1_LNX_ngc,B1_LNY_ngc,B1U,B2_LNX_ngc,B2_LNY_ngc,B2U); 

    Flux_Reconstruction_RP(AR_scheme,var,B1_LNX_ngc,B1_LNY_ngc,GC,B1U,B1FU,B1GU,dt,dx,dy);
    Space_Discrete_Item(var,B1_LNX_ngc,B1_LNY_ngc,GC,dx,dy,B1FU,B1GU,B1S,B1SI);

    #pragma omp parallel for collapse(3)
    for (int i = GC; i <= B1_LNX_ngc - GC - 1; i++) {
        for (int j = GC; j <= B1_LNY_ngc - GC - 1; j++) {
            for (int k = 0; k < var; k++) {
                B1U[k][i][j] = B1U[k][i][j] + dt * B1SI[k][i][j];
            }
        }
    }

    Flux_Reconstruction_RP(AR_scheme,var,B2_LNX_ngc,B2_LNY_ngc,GC,B2U,B2FU,B2GU,dt,dx,dy);
    Space_Discrete_Item(var,B2_LNX_ngc,B2_LNY_ngc,GC,dx,dy,B2FU,B2GU,B2S,B2SI);

    #pragma omp parallel for collapse(3)
    for (int i = GC; i <= B2_LNX_ngc - GC - 1; i++) {
        for (int j = GC; j <= B2_LNY_ngc - GC - 1; j++) {
            for (int k = 0; k < var; k++) {
                B2U[k][i][j] = B2U[k][i][j] + dt * B2SI[k][i][j];
            }
        }
    }

    //更新每个Blcok的值
    for (int k = 0; k < var; k++) {
        for (int i = GC; i < B1_LNX_ngc-GC; i++) {
            for (int j = GC; j < B1_LNY_ngc-GC; j++) {
                y[k][i][j + B2_LNY] = B1U[k][i][j];
            }
        }
    }

    for (int k = 0; k < var; k++) {
        for (int i = GC; i < B2_LNX_ngc-GC; i++) {
            for (int j = GC; j < B2_LNY_ngc-GC; j++) {
                y[k][i][j] = B2U[k][i][j];
            }
        }
    }
     
    
    free(B1U), free(B1SI), free(B1FU), free(B1GU), free(B1S);
    free(B2U), free(B2SI), free(B2FU), free(B2GU), free(B2S);
        
}



void OutputConservationErrors_file(int time_step, double Conser_1[3], double Conser_2[3]) {
    // 创建结果目录（如果不存在）
    system("mkdir -p /mnt/d/Desktop/RP_FVM/data");
    
    // 打开文件（追加模式）
    FILE* file = fopen("/mnt/d/Desktop/RP_FVM/data/Conser_Error.dat", "a");  // 修正文件名拼写
    if (file == NULL) {
        printf("无法打开文件 /mnt/d/Desktop/RP_FVM/data/Conser_Error.dat\n");
        return;
    }
    
    // 计算误差（使用更高精度的中间变量）
    volatile double mass_error = Conser_2[0] - Conser_1[0];      // volatile防止编译器优化
    volatile double momentum_error = Conser_2[1] - Conser_1[1];
    volatile double energy_error = Conser_2[2] - Conser_1[2];
    
    // 如果是第一个时间步，写入文件头部
    if (time_step == 0) {
        fprintf(file, "Variables = \t TimeStep \t  Mass_Error  \t Momentum_Error  \t Energy_Error \n");
        fprintf(file, "Zone T=\"Conservation Errors\"\n");
    }
    
    // 写入数据（强制16位小数，使用科学计数法）
    fprintf(file, "%d\t%.16e\t%.16e\t%.16e\n", 
            time_step, 
            (double)mass_error, 
            (double)momentum_error, 
            (double)energy_error);
    
    // 立即刷新缓冲区确保数据写入
    fflush(file);
    
    // 关闭文件
    fclose(file);
    
}





void OutputFluxData_file(int k, int rows, double dx, double* x,double (*z)[rows]) {
    char filename[100];
    sprintf(filename, "D:/Desktop/Single_Med_data/flux_output_%d.dat", k);
    // 打开文件
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("无法打开文件 %s\n", filename);
        return; // 返回错误代码
    }
    // 写入文件头部信息
    fprintf(file, "variables =  'x' 'flux_rho' 'flux_rhou' 'flux_E' \n");
    // 写入计算结果
    for (int i = 2; i < rows-2 ; i++) {
        fprintf(file,"%.8f\t%.8f\t%.8f\t%.8f\n",x[i-2]+0.5*dx, z[0][i],z[1][i], z[2][i]);
    }
    fclose(file);
    printf("numerical flux ：The %d th calculation ended\n", k);
}

#endif 
