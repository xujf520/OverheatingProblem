#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "function.h"
#include "Characteriz.h"
#include "Golbal.h"

static inline double TVD_minmod_L();
static inline double TVD_minmod_R();
static inline double TVD_vanleer_L();
static inline double TVD_vanleer_R();
static inline double TVD_VanAlbada_L();
static inline double TVD_VanAlbada_R();
static inline double OED_TVD_VanAlbada_L();
static inline double OED_TVD_VanAlbada_R();
static inline double WENO3_L();
static inline double WENO3_R();
static inline double WENO5_L();
static inline double WENO5_R();
static inline double WENO5Z_L();
static inline double WENO5Z_R();

/*                               ************************************                               */
/*                               ************************************                               */
/*                                    Reconstruction Scheme                                         */
/*                               ************************************                               */
/*                               ************************************                               */


/*                                      ******************                                          */
/*                                            Godunov                                               */
/*                                      ******************                                          */
/**
 * First-order Godunov reconstruction (piecewise constant)
 * Reconstructs left and right states at cell interfaces
 * Can perform reconstruction in characteristic variables for better accuracy
 * 
 * @param dir Direction of reconstruction: 1 for x-direction, 2 for y-direction
 * @param var Number of variables (typically 4 for Euler equations)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells
 * @param y Input conservative variables array [var][rows][cols]
 * @param conserl Output: left reconstructed states at cell interfaces
 * @param conserr Output: right reconstructed states at cell interfaces
 * @param delta_x Grid spacing (used for characteristic decomposition but not for 1st order)
 */
static inline void Reconstruction_Godunov(int dir, int var, int rows, int cols, int GC, 
                                          double (*y)[rows][cols],
                                          double (*conserl)[rows][cols], 
                                          double (*conserr)[rows][cols], 
                                          double delta_x) {
    double epsilo = 1e-6;  // Small constant (not used in Godunov but kept for compatibility)
    
    // Check if characteristic decomposition is enabled
    if (Characteriz) {
        // Dynamic memory allocation for characteristic decomposition arrays
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));        // Primitive variables
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Left eigenvectors
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Right eigenvectors
        
        // Check memory allocation
        if (!Pri || !Eigen_L || !Eigen_R) {
            fprintf(stderr, "Memory allocation failed in Reconstruction_Godunov\n");
            free(Pri); free(Eigen_L); free(Eigen_R);
            return;
        }

        // Convert conservative variables to primitive variables
        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            // X-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=1.0, ny=0.0)
            Compute_Eigen_2DX(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in x-direction
            #pragma omp parallel for
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    // Local arrays for characteristic decomposition
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Single characteristic variable (left)
                    double Characteristic_Variable_R;  // Single characteristic variable (right)
                    
                    // Project conservative variables onto characteristic space
                    // Mathematical operation: w = L * u, where:
                    // w = characteristic variables (wave strengths)
                    // L = left eigenvector matrix (size: 4×4)
                    // u = conservative variables (size: 4×1)
                    for (int k = 0; k < var; k++) {
                        Characteristic_Variable_L = 0.0;
                        Characteristic_Variable_R = 0.0;
                        
                        // Matrix-vector multiplication: w_k = Σ_{m=0}^3 (L_km * u_m)
                        // For each characteristic field k, compute weighted sum of all conservative variables
                        for (int m = 0; m < var; m++) {
                            Characteristic_Variable_L += y[m][i][j] * Eigen_L[k][m][i][j];
                            Characteristic_Variable_R += y[m][i + 1][j] * Eigen_L[k][m][i][j];
                        }
                        
                        // Reconstruct conservative variables from characteristic space
                        // Mathematical operation: u' = R * w, where:
                        // u' = reconstructed conservative variables
                        // R = right eigenvector matrix (size: 4×4)
                        // Note: R[m][k] accesses element at row m, column k of eigenvector matrix
                        for (int m = 0; m < var; m++) {
                            // Chara_L[k][m] = contribution of k-th characteristic variable to m-th conservative variable
                            // = R_mk * w_k, where w_k = Characteristic_Variable_L
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Reconstruct states at interface (Godunov: simple average in characteristic space)
                    // For each conservative variable, sum contributions from all characteristic fields:
                    // conserl_k = Σ_{m=0}^3 (Chara_L[m][k]) = Σ_{m=0}^3 (R_km * w_m_L)
                    // conserr_k = Σ_{m=0}^3 (Chara_R[m][k]) = Σ_{m=0}^3 (R_km * w_m_R)
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                        
                        // Complete matrix multiplication: u' = Σ (characteristic contributions)
                        // This implements the full transformation: u' = R * w = R * (L * u) = I * u = u
                        // In theory, R * L = I (identity matrix), so this should be an identity transform
                        for (int m = 0; m < var; m++) {
                            // Index mapping:
                            // k: conservative variable index (0=density, 1=x-momentum, 2=y-momentum, 3=energy)
                            // m: characteristic field index (0=entropy wave, 1=acoustic wave(-), 2=shear wave, 3=acoustic wave(+))
                            conserl[k][i][j] += Chara_L[m][k];  // Left state: u_L = R * w_L
                            conserr[k][i][j] += Chara_R[m][k];  // Right state: u_R = R * w_R
                        }
                    }
                }
            }
        }
        else if (dir == 2) {
            // Y-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=0.0, ny=1.0)
            Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in y-direction
            #pragma omp parallel for
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    // Local arrays for characteristic decomposition
                    double Chara_L[4][4];
                    double Chara_R[4][4];
                    double Characteristic_Variable_L;
                    double Characteristic_Variable_R;
                    
                    // Project conservative variables onto characteristic space
                    // Same mathematical operations as x-direction, but applied to y-direction interfaces
                    for (int k = 0; k < var; k++) {
                        Characteristic_Variable_L = 0.0;
                        Characteristic_Variable_R = 0.0;
                        
                        // Matrix-vector multiplication: w = L * u
                        for (int m = 0; m < var; m++) {
                            Characteristic_Variable_L += y[m][i][j] * Eigen_L[k][m][i][j];
                            Characteristic_Variable_R += y[m][i][j + 1] * Eigen_L[k][m][i][j];
                        }
                        
                        // Matrix-vector multiplication: u' = R * w
                        // Note: Eigen_R[m][k] accesses the (m,k) element of right eigenvector matrix
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Reconstruct states at interface
                    // Complete the transformation back to conservative variables
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                        
                        // Final summation: u' = Σ (characteristic contributions)
                        // This should result in the original conservative variables
                        // since R * L = I (if eigenvectors are properly normalized)
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];
                            conserr[k][i][j] += Chara_R[m][k];
                        }
                    }
                }
            }
        }
        
        // Free allocated memory for characteristic decomposition
        free(Pri); 
        free(Eigen_L); 
        free(Eigen_R);
    } 
    else {
        // Standard reconstruction without characteristic decomposition
        // Godunov scheme: piecewise constant reconstruction (1st order)
        if (dir == 1) {
            // X-direction reconstruction (piecewise constant)
            // For 1st order Godunov: u_{i+1/2}^- = u_i, u_{i+1/2}^+ = u_{i+1}
            #pragma omp parallel for collapse(3)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Left state at interface i+1/2: simply cell-centered value at i
                        conserl[k][i][j] = y[k][i][j];
                        
                        // Right state at interface i+1/2: simply cell-centered value at i+1
                        conserr[k][i][j] = y[k][i + 1][j];
                    }
                }
            }
        } 
        else if (dir == 2) {
            // Y-direction reconstruction (piecewise constant)
            // For 1st order Godunov: u_{j+1/2}^- = u_j, u_{j+1/2}^+ = u_{j+1}
            #pragma omp parallel for collapse(3)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Left state at interface j+1/2: cell-centered value at j
                        conserl[k][i][j] = y[k][i][j];
                        
                        // Right state at interface j+1/2: cell-centered value at j+1
                        conserr[k][i][j] = y[k][i][j + 1];
                    }
                }
            }
        }
    }


}

/*                                      ******************                                          */
/*                                            TVD                                                   */
/*                                      ******************                                          */
/**
 * Second-order TVD (Total Variation Diminishing) reconstruction
 * Reconstructs left and right states at cell interfaces using TVD limiters
 * Can perform reconstruction in characteristic variables for better numerical stability
 * 
 * @param dir Direction of reconstruction: 1 for x-direction, 2 for y-direction
 * @param var Number of variables (typically 4 for Euler equations)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells
 * @param y Input conservative variables array [var][rows][cols]
 * @param conserl Output: left reconstructed states at cell interfaces
 * @param conserr Output: right reconstructed states at cell interfaces
 * @param delta_x Grid spacing in x-direction
 * @param delta_y Grid spacing in y-direction
 */
static inline void TVD_Reconstruction(int dir, int var, int rows, int cols, int GC, 
                                      double (*y)[rows][cols],
                                      double (*conserl)[rows][cols], 
                                      double (*conserr)[rows][cols], 
                                      double delta_x, double delta_y) {
    
    // Check if characteristic decomposition is enabled
    if (Characteriz) {
        // Dynamic memory allocation for characteristic decomposition arrays
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));        // Primitive variables
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Left eigenvectors
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Right eigenvectors
        
        // Check memory allocation
        if (!Pri || !Eigen_L || !Eigen_R) {
            fprintf(stderr, "Memory allocation failed in TVD_Reconstruction\n");
            free(Pri); free(Eigen_L); free(Eigen_R);
            return;
        }

        // Convert conservative variables to primitive variables
        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            // X-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=1.0, ny=0.0)
            Compute_Eigen_2DX(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in x-direction
            #pragma omp parallel for collapse(2)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    // Local arrays for characteristic decomposition
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Single characteristic variable (left)
                    double Characteristic_Variable_R;  // Single characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [i-2, i-1, i, i+1, i+2, i+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil
                        // uu array indices: 0=i-2, 1=i-1, 2=i, 3=i+1, 4=i+2, 5=i+3
                        for (int nn = i - 2; nn <= i + 3; nn++) {
                            uu[nn - i + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            for (int m = 0; m < var; m++) {
                                uu[nn - i + 2] += y[m][nn][j] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply TVD reconstruction in characteristic space
                        // Uses van Leer limiter for better accuracy than minmod
                        Characteristic_Variable_L = TVD_minmod_L(&uu[2], delta_x);
                        Characteristic_Variable_R = TVD_minmod_R(&uu[2], delta_x);

                        // Transform reconstructed characteristic variables back to conservative space
                        // u' = R * w, where w is the reconstructed characteristic variable
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                    
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];  // Left state: u_L = R * w_L
                            conserr[k][i][j] += Chara_R[m][k];  // Right state: u_R = R * w_R
                        }
                    }
                }
            }
        }
        else if (dir == 2) {
            // Y-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=0.0, ny=1.0)
            //Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in y-direction
            #pragma omp parallel for collapse(2)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Characteristic variable (left)
                    double Characteristic_Variable_R;  // Characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [j-2, j-1, j, j+1, j+2, j+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil in y-direction
                        for (int nn = j - 2; nn <= j + 3; nn++) {
                            uu[nn - j + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            for (int m = 0; m < var; m++) {
                                uu[nn - j + 2] += y[m][i][nn] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply TVD reconstruction in characteristic space
                        Characteristic_Variable_L = TVD_minmod_L(&uu[2], delta_y);
                        Characteristic_Variable_R = TVD_minmod_R(&uu[2], delta_y);

                        // Transform back to conservative space: u' = R * w
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                        
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];
                            conserr[k][i][j] += Chara_R[m][k];
                        }
                    }
                }
            }
        }
        
        // Free allocated memory for characteristic decomposition
        free(Pri); 
        free(Eigen_L); 
        free(Eigen_R);
    } 
    else {
        // Standard reconstruction without characteristic decomposition
        // Direct TVD reconstruction in conservative variable space
        
        if (dir == 1) {
            // X-direction reconstruction using minmod limiter
            // 5-point stencil: [i-2, i-1, i, i+1, i+2, i+3] for boundary handling
            
            #pragma omp parallel for collapse(3)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 5-point stencil for conservative variable k
                        double fu[6];  // Stencil: fu[0]=i-2, fu[1]=i-1, fu[2]=i, fu[3]=i+1, fu[4]=i+2, fu[5]=i+3
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i - 2 + nn][j];
                        }
                        
                        // Apply TVD reconstruction using minmod limiter
                        // &fu[2] points to the central 3-point stencil [i-1, i, i+1]
                        conserl[k][i][j] = TVD_minmod_L(&fu[2], delta_x);  // Left state at interface i+1/2
                        conserr[k][i][j] = TVD_minmod_R(&fu[2], delta_x);  // Right state at interface i+1/2
                    }
                }
            }
        } 
        else if (dir == 2) {
            // Y-direction reconstruction using VanAlbda limiter
            
            #pragma omp parallel for collapse(3)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 5-point stencil for conservative variable k in y-direction
                        double fu[6];  // Stencil: fu[0]=j-2, fu[1]=j-1, fu[2]=j, fu[3]=j+1, fu[4]=j+2, fu[5]=j+3
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i][j - 2 + nn];
                        }
                        
                        // Apply TVD reconstruction using VanAlbda limiter
                        // &fu[2] points to the central 3-point stencil [j-1, j, j+1]
                        conserl[k][i][j] = TVD_minmod_L(&fu[2], delta_y);  // Left state at interface j+1/2
                        conserr[k][i][j] = TVD_minmod_R(&fu[2], delta_y);  // Right state at interface j+1/2
                    }
                }
            }
        }
    }
}


static inline void TVD_Reconstruction_OED(int dir, int var, int rows, int cols, int GC, 
                                      double (*y)[rows][cols],
                                      double (*conserl)[rows][cols], 
                                      double (*conserr)[rows][cols], 
                                      double delta_x, double delta_y) {
    
    // Check if characteristic decomposition is enabled
    if (Characteriz) {
        // Dynamic memory allocation for characteristic decomposition arrays
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));        // Primitive variables
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Left eigenvectors
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Right eigenvectors
        
        // Check memory allocation
        if (!Pri || !Eigen_L || !Eigen_R) {
            fprintf(stderr, "Memory allocation failed in TVD_Reconstruction\n");
            free(Pri); free(Eigen_L); free(Eigen_R);
            return;
        }

        // Convert conservative variables to primitive variables
        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            // X-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=1.0, ny=0.0)
            Compute_Eigen_2DX(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in x-direction
            #pragma omp parallel for collapse(2)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    // Local arrays for characteristic decomposition
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Single characteristic variable (left)
                    double Characteristic_Variable_R;  // Single characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [i-2, i-1, i, i+1, i+2, i+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil
                        // uu array indices: 0=i-2, 1=i-1, 2=i, 3=i+1, 4=i+2, 5=i+3
                        for (int nn = i - 2; nn <= i + 3; nn++) {
                            uu[nn - i + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            for (int m = 0; m < var; m++) {
                                uu[nn - i + 2] += y[m][nn][j] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply TVD reconstruction in characteristic space
                        // Uses van Leer limiter for better accuracy than minmod
                        Characteristic_Variable_L = TVD_VanAlbada_L(&uu[2], delta_x);
                        Characteristic_Variable_R = TVD_VanAlbada_R(&uu[2], delta_x);

                        // Transform reconstructed characteristic variables back to conservative space
                        // u' = R * w, where w is the reconstructed characteristic variable
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                    
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];  // Left state: u_L = R * w_L
                            conserr[k][i][j] += Chara_R[m][k];  // Right state: u_R = R * w_R
                        }
                    }
                }
            }
        }
        else if (dir == 2) {
            // Y-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=0.0, ny=1.0)
            //Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in y-direction
            #pragma omp parallel for collapse(2)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Characteristic variable (left)
                    double Characteristic_Variable_R;  // Characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [j-2, j-1, j, j+1, j+2, j+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil in y-direction
                        for (int nn = j - 2; nn <= j + 3; nn++) {
                            uu[nn - j + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            for (int m = 0; m < var; m++) {
                                uu[nn - j + 2] += y[m][i][nn] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        //奇偶失联测试代码
                        double perturbation = 2e-3;
                        Characteristic_Variable_L = OED_TVD_VanAlbada_L(&uu[2], delta_y-perturbation, delta_y+perturbation);
                        Characteristic_Variable_R = OED_TVD_VanAlbada_R(&uu[2], delta_y+perturbation, delta_y-perturbation);
                        /*if (j % 2 == 0){
                            // Apply TVD reconstruction in characteristic space
                            Characteristic_Variable_L = OED_TVD_VanAlbada_L(&uu[2], delta_y+perturbation, delta_y-perturbation);
                            Characteristic_Variable_R = OED_TVD_VanAlbada_R(&uu[2], delta_y);
                        }
                        else{
                            // Apply TVD reconstruction in characteristic space
                            Characteristic_Variable_L = TVD_VanAlbada_L(&uu[2], delta_y);
                            Characteristic_Variable_R = TVD_VanAlbada_R(&uu[2], delta_y);
                        }*/
                        
                        

                        // Transform back to conservative space: u' = R * w
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                        
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];
                            conserr[k][i][j] += Chara_R[m][k];
                        }
                    }
                }
            }
        }
        
        // Free allocated memory for characteristic decomposition
        free(Pri); 
        free(Eigen_L); 
        free(Eigen_R);
    } 
    else {
        // Standard reconstruction without characteristic decomposition
        // Direct TVD reconstruction in conservative variable space
        
        if (dir == 1) {
            // X-direction reconstruction using minmod limiter
            // 5-point stencil: [i-2, i-1, i, i+1, i+2, i+3] for boundary handling
            
            #pragma omp parallel for collapse(3)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 5-point stencil for conservative variable k
                        double fu[6];  // Stencil: fu[0]=i-2, fu[1]=i-1, fu[2]=i, fu[3]=i+1, fu[4]=i+2, fu[5]=i+3
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i - 2 + nn][j];
                        }
                        
                        // Apply TVD reconstruction using minmod limiter
                        // &fu[2] points to the central 3-point stencil [i-1, i, i+1]
                        conserl[k][i][j] = TVD_minmod_L(&fu[2], delta_x);  // Left state at interface i+1/2
                        conserr[k][i][j] = TVD_minmod_R(&fu[2], delta_x);  // Right state at interface i+1/2
                    }
                }
            }
        } 
        else if (dir == 2) {
            // Y-direction reconstruction using minmod limiter
            
            #pragma omp parallel for collapse(3)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 5-point stencil for conservative variable k in y-direction
                        double fu[6];  // Stencil: fu[0]=j-2, fu[1]=j-1, fu[2]=j, fu[3]=j+1, fu[4]=j+2, fu[5]=j+3
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i][j - 2 + nn];
                        }
                        
                        // Apply TVD reconstruction using minmod limiter
                        // &fu[2] points to the central 3-point stencil [j-1, j, j+1]
                        //  = TVD_minmod_L(&fu[2], delta_y);  // Left state at interface j+1/2
                        //  = TVD_minmod_R(&fu[2], delta_y);  // Right state at interface j+1/2
                        double perturbation = 0.0;
                        //conserl[k][i][j] = OED_TVD_VanAlbada_L(&fu[2], delta_y-perturbation, delta_y+perturbation);
                        //conserr[k][i][j] = OED_TVD_VanAlbada_R(&fu[2], delta_y+perturbation, delta_y-perturbation);

                        // 从边界来开始计算，在j = 3时，dy——j = delta_y-perturbation，此时dy——j+1 = delta_y+perturbation
                        if (j % 2 == 0){
                            // Apply TVD reconstruction in characteristic space
                            conserl[k][i][j] = OED_TVD_VanAlbada_L(&fu[2], delta_y+perturbation, delta_y-perturbation);
                            conserr[k][i][j] = OED_TVD_VanAlbada_R(&fu[2], delta_y-perturbation, delta_y+perturbation);
                            
                        }
                        else{
                            // Apply TVD reconstruction in characteristic space
                            conserl[k][i][j] = OED_TVD_VanAlbada_L(&fu[2], delta_y-perturbation, delta_y+perturbation);
                            conserr[k][i][j] = OED_TVD_VanAlbada_R(&fu[2], delta_y+perturbation, delta_y-perturbation);
                        }
                    }
                }
            }
        }
    }
}

/*                                      ******************                                          */
/*                                               WENO                                               */
/*                                      ******************                                          */

/**
 * Third-order WENO (Weighted Essentially Non-Oscillatory) reconstruction
 * Reconstructs left and right states at cell interfaces using WENO-3 scheme
 * Can perform reconstruction in characteristic variables for better shock capturing
 * 
 * @param dir Direction of reconstruction: 1 for x-direction, 2 for y-direction
 * @param var Number of variables (typically 4 for Euler equations)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells
 * @param y Input conservative variables array [var][rows][cols]
 * @param conserl Output: left reconstructed states at cell interfaces
 * @param conserr Output: right reconstructed states at cell interfaces
 */
static inline void WENO3_Reconstruction(int dir, int var, int rows, int cols, int GC, 
                                      double (*y)[rows][cols],
                                      double (*conserl)[rows][cols], 
                                      double (*conserr)[rows][cols]) {
    
    // Check if characteristic decomposition is enabled
    if (Characteriz) {
        // Dynamic memory allocation for characteristic decomposition arrays
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));        // Primitive variables
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Left eigenvectors
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Right eigenvectors
        
        // Check memory allocation
        if (!Pri || !Eigen_L || !Eigen_R) {
            fprintf(stderr, "Memory allocation failed in WENO3_Reconstruction\n");
            free(Pri); free(Eigen_L); free(Eigen_R);
            return;
        }

        // Convert conservative variables to primitive variables
        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            // X-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=1.0, ny=0.0)
            Compute_Eigen_2DX(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in x-direction
            #pragma omp parallel for
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    // Local arrays for characteristic decomposition
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Single characteristic variable (left)
                    double Characteristic_Variable_R;  // Single characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [i-2, i-1, i, i+1, i+2, i+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil
                        // uu array indices: 0=i-2, 1=i-1, 2=i, 3=i+1, 4=i+2, 5=i+3 (extended for WENO boundary handling)
                        for (int nn = i - 2; nn <= i + 3; nn++) {
                            uu[nn - i + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            // w_k = Σ_{m=0}^{3} (L_km * u_m) for each grid point nn
                            for (int m = 0; m < var; m++) {
                                uu[nn - i + 2] += y[m][nn][j] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply WENO-3 reconstruction in characteristic space
                        // Uses 3-point stencil centered at positions 1,2,3 of uu array (indices: i-1, i, i+1)
                        Characteristic_Variable_L = WENO3_L(&uu[2]);  // Left interface value w_{i+1/2}^-
                        Characteristic_Variable_R = WENO3_R(&uu[2]);  // Right interface value w_{i+1/2}^+

                        // Transform reconstructed characteristic variables back to conservative space
                        // u' = R * w, where w is the reconstructed characteristic variable
                        for (int m = 0; m < var; m++) {
                            // Chara_L[k][m] = contribution of k-th characteristic field to m-th conservative variable
                            // = R_mk * w_k, where w_k = Characteristic_Variable_L
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                    
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        // u_k = Σ_{m=0}^{3} (Chara_L[m][k]) = Σ_{m=0}^{3} (R_km * w_m_L)
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];  // Left state: u_L = R * w_L
                            conserr[k][i][j] += Chara_R[m][k];  // Right state: u_R = R * w_R
                        }
                    }
                }
            }
        }
        else if (dir == 2) {
            // Y-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=0.0, ny=1.0)
            Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in y-direction
            #pragma omp parallel for
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Characteristic variable (left)
                    double Characteristic_Variable_R;  // Characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [j-2, j-1, j, j+1, j+2, j+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil in y-direction
                        for (int nn = j - 2; nn <= j + 3; nn++) {
                            uu[nn - j + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            for (int m = 0; m < var; m++) {
                                uu[nn - j + 2] += y[m][i][nn] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply WENO-3 reconstruction in characteristic space
                        Characteristic_Variable_L = WENO3_L(&uu[2]);  // Left interface value w_{j+1/2}^-
                        Characteristic_Variable_R = WENO3_R(&uu[2]);  // Right interface value w_{j+1/2}^+

                        // Transform back to conservative space: u' = R * w
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                        
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];
                            conserr[k][i][j] += Chara_R[m][k];
                        }
                    }
                }
            }
        }
        
        // Free allocated memory for characteristic decomposition
        free(Pri); 
        free(Eigen_L); 
        free(Eigen_R);
    } 
    else {
        // Standard reconstruction without characteristic decomposition
        // Direct WENO-3 reconstruction in conservative variable space
        
        if (dir == 1) {
            // X-direction reconstruction using WENO-3 scheme
            
            #pragma omp parallel for collapse(3)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 5-point stencil for conservative variable k
                        // fu array: [i-2, i-1, i, i+1, i+2, i+3] for boundary handling
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i - 2 + nn][j];
                        }
                        
                        // Apply WENO-3 reconstruction using central 3-point stencil
                        // &fu[2] points to the stencil [i-1, i, i+1] for WENO-3 reconstruction
                        conserl[k][i][j] = WENO3_L(&fu[2]);  // Left state at interface i+1/2
                        conserr[k][i][j] = WENO3_R(&fu[2]);  // Right state at interface i+1/2
                    }
                }
            }
        } 
        else if (dir == 2) {
            // Y-direction reconstruction using WENO-3 scheme
            
            #pragma omp parallel for collapse(3)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 5-point stencil for conservative variable k in y-direction
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i][j - 2 + nn];
                        }
                        
                        // Apply WENO-3 reconstruction
                        conserl[k][i][j] = WENO3_L(&fu[2]);  // Left state at interface j+1/2
                        conserr[k][i][j] = WENO3_R(&fu[2]);  // Right state at interface j+1/2
                    }
                }
            }
        }
    }
}

/**
 * Fifth-order WENO (Weighted Essentially Non-Oscillatory) reconstruction
 * Reconstructs left and right states at cell interfaces using WENO-5 scheme
 * Provides higher-order accuracy while maintaining non-oscillatory properties near discontinuities
 * 
 * @param dir Direction of reconstruction: 1 for x-direction, 2 for y-direction
 * @param var Number of variables (typically 4 for Euler equations)
 * @param rows Number of rows in the grid (including ghost cells)
 * @param cols Number of columns in the grid (including ghost cells)
 * @param GC Number of ghost cells
 * @param y Input conservative variables array [var][rows][cols]
 * @param conserl Output: left reconstructed states at cell interfaces
 * @param conserr Output: right reconstructed states at cell interfaces
 */
static inline void WENO5_Reconstruction(int dir, int var, int rows, int cols, int GC, 
                                      double (*y)[rows][cols],
                                      double (*conserl)[rows][cols], 
                                      double (*conserr)[rows][cols]) {
    
    // Check if characteristic decomposition is enabled
    if (Characteriz) {
        // Dynamic memory allocation for characteristic decomposition arrays
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));        // Primitive variables
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Left eigenvectors
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));  // Right eigenvectors
        
        // Check memory allocation
        if (!Pri || !Eigen_L || !Eigen_R) {
            fprintf(stderr, "Memory allocation failed in WENO5_Reconstruction\n");
            free(Pri); free(Eigen_L); free(Eigen_R);
            return;
        }

        // Convert conservative variables to primitive variables
        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            // X-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=1.0, ny=0.0)
            Compute_Eigen_2DX(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in x-direction
            #pragma omp parallel for
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    // Local arrays for characteristic decomposition
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Single characteristic variable (left)
                    double Characteristic_Variable_R;  // Single characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [i-2, i-1, i, i+1, i+2, i+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil (WENO-5 requires 5 points)
                        // uu array indices: 0=i-2, 1=i-1, 2=i, 3=i+1, 4=i+2, 5=i+3 (6 points for boundary handling)
                        for (int nn = i - 2; nn <= i + 3; nn++) {
                            uu[nn - i + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            // w_k = Σ_{m=0}^{3} (L_km * u_m) for each grid point nn
                            for (int m = 0; m < var; m++) {
                                uu[nn - i + 2] += y[m][nn][j] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply WENO-5 reconstruction in characteristic space
                        // Uses 5-point stencil centered at positions 0-4 of uu array (indices: i-2, i-1, i, i+1, i+2)
                        Characteristic_Variable_L = WENO5Z_L(&uu[2]);  // Left interface value w_{i+1/2}^-
                        Characteristic_Variable_R = WENO5Z_R(&uu[2]);  // Right interface value w_{i+1/2}^+

                        // Transform reconstructed characteristic variables back to conservative space
                        // u' = R * w, where w is the reconstructed characteristic variable
                        for (int m = 0; m < var; m++) {
                            // Chara_L[k][m] = contribution of k-th characteristic field to m-th conservative variable
                            // = R_mk * w_k, where w_k = Characteristic_Variable_L
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
                    
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        // u_k = Σ_{m=0}^{3} (Chara_L[m][k]) = Σ_{m=0}^{3} (R_km * w_m_L)
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];  // Left state: u_L = R * w_L
                            conserr[k][i][j] += Chara_R[m][k];  // Right state: u_R = R * w_R
                        }
                    }
                }
            }
        }
        else if (dir == 2) {
            // Y-direction reconstruction with characteristic decomposition
            
            // Compute eigenvectors for normal direction (nx=0.0, ny=1.0)
            Compute_Eigen_2DY(var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // Process each cell interface in y-direction
            #pragma omp parallel for
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    double Chara_L[4][4];  // Left state in characteristic space
                    double Chara_R[4][4];  // Right state in characteristic space
                    double Characteristic_Variable_L;  // Characteristic variable (left)
                    double Characteristic_Variable_R;  // Characteristic variable (right)
                    double uu[10];  // Array for characteristic variables in 5-point stencil [j-2, j-1, j, j+1, j+2, j+3]
                    
                    // Process each characteristic field separately
                    for (int k = 0; k < var; k++) {
                        // Extract characteristic variables for 5-point stencil in y-direction
                        for (int nn = j - 2; nn <= j + 3; nn++) {
                            uu[nn - j + 2] = 0.0;
                            // Project conservative variables to characteristic space: w = L * u
                            for (int m = 0; m < var; m++) {
                                uu[nn - j + 2] += y[m][i][nn] * Eigen_L[k][m][i][j];
                            }
                        }
                        
                        // Apply WENO-5 reconstruction in characteristic space
                        Characteristic_Variable_L = WENO5Z_L(&uu[2]);  // Left interface value w_{j+1/2}^-
                        Characteristic_Variable_R = WENO5Z_R(&uu[2]);  // Right interface value w_{j+1/2}^+
                    
                        // Transform back to conservative space: u' = R * w
                        for (int m = 0; m < var; m++) {
                            Chara_L[k][m] = Characteristic_Variable_L * Eigen_R[m][k][i][j];
                            Chara_R[k][m] = Characteristic_Variable_R * Eigen_R[m][k][i][j];
                        }
                    }
                    
                    // Combine contributions from all characteristic fields
                    for (int k = 0; k < var; k++) {
                        conserl[k][i][j] = 0.0;
                        conserr[k][i][j] = 0.0;
        
                        // Sum contributions: u = Σ (R * w) for each characteristic field
                        for (int m = 0; m < var; m++) {
                            conserl[k][i][j] += Chara_L[m][k];
                            conserr[k][i][j] += Chara_R[m][k];
                        }
                    }
                }
            }
        }
        
        // Free allocated memory for characteristic decomposition
        free(Pri); 
        free(Eigen_L); 
        free(Eigen_R);
    } 
    else {
        // Standard reconstruction without characteristic decomposition
        // Direct WENO-5 reconstruction in conservative variable space
        
        if (dir == 1) {
            // X-direction reconstruction using WENO-5 scheme
            
            #pragma omp parallel for collapse(3)
            for (int i = GC - 1; i < rows - GC; i++) {
                for (int j = GC; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 6-point stencil for conservative variable k (5 points for WENO-5 + 1 for boundary)
                        // fu array: [i-2, i-1, i, i+1, i+2, i+3]
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i - 2 + nn][j];
                        }
                        
                        // Apply WENO-5 reconstruction using 5-point stencil
                        // &fu[2] points to the stencil [i-2, i-1, i, i+1, i+2] for WENO-5 reconstruction
                        // Note: WENO5_L and WENO5_R functions internally use 5-point stencil
                        conserl[k][i][j] = WENO5_L(&fu[2]);  // Left state at interface i+1/2
                        conserr[k][i][j] = WENO5_R(&fu[2]);  // Right state at interface i+1/2
                    }
                }
            }
        } 
        else if (dir == 2) {
            // Y-direction reconstruction using WENO-5 scheme
            
            #pragma omp parallel for collapse(3)
            for (int i = GC; i < rows - GC; i++) {
                for (int j = GC - 1; j < cols - GC; j++) {
                    for (int k = 0; k < var; k++) {
                        // Extract 6-point stencil for conservative variable k in y-direction
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[k][i][j - 2 + nn];
                        }
                        
                        // Apply WENO-5 reconstruction
                        conserl[k][i][j] = WENO5_L(&fu[2]);  // Left state at interface j+1/2
                        conserr[k][i][j] = WENO5_R(&fu[2]);  // Right state at interface j+1/2
                    }
                }
            }
        }
    }
}
/*                                      ******************                                          */
/*                                      Reconstruction Functions                                    */
/*                                      ******************                                          */

/**
 * TVD reconstruction using minmod limiter (left interface)
 * Minmod limiter: slope = 0 if signs differ, otherwise min(|slope1|, |slope2|) * sign(slope1)
 * 
 * @param f Pointer to array of cell-centered values [i-1, i, i+1]
 * @param delta Grid spacing
 * @return Reconstructed value at left cell interface (i+1/2)
 */
double TVD_minmod_L(double *f, double delta)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i-1, i, i+1]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 1);  // Value at cell i-1
    v2 = *(f + k);      // Value at cell i (center)
    v3 = *(f + k + 1);  // Value at cell i+1
    
    // Compute slopes
    slope1 = (v2 - v1) / delta;  // Left slope
    slope2 = (v3 - v2) / delta;  // Right slope
    
    // Apply minmod limiter
    slope = min_mod(slope1, slope2);
    
    // Reconstruct value at left interface: u_{i+1/2}^- = u_i + 0.5 * slope * Δx
    return v2 + 0.5 * slope * delta;
}

/**
 * TVD reconstruction using minmod limiter (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i, i+1, i+2]
 * @param delta Grid spacing
 * @return Reconstructed value at right cell interface (i+1/2)
 */

double TVD_minmod_R(double *f, double delta)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i, i+1, i+2]
    k = 1;  // Stencil centered at cell i+1
    v1 = *(f + k - 1);  // Value at cell i
    v2 = *(f + k);      // Value at cell i+1 (center)
    v3 = *(f + k + 1);  // Value at cell i+2
    
    // Compute slopes
    slope1 = (v2 - v1) / delta;  // Left slope
    slope2 = (v3 - v2) / delta;  // Right slope
    
    // Apply minmod limiter
    slope = min_mod(slope1, slope2);

    // Reconstruct value at right interface: u_{i+1/2}^+ = u_{i+1} - 0.5 * slope * Δx
    return v2 - 0.5 * slope * delta;
}

/**
 * TVD reconstruction using van Leer limiter (left interface)
 * Van Leer limiter: harmonic mean of slopes (less diffusive than minmod)
 * 
 * @param f Pointer to array of cell-centered values [i-1, i, i+1]
 * @param delta Grid spacing
 * @return Reconstructed value at left cell interface (i+1/2)
 */
double TVD_vanleer_L(double *f, double delta)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i-1, i, i+1]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 1);  // Value at cell i-1
    v2 = *(f + k);      // Value at cell i (center)
    v3 = *(f + k + 1);  // Value at cell i+1
    
    // Compute slopes
    slope1 = (v2 - v1) / delta;  // Left slope
    slope2 = (v3 - v2) / delta;  // Right slope
    
    // Apply van Leer limiter
    slope = van_leer(slope1, slope2);

    // Reconstruct value at left interface: u_{i+1/2}^- = u_i + 0.5 * slope * Δx
    return v2 + 0.5 * slope * delta;
}

/**
 * TVD reconstruction using van Leer limiter (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i, i+1, i+2]
 * @param delta Grid spacing
 * @return Reconstructed value at right cell interface (i+1/2)
 */
double TVD_vanleer_R(double *f, double delta)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i, i+1, i+2]
    k = 1;  // Stencil centered at cell i+1
    v1 = *(f + k - 1);  // Value at cell i
    v2 = *(f + k);      // Value at cell i+1 (center)
    v3 = *(f + k + 1);  // Value at cell i+2
    
    // Compute slopes
    slope1 = (v2 - v1) / delta;  // Left slope
    slope2 = (v3 - v2) / delta;  // Right slope
    
    // Apply van Leer limiter
    slope = van_leer(slope1, slope2);

    // Reconstruct value at right interface: u_{i+1/2}^+ = u_{i+1} - 0.5 * slope * Δx
    return v2 - 0.5 * slope * delta;
}


/**
 * TVD reconstruction using van Leer limiter (left interface)
 * Van Albada limiter: harmonic mean of slopes
 * 
 * @param f Pointer to array of cell-centered values [i-1, i, i+1]
 * @param delta Grid spacing
 * @return Reconstructed value at left cell interface (i+1/2)
 */
double TVD_VanAlbada_L(double *f, double delta)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i-1, i, i+1]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 1);  // Value at cell i-1
    v2 = *(f + k);      // Value at cell i (center)
    v3 = *(f + k + 1);  // Value at cell i+1
    
    // Compute slopes
    slope1 = (v2 - v1) / delta;  // Left slope
    slope2 = (v3 - v2) / delta;  // Right slope
    
    // Apply van Leer limiter
    slope = van_albada(slope1, slope2);

    // Reconstruct value at left interface: u_{i+1/2}^- = u_i + 0.5 * slope * Δx
    return v2 + 0.5 * slope * delta;
}

/**
 * TVD reconstruction using van Leer limiter (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i, i+1, i+2]
 * @param delta Grid spacing
 * @return Reconstructed value at right cell interface (i+1/2)
 */
double TVD_VanAlbada_R(double *f, double delta)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i, i+1, i+2]
    k = 1;  // Stencil centered at cell i+1
    v1 = *(f + k - 1);  // Value at cell i
    v2 = *(f + k);      // Value at cell i+1 (center)
    v3 = *(f + k + 1);  // Value at cell i+2
    
    // Compute slopes
    slope1 = (v2 - v1) / delta;  // Left slope
    slope2 = (v3 - v2) / delta;  // Right slope
    
    // Apply van Leer limiter
    slope = van_albada(slope1, slope2);

    // Reconstruct value at right interface: u_{i+1/2}^+ = u_{i+1} - 0.5 * slope * Δx
    return v2 - 0.5 * slope * delta;
}

double OED_TVD_VanAlbada_L(double *f, double delta_a, double delta_b)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i-1, i, i+1]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 1);  // Value at cell i-1
    v2 = *(f + k);      // Value at cell i (center)
    v3 = *(f + k + 1);  // Value at cell i+1
    
    // Compute slopes
    slope1 = (v2 - v1) / delta_a;  // Left slope
    slope2 = (v3 - v2) / delta_b;  // Right slope
    
    // Apply van Leer limiter
    slope = van_albada(slope1, slope2);

    // Reconstruct value at left interface: u_{i+1/2}^- = u_i + 0.5 * slope * Δx
    return v2 + 0.5 * slope * delta_b;
}

/**
 * TVD reconstruction using van Leer limiter (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i, i+1, i+2]
 * @param delta Grid spacing
 * @return Reconstructed value at right cell interface (i+1/2)
 */
double OED_TVD_VanAlbada_R(double *f, double delta_a, double delta_b)
{
    int k;
    double v1, v2, v3;
    double slope1, slope2, slope;

    // Assign values to v1, v2, v3 for stencil [i, i+1, i+2]
    k = 1;  // Stencil centered at cell i+1
    v1 = *(f + k - 1);  // Value at cell i
    v2 = *(f + k);      // Value at cell i+1 (center)
    v3 = *(f + k + 1);  // Value at cell i+2
    
    // Compute slopes
    slope1 = (v2 - v1) / delta_a;  // Left slope
    slope2 = (v3 - v2) / delta_b;  // Right slope
    
    // Apply van Leer limiter
    slope = van_albada(slope1, slope2);

    // Reconstruct value at right interface: u_{i+1/2}^+ = u_{i+1} - 0.5 * slope * Δx
    return v2 - 0.5 * slope * delta_a;
}

/**
 * WENO-3 reconstruction (left interface)
 * Third-order Weighted Essentially Non-Oscillatory scheme
 * Uses 3-point stencil: [i-1, i, i+1]
 * 
 * @param f Pointer to array of cell-centered values [i-1, i, i+1]
 * @return Reconstructed value at left cell interface (i+1/2)
 */
static inline double WENO3_L(double *f)
{
    int k;
    double v1, v2, v3;
    double s1, s2;
    double a1, a2, a3, w1, w2, w3;
    double epsilon = 1e-6;  // Small constant to avoid division by zero

    // Assign values to v1, v2, v3 for stencil [i-1, i, i+1]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 1);  // Value at cell i-1
    v2 = *(f + k);      // Value at cell i (center)
    v3 = *(f + k + 1);  // Value at cell i+1

    // Compute smoothness indicators (measure of solution variation)
    s1 = (v2 - v1) * (v2 - v1);  // Variation between i-1 and i
    s2 = (v3 - v2) * (v3 - v2);  // Variation between i and i+1

    // Compute nonlinear weights
    a1 = (1.0/3.0) / pow(epsilon + s1, 2);
    a2 = (2.0/3.0) / pow(epsilon + s2, 2);

    w1 = a1 / (a1 + a2);
    w2 = a2 / (a1 + a2);

    // Check for negative weights (should not happen)
    if (w1 < 0.0 || w2 < 0.0)
        printf("Negative weights appear in WENO!!!\n");
    
    // Return weighted average of second-order reconstructions
    // First polynomial: -0.5*v1 + 1.5*v2 (left-biased)
    // Second polynomial: 0.5*v2 + 0.5*v3 (centered)
    return w1 * (-0.5 * v1 + 1.5 * v2) + w2 * (0.5 * v2 + 0.5 * v3);
}

/**
 * WENO-3 reconstruction (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i+2, i+1, i]
 * @return Reconstructed value at right cell interface (i+1/2)
 */
static inline double WENO3_R(double *f)
{
    int k;
    double v1, v2, v3;
    double s1, s2;
    double a1, a2, a3, w1, w2, w3;
    double epsilon = 1e-6;

    // Assign values to v1, v2, v3 for stencil [i+2, i+1, i]
    k = 1;  // Stencil centered at cell i+1 (mirrored for right interface)
    v1 = *(f + k + 1);  // Value at cell i+2
    v2 = *(f + k);      // Value at cell i+1 (center)
    v3 = *(f + k - 1);  // Value at cell i

    // Compute smoothness indicators
    s1 = (v2 - v1) * (v2 - v1);  // Variation between i+2 and i+1
    s2 = (v3 - v2) * (v3 - v2);  // Variation between i+1 and i

    // Compute nonlinear weights
    a1 = (1.0/3.0) / pow(epsilon + s1, 2);
    a2 = (2.0/3.0) / pow(epsilon + s2, 2);

    w1 = a1 / (a1 + a2);
    w2 = a2 / (a1 + a2);

    // Check for negative weights
    if (w1 < 0.0 || w2 < 0.0)
        printf("Negative weights appear in WENO!!!\n");
    
    // Return weighted average of second-order reconstructions
    return w1 * (-0.5 * v1 + 1.5 * v2) + w2 * (0.5 * v2 + 0.5 * v3);
}

/**
 * WENO-5 reconstruction (left interface)
 * Fifth-order Weighted Essentially Non-Oscillatory scheme
 * Uses 5-point stencil: [i-2, i-1, i, i+1, i+2]
 * 
 * @param f Pointer to array of cell-centered values [i-2, i-1, i, i+1, i+2]
 * @return Reconstructed value at left cell interface (i+1/2)
 */
static inline double WENO5_L(double *f)
{
    int k;
    double v1, v2, v3, v4, v5;
    double s1, s2, s3;
    double a1, a2, a3, w1, w2, w3;
    double epsilon = 1.0e-6;

    // Assign values to v1, v2, v3, v4, v5 for stencil [i-2, i-1, i, i+1, i+2]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 2);  // Value at cell i-2
    v2 = *(f + k - 1);  // Value at cell i-1
    v3 = *(f + k);      // Value at cell i (center)
    v4 = *(f + k + 1);  // Value at cell i+1
    v5 = *(f + k + 2);  // Value at cell i+2

    // Compute smoothness indicators (Jiang & Shu, 1996)
    s1 = 13.0/12.0 * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3) 
       + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) * (v1 - 4.0 * v2 + 3.0 * v3);
    
    s2 = 13.0/12.0 * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4) 
       + 0.25 * (v2 - v4) * (v2 - v4);
    
    s3 = 13.0/12.0 * (v3 - 2.0 * v4 + v5) * (v3 - 2.0 * v4 + v5) 
       + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) * (3.0 * v3 - 4.0 * v4 + v5);

    // Compute nonlinear weights
    a1 = 0.1/(epsilon + s1)/(epsilon + s1);
	a2 = 0.6/(epsilon + s2)/(epsilon + s2);
	a3 = 0.3/(epsilon + s3)/(epsilon + s3);

    // Normalize weights
    w1 = a1 / (a1 + a2 + a3);
    w2 = a2 / (a1 + a2 + a3);
    w3 = a3 / (a1 + a2 + a3);
    
    // Check for negative weights
    if (w1 < 0.0 || w2 < 0.0 || w3 < 0.0)
        printf("Negative weights appear in WENO!!!\n");
    
    // Return weighted average of third-order reconstructions
    // First polynomial: (2v1 - 7v2 + 11v3)/6 (left-most)
    // Second polynomial: (-v2 + 5v3 + 2v4)/6 (centered)
    // Third polynomial: (2v3 + 5v4 - v5)/6 (right-most)
    return w1 * (2.0 * v1 - 7.0 * v2 + 11.0 * v3) / 6.0
         + w2 * (-v2 + 5.0 * v3 + 2.0 * v4) / 6.0
         + w3 * (2.0 * v3 + 5.0 * v4 - v5) / 6.0;
}

/**
 * WENO-5 reconstruction (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i+3, i+2, i+1, i, i-1]
 * @return Reconstructed value at right cell interface (i+1/2)
 */

static inline double WENO5_R(double *f)
{
    int k;
    double v1, v2, v3, v4, v5;
    double s1, s2, s3;
    double a1, a2, a3, w1, w2, w3;
    double epsilon = 1.0e-6;

    // Assign values to v1, v2, v3, v4, v5 for stencil [i+3, i+2, i+1, i, i-1]
    k = 1;  // Stencil centered at cell i+1 (mirrored for right interface)
    v1 = *(f + k + 2);  // Value at cell i+3
    v2 = *(f + k + 1);  // Value at cell i+2
    v3 = *(f + k);      // Value at cell i+1 (center)
    v4 = *(f + k - 1);  // Value at cell i
    v5 = *(f + k - 2);  // Value at cell i-1

    // Compute smoothness indicators (mirrored stencil)
    s1 = 13.0/12.0 * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3) 
       + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) * (v1 - 4.0 * v2 + 3.0 * v3);
    
    s2 = 13.0/12.0 * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4) 
       + 0.25 * (v2 - v4) * (v2 - v4);
    
    s3 = 13.0/12.0 * (v3 - 2.0 * v4 + v5) * (v3 - 2.0 * v4 + v5) 
       + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) * (3.0 * v3 - 4.0 * v4 + v5);

    // Compute nonlinear weights
    a1 = 0.1/(epsilon + s1)/(epsilon + s1);
	a2 = 0.6/(epsilon + s2)/(epsilon + s2);
	a3 = 0.3/(epsilon + s3)/(epsilon + s3);

    // Normalize weights
    w1 = a1 / (a1 + a2 + a3);
    w2 = a2 / (a1 + a2 + a3);
    w3 = a3 / (a1 + a2 + a3);
    
    // Check for negative weights
    if (w1 < 0.0 || w2 < 0.0 || w3 < 0.0)
        printf("Negative weights appear in WENO!!!\n");
    
    // Return weighted average of third-order reconstructions
    return w1 * (2.0 * v1 - 7.0 * v2 + 11.0 * v3) / 6.0
         + w2 * (-v2 + 5.0 * v3 + 2.0 * v4) / 6.0
         + w3 * (2.0 * v3 + 5.0 * v4 - v5) / 6.0;
}


/**
 * WENO-5 reconstruction (left interface)
 * Fifth-order Weighted Essentially Non-Oscillatory scheme
 * Uses 5-point stencil: [i-2, i-1, i, i+1, i+2]
 * 
 * @param f Pointer to array of cell-centered values [i-2, i-1, i, i+1, i+2]
 * @return Reconstructed value at left cell interface (i+1/2)
 */
static inline double WENO5Z_L(double *f)
{
    int k;
    double v1, v2, v3, v4, v5;
    double s1, s2, s3;
    double a1, a2, a3, w1, w2, w3;
    double epsilon = 1.0e-15;

    // Assign values to v1, v2, v3, v4, v5 for stencil [i-2, i-1, i, i+1, i+2]
    k = 0;  // Stencil centered at cell i
    v1 = *(f + k - 2);  // Value at cell i-2
    v2 = *(f + k - 1);  // Value at cell i-1
    v3 = *(f + k);      // Value at cell i (center)
    v4 = *(f + k + 1);  // Value at cell i+1
    v5 = *(f + k + 2);  // Value at cell i+2

    // Compute smoothness indicators (Jiang & Shu, 1996)
    s1 = 13.0/12.0 * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3) 
       + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) * (v1 - 4.0 * v2 + 3.0 * v3);
    
    s2 = 13.0/12.0 * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4) 
       + 0.25 * (v2 - v4) * (v2 - v4);
    
    s3 = 13.0/12.0 * (v3 - 2.0 * v4 + v5) * (v3 - 2.0 * v4 + v5) 
       + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) * (3.0 * v3 - 4.0 * v4 + v5);


    double tau_5 = fabs(s1-s3);

    // Compute nonlinear weights
    a1 = 0.1*(1.0  + tau_5 /(s1 +epsilon));
	a2 = 0.6*(1.0  + tau_5 /(s2 +epsilon));
	a3 = 0.3*(1.0  + tau_5 /(s3 +epsilon));

    // Normalize weights
    w1 = a1 / (a1 + a2 + a3);
    w2 = a2 / (a1 + a2 + a3);
    w3 = a3 / (a1 + a2 + a3);
    
    // Check for negative weights
    if (w1 < 0.0 || w2 < 0.0 || w3 < 0.0)
        printf("Negative weights appear in WENO!!!\n");
    
    // Return weighted average of third-order reconstructions
    // First polynomial: (2v1 - 7v2 + 11v3)/6 (left-most)
    // Second polynomial: (-v2 + 5v3 + 2v4)/6 (centered)
    // Third polynomial: (2v3 + 5v4 - v5)/6 (right-most)
    return w1 * (2.0 * v1 - 7.0 * v2 + 11.0 * v3) / 6.0
         + w2 * (-v2 + 5.0 * v3 + 2.0 * v4) / 6.0
         + w3 * (2.0 * v3 + 5.0 * v4 - v5) / 6.0;
}

/**
 * WENO-5 reconstruction (right interface)
 * 
 * @param f Pointer to array of cell-centered values [i+3, i+2, i+1, i, i-1]
 * @return Reconstructed value at right cell interface (i+1/2)
 */

static inline double WENO5Z_R(double *f)
{
    int k;
    double v1, v2, v3, v4, v5;
    double s1, s2, s3;
    double a1, a2, a3, w1, w2, w3;
    double epsilon = 1.0e-15;

    // Assign values to v1, v2, v3, v4, v5 for stencil [i+3, i+2, i+1, i, i-1]
    k = 1;  // Stencil centered at cell i+1 (mirrored for right interface)
    v1 = *(f + k + 2);  // Value at cell i+3
    v2 = *(f + k + 1);  // Value at cell i+2
    v3 = *(f + k);      // Value at cell i+1 (center)
    v4 = *(f + k - 1);  // Value at cell i
    v5 = *(f + k - 2);  // Value at cell i-1

    // Compute smoothness indicators (mirrored stencil)
    s1 = 13.0/12.0 * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3) 
       + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) * (v1 - 4.0 * v2 + 3.0 * v3);
    
    s2 = 13.0/12.0 * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4) 
       + 0.25 * (v2 - v4) * (v2 - v4);
    
    s3 = 13.0/12.0 * (v3 - 2.0 * v4 + v5) * (v3 - 2.0 * v4 + v5) 
       + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) * (3.0 * v3 - 4.0 * v4 + v5);

    // Compute nonlinear weights
    double tau_5 = fabs(s1-s3);
    // Compute nonlinear weights
    a1 = 0.1*(1.0  + tau_5 /(s1 +epsilon));
	a2 = 0.6*(1.0  + tau_5 /(s2 +epsilon));
	a3 = 0.3*(1.0  + tau_5 /(s3 +epsilon));

    // Normalize weights
    w1 = a1 / (a1 + a2 + a3);
    w2 = a2 / (a1 + a2 + a3);
    w3 = a3 / (a1 + a2 + a3);
    
    // Check for negative weights
    if (w1 < 0.0 || w2 < 0.0 || w3 < 0.0)
        printf("Negative weights appear in WENO!!!\n");
    
    // Return weighted average of third-order reconstructions
    return w1 * (2.0 * v1 - 7.0 * v2 + 11.0 * v3) / 6.0
         + w2 * (-v2 + 5.0 * v3 + 2.0 * v4) / 6.0
         + w3 * (2.0 * v3 + 5.0 * v4 - v5) / 6.0;
}




#endif  