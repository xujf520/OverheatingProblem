#ifndef FLUX_H
#define FLUX_H

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "Characteriz.h"
#include "Golbal.h"
#include "function.h"

// 定义全局变量来累加 p_star - 2.926650
static double p_star_accumulator = 0.0;
static double  Flux_test = 0.0;

//提前声明部分函数
double ExactRieamnna_pstar(double rhol, double rhor, double ul,double ur, double pl, \
							double pr, double gammal,double gammar ,double tol, double maxit);
double ExactRieamnna_ustar(double p_star, double rhol, double rhor, double ul,double ur, double pl,double pr, double gammal,double gammar);

//将中间变量第二次还原
static inline void restore2(double*u1,double*u2,double*u3,int n,double*roarr,double*uarr,double*parr){
    int i;
    for(i = 0;i < n;i++){
        roarr[i] = u1[i];
    }
    for(i = 0;i < n;i++){
        uarr[i] = u2[i]/u1[i];
    }
    for(i = 0;i < n;i++){
        parr[i] = (u3[i]-(0.5*pow(u2[i],2))/u1[i])*0.4;
    }
}






/*                               ************************************                               */
/*                               ************************************                               */
/*                              Riemann Solver：Exact and Approximate                               */
/*                               ************************************                               */
/*                               ************************************                               */

/*                                      ******************                                          */
/*                                   在固定参考系下的黎曼解法器                                        */
/*                                      ******************                                          */
/**
 * LLF (Local Lax-Friedrichs) Flux Calculator for Euler Equations (2D)
 * Computes numerical fluxes using Local Lax-Friedrichs (Rusanov) scheme
 * 
 * @param dir Direction: 1 for x-direction fluxes, 2 for y-direction fluxes
 * @param var Number of variables (currently fixed at 4)
 * @param rows Number of grid points in y-direction
 * @param cols Number of grid points in x-direction
 * @param GC Number of ghost cells
 * @param L_state Left state conservative variables at cell interfaces
 * @param R_state Right state conservative variables at cell interfaces
 * @param Flux Output flux arrays
 */
static inline void LLF_Flux(int dir, int var, int rows, int cols, int GC, 
                            double (*L_state)[rows][cols], 
                            double (*R_state)[rows][cols], 
                            double (*Flux)[rows][cols]) {
    
    // Precompute constant for efficiency (gamma-1 is frequently used)
    const double gamma_minus_one = M_gamma - 1.0;
    
    if (dir == 1) {
        // X-direction flux calculation - fully parallelized
        // LLF flux formula: F_LLF = 0.5*(F_L + F_R) - 0.5*α*(U_R - U_L)
        // where α = max(|u_L|+a_L, |u_R|+a_R) is the maximum wave speed
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j < rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                // Read left and right conservative variables (density, momentum, energy)
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables using division optimization
                // Instead of dividing multiple times, compute reciprocal once
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;      // x-velocity left
                const double u_R = rhou_R * inv_rho_R;      // x-velocity right
                const double v_L = rhov_L * inv_rho_L;      // y-velocity left
                const double v_R = rhov_R * inv_rho_R;      // y-velocity right
                
                // Compute squares of velocities for kinetic energy calculation
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law: p = (ρe - 0.5ρ(u²+v²)) * (γ-1)
                // This converts total energy to pressure
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_one;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_one;
                
                // Compute speed of sound: a = sqrt(γ * p / ρ)
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);

                // Compute total enthalpy: H = 0.5(u²+v²) + γ/(γ-1) * p/ρ
                // Enthalpy is needed for energy flux calculation
                const double gamma_ratio = M_gamma / gamma_minus_one;
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in x-direction:
                // F(ρ) = ρu
                // F(ρu) = ρu² + p
                // F(ρv) = ρuv
                // F(ρe) = ρHu
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;
            
                // Calculate maximum wave speed for numerical dissipation
                // α = max(|u_L|+a_L, |u_R|+a_R) ensures numerical stability
                const double splus = max_of_two(fabs(u_L) + a_L, fabs(u_R) + a_R);

                // Compute LLF flux using central flux plus dissipation term:
                // F_LLF = 0.5*(F_L + F_R) - 0.5*α*(U_R - U_L)
                Flux[0][j][k] = 0.5 * (rho_FL + rho_FR) - 0.5 * splus * (rho_R - rho_L); 
                Flux[1][j][k] = 0.5 * (rhou_FL + rhou_FR) - 0.5 * splus * (rhou_R - rhou_L); 
                Flux[2][j][k] = 0.5 * (rhov_FL + rhov_FR) - 0.5 * splus * (rhov_R - rhov_L);
                Flux[3][j][k] = 0.5 * (rhoe_FL + rhoe_FR) - 0.5 * splus * (rhoe_R - rhoe_L);
            }
        }
    }
    else if (dir == 2) {
        // Y-direction flux calculation - fully parallelized
        // Same LLF formula but using y-direction velocities and fluxes
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                // Velocity squares for kinetic energy
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure from energy equation
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_one;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_one;

                // Compute speed of sound in y-direction
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);

                // Compute total enthalpy for energy flux
                const double gamma_ratio = M_gamma / gamma_minus_one;
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in y-direction:
                // G(ρ) = ρv
                // G(ρu) = ρuv
                // G(ρv) = ρv² + p
                // G(ρe) = ρHv
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;
            
                // Maximum wave speed for y-direction: α = max(|v_L|+a_L, |v_R|+a_R)
                const double splus = max_of_two(fabs(v_L) + a_L, fabs(v_R) + a_R);

                // Compute LLF flux in y-direction
                Flux[0][j][k] = 0.5 * (rho_GL + rho_GR) - 0.5 * splus * (rho_R - rho_L); 
                Flux[1][j][k] = 0.5 * (rhou_GL + rhou_GR) - 0.5 * splus * (rhou_R - rhou_L); 
                Flux[2][j][k] = 0.5 * (rhov_GL + rhov_GR) - 0.5 * splus * (rhov_R - rhov_L);
                Flux[3][j][k] = 0.5 * (rhoe_GL + rhoe_GR) - 0.5 * splus * (rhoe_R - rhoe_L);
            }
        }
    }
}



/**
 * HLL Flux Calculator for Euler Equations (2D)
 * 
 * Computes numerical fluxes using HLL (Harten-Lax-van Leer) approximate Riemann solver
 * Supports both x-direction (dir=1) and y-direction (dir=2) flux calculations
 * 
 * @param dir Direction of flux calculation (1 for x-direction, 2 for y-direction)
 * @param var Number of variables (not used in current implementation but kept for compatibility)
 * @param rows Number of rows in the computational domain
 * @param cols Number of columns in the computational domain
 * @param GC Number of ghost cells
 * @param L_state Left state variables (conservative form)
 * @param R_state Right state variables (conservative form)
 * @param Flux Output flux arrays
 */
static inline void HLL_Flux(int dir, int var, int rows, int cols, int GC, 
                            double (*L_state)[rows][cols], 
                            double (*R_state)[rows][cols], 
                            double (*Flux)[rows][cols]) {

    // Check for NULL pointers
    if (L_state == NULL || R_state == NULL || Flux == NULL) {
        fprintf(stderr, "Error: NULL pointer passed to HLL_Flux function.\n");
        return;
    }
    
    if (dir == 1) {
        // x-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int i = GC-1; i < rows-GC; i++) {
            for (int j = GC; j < cols-GC; j++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][i][j];
                const double rho_R = R_state[0][i][j];
                const double rhou_L = L_state[1][i][j];
                const double rhou_R = R_state[1][i][j];
                const double rhov_L = L_state[2][i][j];
                const double rhov_R = R_state[2][i][j];
                const double rhoe_L = L_state[3][i][j];
                const double rhoe_R = R_state[3][i][j];
                
                // Error checking for density
                if (rho_L <= 0.0 || rho_R <= 0.0) {
                    fprintf(stderr, "Error: Non-positive density detected at cell (%d, %d): rho_L=%e, rho_R=%e\n", 
                            i, j, rho_L, rho_R);
                    continue; // Skip this cell or handle error appropriately
                }
                
                if (isnan(rho_L) || isnan(rho_R)) {
                    fprintf(stderr, "Error: NaN density detected at cell (%d, %d)\n", i, j);
                    continue;
                }

                // Compute primitive variables (using division optimization)
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                // Compute squares of velocities
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law
                double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * (M_gamma - 1.0);
                double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * (M_gamma - 1.0);

                // Check for negative pressure
                if (p_L <= 0.0 || p_R <= 0.0) {
                    fprintf(stderr, "Warning: Non-positive pressure at cell (%d, %d): p_L=%e, p_R=%e\n", 
                            i, j, p_L, p_R);
                    // Continue with calculation but note this may affect stability
                }

                // Compute total enthalpy H
                const double gamma_ratio = M_gamma / (M_gamma - 1.0);
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                //Calculate The Left And Right Sound Speed Values
                const double a_L = sqrt((M_gamma - 1.0) * (H_L - 0.5 * (u_L_sq + v_L_sq)));
                const double a_R = sqrt((M_gamma - 1.0) * (H_R - 0.5 * (u_R_sq + v_R_sq)));

                // Compute left and right flux vectors in x-direction
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;
            
                // Compute Roe averages for intermediate state
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;

                // Compute approximate wave speeds using Roe-averaged variables, following Einfeldt's approach to Roe eigenvalue estimation.
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;

                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R /(sum_sqrt*sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L *a_L*a_L + sqrt_rho_R * a_R*a_R)*inv_sum_sqrt + eta_sq*(u_R-u_L)*(u_R-u_L));
                
                // Check for valid wave speeds
                if (isnan(sqrt_d_sq)) {
                    fprintf(stderr, "Error: Invalid Roe-averaged speed of sound at cell (%d, %d)\n", i, j);
                    continue;
                }
                
                const double sleft = ubar - sqrt_d_sq;
                const double sright = ubar + sqrt_d_sq;

                // Determine HLL numerical flux based on wave speeds
                double rho_F, rhou_F, rhov_F, rhoe_F;
                
                if (sleft >= 0.0) {
                    // Entirely left-going waves
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if (sright <= 0.0) {
                    // Entirely right-going waves
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }
                else { 
                    // Intermediate region (sleft < 0 && sright > 0)
                    const double inv_sdiff = 1.0 / (sright - sleft);
                    const double sleft_sright = sleft * sright;
                    
                    rho_F = (sright * rho_FL - sleft * rho_FR + sleft_sright * (rho_R - rho_L)) * inv_sdiff;
                    rhou_F = (sright * rhou_FL - sleft * rhou_FR + sleft_sright * (rhou_R - rhou_L)) * inv_sdiff;
                    rhov_F = (sright * rhov_FL - sleft * rhov_FR + sleft_sright * (rhov_R - rhov_L)) * inv_sdiff;
                    rhoe_F = (sright * rhoe_FL - sleft * rhoe_FR + sleft_sright * (rhoe_R - rhoe_L)) * inv_sdiff;
                }
                
                // Store computed fluxes
                Flux[0][i][j] = rho_F; 
                Flux[1][i][j] = rhou_F; 
                Flux[2][i][j] = rhov_F;
                Flux[3][i][j] = rhoe_F;
            }
        }
    }
    else if (dir == 2) {
        // y-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int i = GC; i < rows-GC; i++) {
            for (int j = GC-1; j < cols-GC; j++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][i][j];
                const double rho_R = R_state[0][i][j];
                const double rhou_L = L_state[1][i][j];
                const double rhou_R = R_state[1][i][j];
                const double rhov_L = L_state[2][i][j];
                const double rhov_R = R_state[2][i][j];
                const double rhoe_L = L_state[3][i][j];
                const double rhoe_R = R_state[3][i][j];
                
                // Error checking for density
                if (rho_L <= 0.0 || rho_R <= 0.0) {
                    fprintf(stderr, "Error: Non-positive density detected at cell (%d, %d): rho_L=%e, rho_R=%e\n", 
                            i, j, rho_L, rho_R);
                    continue;
                }
                
                if (isnan(rho_L) || isnan(rho_R)) {
                    fprintf(stderr, "Error: NaN density detected at cell (%d, %d)\n", i, j);
                    continue;
                }

                // Compute primitive variables (using division optimization)
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                // Compute squares of velocities
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law
                double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * (M_gamma - 1.0);
                double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * (M_gamma - 1.0);
                
                // Check for negative pressure
                if (p_L <= 0.0 || p_R <= 0.0) {
                    fprintf(stderr, "Warning: Non-positive pressure at cell (%d, %d): p_L=%e, p_R=%e\n", 
                            i, j, p_L, p_R);
                }

                // Compute total enthalpy H
                const double gamma_ratio = M_gamma / (M_gamma - 1.0);
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);
                const double a_L = sqrt((M_gamma - 1.0) * (H_L - 0.5 * (u_L_sq + v_L_sq)));
                const double a_R = sqrt((M_gamma - 1.0) * (H_R - 0.5 * (u_R_sq + v_R_sq)));

                // Compute left and right flux vectors in y-direction
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;
            
                // Compute Roe averages for intermediate state
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;

                // Compute approximate wave speeds using Roe-averaged variables, following Einfeldt's approach to Roe eigenvalue estimation.
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;

                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R /(sum_sqrt*sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L *a_L*a_L + sqrt_rho_R * a_R*a_R)*inv_sum_sqrt + eta_sq*(v_R-v_L)*(v_R-v_L));
                
                // Check for valid wave speeds
                if (isnan(sqrt_d_sq)) {
                    fprintf(stderr, "Error: Invalid Roe-averaged speed of sound at cell (%d, %d)\n", i, j);
                    continue;
                }
                
                const double sleft = vbar - sqrt_d_sq;
                const double sright = vbar + sqrt_d_sq;

                // Determine HLL numerical flux based on wave speeds
                double rho_G, rhou_G, rhov_G, rhoe_G;
                
                if (sleft >= 0.0) {
                    // Entirely left-going waves
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if (sright <= 0.0) {
                    // Entirely right-going waves
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }
                else { 
                    // Intermediate region (sleft < 0 && sright > 0)
                    const double inv_sdiff = 1.0 / (sright - sleft);
                    const double sleft_sright = sleft * sright;
                    
                    rho_G = (sright * rho_GL - sleft * rho_GR + sleft_sright * (rho_R - rho_L)) * inv_sdiff;
                    rhou_G = (sright * rhou_GL - sleft * rhou_GR + sleft_sright * (rhou_R - rhou_L)) * inv_sdiff;
                    rhov_G = (sright * rhov_GL - sleft * rhov_GR + sleft_sright * (rhov_R - rhov_L)) * inv_sdiff;
                    rhoe_G = (sright * rhoe_GL - sleft * rhoe_GR + sleft_sright * (rhoe_R - rhoe_L)) * inv_sdiff;
                }
                
                // Store computed fluxes
                Flux[0][i][j] = rho_G; 
                Flux[1][i][j] = rhou_G; 
                Flux[2][i][j] = rhov_G;
                Flux[3][i][j] = rhoe_G;
            }
        }
    }
}

/**
 * HLLC (Harten-Lax-van Leer Contact) Flux Calculator for Euler Equations (2D)
 * Computes numerical fluxes using HLLC approximate Riemann solver with Einfeldt's wave speed estimate
 * HLLC solver captures contact discontinuities better than HLL by including intermediate wave
 * 
 * @param dir Direction: 1 for x-direction fluxes, 2 for y-direction fluxes
 * @param var Number of variables (currently fixed at 4)
 * @param rows Number of grid points in y-direction
 * @param cols Number of grid points in x-direction
 * @param GC Number of ghost cells
 * @param L_state Left state conservative variables at cell interfaces
 * @param R_state Right state conservative variables at cell interfaces
 * @param Flux Output flux arrays
 */
static inline void HLLC_Flux(int dir, int var, int rows, int cols, int GC, 
                             double (*L_state)[rows][cols], 
                             double (*R_state)[rows][cols], 
                             double (*Flux)[rows][cols]) {
    
    // Precompute frequently used constants
    const double gamma_minus_1 = M_gamma - 1.0;
    const double gamma_ratio = M_gamma / gamma_minus_1;  // γ/(γ-1) for enthalpy calculation

    if (dir == 1) {
        // X-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j < rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables using division optimization
                // Compute reciprocal once to avoid multiple divisions
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;      // x-velocity left state
                const double u_R = rhou_R * inv_rho_R;      // x-velocity right state
                const double v_L = rhov_L * inv_rho_L;      // y-velocity left state
                const double v_R = rhov_R * inv_rho_R;      // y-velocity right state
                
                // Compute velocity squares for kinetic energy calculation
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law: p = (ρe - 0.5ρ(u²+v²)) * (γ-1)
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;
                
                // Compute speed of sound: a = sqrt(γ * p / ρ)
                //const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                //const double a_R = sqrt(M_gamma * p_R * inv_rho_R);

                // Compute total enthalpy: H = 0.5(u²+v²) + γ/(γ-1) * p/ρ
                // Enthalpy is conserved across shocks and used in flux calculation
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);


                //Calculate The Left And Right Sound Speed Values
                const double a_L = sqrt((M_gamma - 1.0) * (H_L - 0.5 * (u_L_sq + v_L_sq)));
                const double a_R = sqrt((M_gamma - 1.0) * (H_R - 0.5 * (u_R_sq + v_R_sq)));

                // Compute Euler fluxes in x-direction:
                // F(ρ) = ρu, F(ρu) = ρu² + p, F(ρv) = ρuv, F(ρe) = ρHu
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;
            
                // Compute Roe-averaged variables for intermediate state estimation
                // Roe averages provide a consistent linearization of the Euler equations
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;

                // Compute approximate wave speeds using Roe-averaged variables, 
                // following Einfeldt's approach to Roe eigenvalue estimation
                // This provides more robust wave speed estimates, especially for strong shocks
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;
                
                // Compute Einfeldt's correction term for better wave speed estimation
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R / (sum_sqrt * sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (u_R - u_L) * (u_R - u_L));
                
                // Check for valid wave speeds (avoid numerical issues)
                if (isnan(sqrt_d_sq)) {
                    // Fall back to simple Roe average if Einfeldt's method fails
                    const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * (ubar_sq + vbar_sq)));
                    sqrt_d_sq = cbar;
                }
                
                // Wave speeds: s_left = ubar - â, s_right = ubar + â
                // where â is the modified sound speed from Einfeldt's method
                const double sleft = ubar - sqrt_d_sq;
                const double sright = ubar + sqrt_d_sq;
            
                // Compute intermediate wave speed s_star and pressure p_star
                // s_star represents the speed of the contact discontinuity
                double rho_F, rhou_F, rhov_F, rhoe_F;
                
                // Handle potential division by zero in s_star calculation
                const double denom = rho_L * (sleft - u_L) - rho_R * (sright - u_R);
                
                if (fabs(denom) < 1e-12) {
                    // Special case: fall back to HLL flux when denominator is near zero
                    // This provides numerical stability in degenerate cases
                    const double inv_sdiff = 1.0 / (sright - sleft);
                    const double sleft_sright = sleft * sright;
                    
                    rho_F = (sright * rho_FL - sleft * rho_FR + sleft_sright * (rho_R - rho_L)) * inv_sdiff;
                    rhou_F = (sright * rhou_FL - sleft * rhou_FR + sleft_sright * (rhou_R - rhou_L)) * inv_sdiff;
                    rhov_F = (sright * rhov_FL - sleft * rhov_FR + sleft_sright * (rhov_R - rhov_L)) * inv_sdiff;
                    rhoe_F = (sright * rhoe_FL - sleft * rhoe_FR + sleft_sright * (rhoe_R - rhoe_L)) * inv_sdiff;
                }
                else {
                    // Compute intermediate wave speed (contact discontinuity speed)
                    // Toro's formula: s_star = [p_R - p_L + ρ_L u_L(sleft - u_L) - ρ_R u_R(sright - u_R)] / denom
                    const double s_star = (p_R - p_L + rho_L * u_L * (sleft - u_L) 
                                         - rho_R * u_R * (sright - u_R)) / denom;
                    
                    // Compute intermediate pressure (constant across contact discontinuity)
                    const double p_star = p_L + rho_L * (sleft - u_L) * (s_star - u_L);
                    
                    // Compute intermediate state densities (star region)
                    const double u_stat_L = rho_L * (sleft - u_L) / (sleft - s_star);
                    const double u_stat_R = rho_R * (sright - u_R) / (sright - s_star);

                    // Compute intermediate state energy fluxes
                    const double fe_star_L = rhoe_L * inv_rho_L + (s_star - u_L) 
                                           * (s_star + p_L / (rho_L * (sleft - u_L)));
                    const double fe_star_R = rhoe_R * inv_rho_R + (s_star - u_R) 
                                           * (s_star + p_R / (rho_R * (sright - u_R)));
                    
                    // Determine HLLC numerical flux based on wave speeds
                    // Four possible cases depending on wave positions relative to zero
                    if (sleft >= 0.0) {
                        // Entirely left-going waves: use left flux
                        rho_F = rho_FL;
                        rhou_F = rhou_FL;
                        rhov_F = rhov_FL;
                        rhoe_F = rhoe_FL;
                    }
                    else if (sright <= 0.0) {
                        // Entirely right-going waves: use right flux
                        rho_F = rho_FR;
                        rhou_F = rhou_FR;
                        rhov_F = rhov_FR;
                        rhoe_F = rhoe_FR;
                    }
                    else if (s_star >= 0.0) {
                        // Left and star region: sleft < 0 ≤ s_star
                        rho_F = rho_FL + sleft * (u_stat_L - rho_L);
                        rhou_F = rhou_FL + sleft * (u_stat_L * s_star - rhou_L);
                        rhov_F = rhov_FL + sleft * (u_stat_L * v_L - rhov_L);
                        rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L - rhoe_L);
                    }
                    else {
                        // Star and right region: s_star < 0 < sright
                        rho_F = rho_FR + sright * (u_stat_R - rho_R);
                        rhou_F = rhou_FR + sright * (u_stat_R * s_star - rhou_R);
                        rhov_F = rhov_FR + sright * (u_stat_R * v_R - rhov_R);
                        rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R - rhoe_R);
                    }
                }

                // Store computed fluxes in output array
                Flux[0][j][k] = rho_F; 
                Flux[1][j][k] = rhou_F; 
                Flux[2][j][k] = rhov_F;
                Flux[3][j][k] = rhoe_F;
            }
        }
    }
    else if (dir == 2) {
        // Y-direction flux calculation - fully parallelized
        // Same algorithm as x-direction but with velocities and fluxes rotated
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                // Velocity squares for kinetic energy
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;

                // Compute speed of sound
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);

                // Compute total enthalpy
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in y-direction:
                // G(ρ) = ρv, G(ρu) = ρuv, G(ρv) = ρv² + p, G(ρe) = ρHv
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;
            
                // Compute Roe-averaged variables
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;

                // Compute approximate wave speeds using Roe-averaged variables,
                // following Einfeldt's approach for y-direction
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;
                
                // Einfeldt's correction term for y-direction wave speeds
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R / (sum_sqrt * sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (v_R - v_L) * (v_R - v_L));
                
                // Check for valid wave speeds
                if (isnan(sqrt_d_sq)) {
                    // Fall back to simple Roe average
                    const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * (ubar_sq + vbar_sq)));
                    sqrt_d_sq = cbar;
                }
                
                // Wave speeds for y-direction: based on vbar instead of ubar
                const double sleft = vbar - sqrt_d_sq;
                const double sright = vbar + sqrt_d_sq;
            
                // Compute intermediate wave speed and pressure for y-direction
                double rho_G, rhou_G, rhov_G, rhoe_G;
                
                // Handle potential division by zero
                const double denom = rho_L * (sleft - v_L) - rho_R * (sright - v_R);
                
                if (fabs(denom) < 1e-12) {
                    // Fall back to HLL flux
                    const double inv_sdiff = 1.0 / (sright - sleft);
                    const double sleft_sright = sleft * sright;
                    
                    rho_G = (sright * rho_GL - sleft * rho_GR + sleft_sright * (rho_R - rho_L)) * inv_sdiff;
                    rhou_G = (sright * rhou_GL - sleft * rhou_GR + sleft_sright * (rhou_R - rhou_L)) * inv_sdiff;
                    rhov_G = (sright * rhov_GL - sleft * rhov_GR + sleft_sright * (rhov_R - rhov_L)) * inv_sdiff;
                    rhoe_G = (sright * rhoe_GL - sleft * rhoe_GR + sleft_sright * (rhoe_R - rhoe_L)) * inv_sdiff;
                }
                else {
                    // Compute contact discontinuity speed for y-direction
                    const double s_star = (p_R - p_L + rho_L * v_L * (sleft - v_L) 
                                         - rho_R * v_R * (sright - v_R)) / denom;
                    
                    const double p_star = p_L + rho_L * (sleft - v_L) * (s_star - v_L);
                    const double v_stat_L = rho_L * (sleft - v_L) / (sleft - s_star);
                    const double v_stat_R = rho_R * (sright - v_R) / (sright - s_star);

                    const double fe_star_L = rhoe_L * inv_rho_L + (s_star - v_L) 
                                           * (s_star + p_L / (rho_L * (sleft - v_L)));
                    const double fe_star_R = rhoe_R * inv_rho_R + (s_star - v_R) 
                                           * (s_star + p_R / (rho_R * (sright - v_R)));
                    
                    // Determine HLLC numerical flux based on wave speeds
                    if (sleft >= 0.0) {
                        rho_G = rho_GL;
                        rhou_G = rhou_GL;
                        rhov_G = rhov_GL;
                        rhoe_G = rhoe_GL;
                    }
                    else if (sright <= 0.0) {
                        rho_G = rho_GR;
                        rhou_G = rhou_GR;
                        rhov_G = rhov_GR;
                        rhoe_G = rhoe_GR;
                    }
                    else if (s_star >= 0.0) {
                        // Left and star region
                        rho_G = rho_GL + sleft * (v_stat_L - rho_L);
                        rhou_G = rhou_GL + sleft * (v_stat_L * u_L - rhou_L);
                        rhov_G = rhov_GL + sleft * (v_stat_L * s_star - rhov_L);
                        rhoe_G = rhoe_GL + sleft * (v_stat_L * fe_star_L - rhoe_L);
                    }
                    else {
                        // Star and right region
                        rho_G = rho_GR + sright * (v_stat_R - rho_R);
                        rhou_G = rhou_GR + sright * (v_stat_R * u_R - rhou_R);
                        rhov_G = rhov_GR + sright * (v_stat_R * s_star - rhov_R);
                        rhoe_G = rhoe_GR + sright * (v_stat_R * fe_star_R - rhoe_R);
                    }
                }

                // Store computed fluxes
                Flux[0][j][k] = rho_G; 
                Flux[1][j][k] = rhou_G; 
                Flux[2][j][k] = rhov_G;
                Flux[3][j][k] = rhoe_G;
            }
        }
    }
}


/**
 * Roe Flux Calculator for Euler Equations (2D)
 * Computes numerical fluxes using Roe's approximate Riemann solver with entropy fix
 * Roe solver provides exact resolution of isolated shock and contact discontinuities
 * 
 * @param dir Direction: 1 for x-direction fluxes, 2 for y-direction fluxes
 * @param var Number of variables (currently fixed at 4)
 * @param rows Number of grid points in y-direction
 * @param cols Number of grid points in x-direction
 * @param GC Number of ghost cells
 * @param L_state Left state conservative variables at cell interfaces
 * @param R_state Right state conservative variables at cell interfaces
 * @param flux Output flux arrays
 */
static inline void Roe_Flux(int dir, int var, int rows, int cols, int GC, 
                            double (*L_state)[rows][cols], 
                            double (*R_state)[rows][cols], 
                            double (*flux)[rows][cols]) {
    
    // Parameters for entropy fix (prevents expansion shocks)
    const double epsilon = 1e-6;
    const double epsilon2 = epsilon * epsilon;
    
    // Precompute frequently used constants
    const double gamma_minus_1 = M_gamma - 1.0;
    const double gamma_ratio = M_gamma / gamma_minus_1;  // γ/(γ-1)
    
    if (dir == 1) { // x-direction flux calculation
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j < rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables using division optimization
                // Compute reciprocal once to avoid multiple divisions
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;      // x-velocity left state
                const double u_R = rhou_R * inv_rho_R;      // x-velocity right state
                const double v_L = rhov_L * inv_rho_L;      // y-velocity left state
                const double v_R = rhov_R * inv_rho_R;      // y-velocity right state
                
                // Compute velocity squares for kinetic energy calculation
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law: p = (ρe - 0.5ρ(u²+v²)) * (γ-1)
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;

                 // Compute speed of sound
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);


                // Compute total enthalpy: H = 0.5(u²+v²) + γ/(γ-1) * p/ρ
                // Enthalpy is conserved across shocks and used in flux calculation
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in x-direction:
                // F(ρ) = ρu, F(ρu) = ρu² + p, F(ρv) = ρuv, F(ρe) = ρHu
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;

                // Compute Roe-averaged variables for the intermediate state
                // Roe averages satisfy the "property U" (conservation and consistency)
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double rhobar = sqrt_rho_L * sqrt_rho_R;          // Roe-averaged density
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;  // Roe-averaged x-velocity
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;  // Roe-averaged y-velocity
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;  // Roe-averaged enthalpy
                const double q2bar = ubar * ubar + vbar * vbar;         // Roe-averaged velocity magnitude squared
                const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * q2bar));  // Roe-averaged speed of sound

                // Compute Einfeldt's correction term for better wave speed estimation
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R / (sum_sqrt * sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (u_R - u_L) * (u_R - u_L));


                // Compute eigenvalues of the Roe matrix
                // For x-direction: λ₁ = ubar-cbar, λ₂ = ubar, λ₃ = ubar, λ₄ = ubar+cbar
                double eigenvalues[4];
                eigenvalues[0] = ubar - sqrt_d_sq;
                eigenvalues[1] = ubar;
                eigenvalues[2] = ubar;
                eigenvalues[3] = ubar + sqrt_d_sq;

                // Apply entropy fix (Harten's entropy correction)
                // Prevents expansion shocks and ensures entropy condition is satisfied
                const double two_epsilon = 2.0 * epsilon;
                for (int i = 0; i < 4; i++) {
                    const double abs_lambda = fabs(eigenvalues[i]);
                    if (abs_lambda < epsilon) {
                        // Smooth transition near zero eigenvalues
                        eigenvalues[i] = (eigenvalues[i] * eigenvalues[i] + epsilon2) / two_epsilon;
                    } else {
                        eigenvalues[i] = abs_lambda;
                    }
                }

                // Compute wave strengths (α coefficients) for characteristic decomposition
                const double drho = rho_R - rho_L;      // Density jump
                const double du = u_R - u_L;            // x-velocity jump
                const double dv = v_R - v_L;            // y-velocity jump
                const double dp = p_R - p_L;            // Pressure jump
                
                const double cbar_sq = cbar * cbar;
                const double inv_2cbar_sq = 1.0 / (2.0 * cbar_sq);
                const double inv_cbar_sq = 1.0 / cbar_sq;

                // Wave strengths (α) from the jump in characteristic variables
                const double alpha1 = (dp - rhobar * cbar * du) * inv_2cbar_sq;  // Left acoustic wave
                const double alpha2 = drho - dp * inv_cbar_sq;                   // Entropy wave
                const double alpha3 = rhobar * dv;                               // Shear wave
                const double alpha4 = (dp + rhobar * cbar * du) * inv_2cbar_sq;  // Right acoustic wave

                // Right eigenvectors of the Roe matrix for x-direction
                // These form a basis for the characteristic decomposition
                const double eigenvectors[4][4] = {
                    {1.0, ubar - cbar, vbar, Hbar - ubar * cbar},  // Left acoustic wave
                    {1.0, ubar, vbar, 0.5 * q2bar},                // Entropy wave
                    {0.0, 0.0, 1.0, vbar},                         // Shear wave
                    {1.0, ubar + cbar, vbar, Hbar + ubar * cbar}   // Right acoustic wave
                };

                // Compute Roe numerical flux: F_Roe = 0.5*(F_L + F_R) - 0.5*Σ|λ|αR
                // Central part (average of left and right fluxes)
                double rho_F = 0.5 * (rho_FL + rho_FR);
                double rhou_F = 0.5 * (rhou_FL + rhou_FR);
                double rhov_F = 0.5 * (rhov_FL + rhov_FR);
                double rhoe_F = 0.5 * (rhoe_FL + rhoe_FR);

                // Add dissipation term (characteristic-based upwinding)
                const double wave_strengths[4] = {alpha1, alpha2, alpha3, alpha4};
                
                // Loop over all characteristic waves
                for (int i = 0; i < 4; i++) {
                    const double term = 0.5 * eigenvalues[i] * wave_strengths[i];
                    rho_F -= term * eigenvectors[i][0];
                    rhou_F -= term * eigenvectors[i][1];
                    rhov_F -= term * eigenvectors[i][2];
                    rhoe_F -= term * eigenvectors[i][3];
                }

                // Store computed fluxes in output array
                flux[0][j][k] = rho_F;
                flux[1][j][k] = rhou_F;
                flux[2][j][k] = rhov_F;
                flux[3][j][k] = rhoe_F;
            }
        }
    }
    else if (dir == 2) { // y-direction flux calculation
        // Similar to x-direction but with velocities and eigenvectors rotated
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                // Velocity squares for kinetic energy
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;

                 // Compute speed of sound
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);


                // Compute total enthalpy
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in y-direction:
                // G(ρ) = ρv, G(ρu) = ρuv, G(ρv) = ρv² + p, G(ρe) = ρHv
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;

                // Compute Roe-averaged variables for y-direction
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double rhobar = sqrt_rho_L * sqrt_rho_R;
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;
                const double q2bar = ubar * ubar + vbar * vbar;
                const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * q2bar));


                // Compute Einfeldt's correction term for better wave speed estimation
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R / (sum_sqrt * sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (v_R - v_L) * (v_R - v_L));

                // Compute eigenvalues for y-direction
                // For y-direction: λ₁ = vbar-cbar, λ₂ = vbar, λ₃ = vbar, λ₄ = vbar+cbar
                double eigenvalues[4];
                eigenvalues[0] = vbar - sqrt_d_sq;
                eigenvalues[1] = vbar;
                eigenvalues[2] = vbar;
                eigenvalues[3] = vbar + sqrt_d_sq;

                // Apply entropy fix
                const double two_epsilon = 2.0 * epsilon;
                for (int i = 0; i < 4; i++) {
                    const double abs_lambda = fabs(eigenvalues[i]);
                    if (abs_lambda < epsilon) {
                        eigenvalues[i] = (eigenvalues[i] * eigenvalues[i] + epsilon2) / two_epsilon;
                    } else {
                        eigenvalues[i] = abs_lambda;
                    }
                }

                // Compute wave strengths for y-direction
                const double drho = rho_R - rho_L;
                const double du = u_R - u_L;
                const double dv = v_R - v_L;
                const double dp = p_R - p_L;
                
                const double cbar_sq = cbar * cbar;
                const double inv_2cbar_sq = 1.0 / (2.0 * cbar_sq);
                const double inv_cbar_sq = 1.0 / cbar_sq;

                // Note: du and dv swapped compared to x-direction
                const double alpha1 = (dp - rhobar * cbar * dv) * inv_2cbar_sq;  // Left acoustic wave
                const double alpha2 = drho - dp * inv_cbar_sq;                   // Entropy wave
                const double alpha3 = rhobar * du;                               // Shear wave
                const double alpha4 = (dp + rhobar * cbar * dv) * inv_2cbar_sq;  // Right acoustic wave

                // Right eigenvectors of the Roe matrix for y-direction
                const double eigenvectors[4][4] = {
                    {1.0, ubar, vbar - cbar, Hbar - vbar * cbar},  // Left acoustic wave
                    {1.0, ubar, vbar, 0.5 * q2bar},                // Entropy wave
                    {0.0, 1.0, 0.0, ubar},                         // Shear wave
                    {1.0, ubar, vbar + cbar, Hbar + vbar * cbar}   // Right acoustic wave
                };

                // Compute Roe numerical flux for y-direction
                double rho_G = 0.5 * (rho_GL + rho_GR);
                double rhou_G = 0.5 * (rhou_GL + rhou_GR);
                double rhov_G = 0.5 * (rhov_GL + rhov_GR);
                double rhoe_G = 0.5 * (rhoe_GL + rhoe_GR);

                // Add dissipation term
                const double wave_strengths[4] = {alpha1, alpha2, alpha3, alpha4};
                
                for (int i = 0; i < 4; i++) {
                    const double term = 0.5 * eigenvalues[i] * wave_strengths[i];
                    rho_G -= term * eigenvectors[i][0];
                    rhou_G -= term * eigenvectors[i][1];
                    rhov_G -= term * eigenvectors[i][2];
                    rhoe_G -= term * eigenvectors[i][3];
                }

                // Store computed fluxes
                flux[0][j][k] = rho_G;
                flux[1][j][k] = rhou_G;
                flux[2][j][k] = rhov_G;
                flux[3][j][k] = rhoe_G;
            }
        }
    }
}

/*static inline void ER_Flux(int dir, int var, int rows, int cols, int GC,
                           double (*x)[rows][cols],
                           double (*y)[rows][cols],
                           double (*z)[rows][cols]) {
    
    const double gamma_minus_1 = M_gamma - 1.0;
    const double gamma_plus_1 = M_gamma + 1.0;
    const double inv_gamma = 1.0 / M_gamma;
    const double inv_gamma_minus_1 = 1.0 / gamma_minus_1;
    const double two_inv_gamma_plus_1 = 2.0 / gamma_plus_1;
    const double gamma_minus_1_over_2gamma = (M_gamma - 1.0) / (2.0 * M_gamma);
    const double two_gamma_over_gamma_minus_1 = 2.0 * M_gamma / gamma_minus_1;
    const double two_over_gamma_minus_1 = 2.0 / gamma_minus_1;
    const double eps_p = 1e-10;
    const double max_p = 1e6;
    
    if (dir == 1) { // x方向通量
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j <= rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                double rho_F = 0.0, u_F = 0.0, p_F = 0.0;
                
                // 读取左右原始变量 (x方向)
                const double rho_l = x[0][j][k];
                const double rho_r = y[0][j][k];
                const double u_l = x[1][j][k];
                const double u_r = y[1][j][k];
                const double p_l = x[3][j][k];  // 注意：压力通常在第4个变量
                const double p_r = y[3][j][k];
                
                const double c_l = sqrt(M_gamma * p_l / rho_l);
                const double c_r = sqrt(M_gamma * p_r / rho_r);
                
                // 调用精确黎曼求解器
                const double p_star = ExactRiemann_pstar(rho_l, rho_r, u_l, u_r, p_l, p_r, M_gamma, M_gamma, eps_p, max_p);
                const double u_star = ExactRiemann_ustar(p_star, rho_l, rho_r, u_l, u_r, p_l, p_r, M_gamma, M_gamma);
                
                const double si = 0.0; // 界面位置，对于通量计算通常为0

                // 左侧波系处理
                if (p_star >= p_l) { // 左激波
                    const double part1 = gamma_plus_1 * p_star + gamma_minus_1 * p_l;
                    const double part2 = gamma_minus_1 * p_star + gamma_plus_1 * p_l;
                    const double rho_star_L = rho_l * part1 / part2;
                    const double A_L = 2.0 / (gamma_plus_1 * rho_l);
                    const double B_L = gamma_minus_1 * p_l / gamma_plus_1;
                    const double S_L = u_l - (1.0 / rho_l) * sqrt((B_L + p_star) / A_L);
                    
                    if (si <= S_L) {
                        rho_F = rho_l;
                        u_F = u_l;
                        p_F = p_l;
                    }
                    else if (S_L < si && si <= u_star) {
                        rho_F = rho_star_L;
                        u_F = u_star;
                        p_F = p_star;
                    }
                }
                else { // 左稀疏波
                    const double p_ratio_L = p_star / p_l;
                    const double rho_star_L = rho_l * pow(p_ratio_L, inv_gamma);
                    const double aL_star = c_l * pow(p_ratio_L, gamma_minus_1_over_2gamma);
                    const double S_HL = u_l - c_l;
                    const double S_TL = u_star - aL_star;
                    
                    if (si <= S_HL) {
                        rho_F = rho_l;
                        u_F = u_l;
                        p_F = p_l;
                    }
                    else if (si > S_HL && si <= S_TL) {
                        const double temp_L = two_inv_gamma_plus_1 + gamma_minus_1 * (u_l - si) / (gamma_plus_1 * c_l);
                        rho_F = rho_l * pow(temp_L, two_over_gamma_minus_1);
                        u_F = 2.0 * (c_l + gamma_minus_1 * u_l * 0.5 + si) / gamma_plus_1;
                        p_F = p_l * pow(temp_L, two_gamma_over_gamma_minus_1);
                    }
                    else if (si > S_TL && si <= u_star) {
                        rho_F = rho_star_L;
                        u_F = u_star;
                        p_F = p_star;
                    }
                }
                
                // 右侧波系处理
                if (p_star >= p_r) { // 右激波
                    const double part1_r = gamma_plus_1 * p_star + gamma_minus_1 * p_r;
                    const double part2_r = gamma_minus_1 * p_star + gamma_plus_1 * p_r;
                    const double rho_star_r = rho_r * part1_r / part2_r;
                    const double A_r = 2.0 / (gamma_plus_1 * rho_r);
                    const double B_r = gamma_minus_1 * p_r / gamma_plus_1;
                    const double S_r = u_r + (1.0 / rho_r) * sqrt((B_r + p_star) / A_r);
                    
                    if (si > u_star && si <= S_r) {
                        rho_F = rho_star_r;
                        u_F = u_star;
                        p_F = p_star;
                    }
                    else if (si > S_r) {
                        rho_F = rho_r;
                        u_F = u_r;
                        p_F = p_r;
                    }
                }
                else { // 右稀疏波
                    const double p_ratio_R = p_star / p_r;
                    const double rho_star_r = rho_r * pow(p_ratio_R, inv_gamma);
                    const double aR_star = c_r * pow(p_ratio_R, gamma_minus_1_over_2gamma);
                    const double S_HR = u_r + c_r;
                    const double S_TR = u_star + aR_star;
                    
                    if (si > u_star && si <= S_TR) {
                        rho_F = rho_star_r;
                        u_F = u_star;
                        p_F = p_star;
                    }
                    else if (si > S_TR && si <= S_HR) {
                        const double temp_R = two_inv_gamma_plus_1 - gamma_minus_1 * (u_r - si) / (gamma_plus_1 * c_r);
                        rho_F = rho_r * pow(temp_R, two_over_gamma_minus_1);
                        u_F = 2.0 * (-c_r + gamma_minus_1 * u_r * 0.5 + si) / gamma_plus_1;
                        p_F = p_r * pow(temp_R, two_gamma_over_gamma_minus_1);
                    }
                    else if (si > S_HR) {
                        rho_F = rho_r;
                        u_F = u_r;
                        p_F = p_r;
                    }
                }
                
                // 计算能量修正项
                const double rhobar = sqrt(rho_l * rho_r);
                const double T_R = p_r / rho_r;
                const double T_L = p_l / rho_l;
                const double s_L_P = u_l - c_l;
                const double s_R_P = u_r + c_r;
                const double denom = s_R_P - s_L_P;
                const double energy_correction = (M_gamma * inv_gamma_minus_1) * rhobar * (T_R - T_L) * fabs(s_L_P * s_R_P / denom);
                
                // 存储x方向通量
                z[0][j][k] = rho_F * u_F;  // 质量通量
                z[1][j][k] = rho_F * u_F * u_F + p_F;  // x动量通量
                z[2][j][k] = 0.0;  // y动量通量（x方向为0）
                z[3][j][k] = ((0.5 * rho_F * u_F * u_F + p_F * inv_gamma_minus_1) + p_F) * u_F - energy_correction;  // 能量通量
            }
        }
    }
    else if (dir == 2) { // y方向通量
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                double rho_G = 0.0, v_G = 0.0, p_G = 0.0;
                
                // 读取上下原始变量 (y方向)
                const double rho_l = x[0][j][k];
                const double rho_r = y[0][j][k];
                const double v_l = x[2][j][k];  // 注意：y方向使用v分量
                const double v_r = y[2][j][k];
                const double p_l = x[3][j][k];
                const double p_r = y[3][j][k];
                
                const double c_l = sqrt(M_gamma * p_l / rho_l);
                const double c_r = sqrt(M_gamma * p_r / rho_r);
                
                // 调用精确黎曼求解器（y方向使用v分量）
                const double p_star = ExactRiemann_pstar(rho_l, rho_r, v_l, v_r, p_l, p_r, M_gamma, M_gamma, eps_p, max_p);
                const double v_star = ExactRiemann_ustar(p_star, rho_l, rho_r, v_l, v_r, p_l, p_r, M_gamma, M_gamma);
                
                const double si = 0.0; // 界面位置，对于通量计算通常为0

                // 下侧波系处理（对应y方向的"左"侧）
                if (p_star >= p_l) { // 下激波
                    const double part1 = gamma_plus_1 * p_star + gamma_minus_1 * p_l;
                    const double part2 = gamma_minus_1 * p_star + gamma_plus_1 * p_l;
                    const double rho_star_L = rho_l * part1 / part2;
                    const double A_L = 2.0 / (gamma_plus_1 * rho_l);
                    const double B_L = gamma_minus_1 * p_l / gamma_plus_1;
                    const double S_L = v_l - (1.0 / rho_l) * sqrt((B_L + p_star) / A_L);
                    
                    if (si <= S_L) {
                        rho_G = rho_l;
                        v_G = v_l;
                        p_G = p_l;
                    }
                    else if (S_L < si && si <= v_star) {
                        rho_G = rho_star_L;
                        v_G = v_star;
                        p_G = p_star;
                    }
                }
                else { // 下稀疏波
                    const double p_ratio_L = p_star / p_l;
                    const double rho_star_L = rho_l * pow(p_ratio_L, inv_gamma);
                    const double aL_star = c_l * pow(p_ratio_L, gamma_minus_1_over_2gamma);
                    const double S_HL = v_l - c_l;
                    const double S_TL = v_star - aL_star;
                    
                    if (si <= S_HL) {
                        rho_G = rho_l;
                        v_G = v_l;
                        p_G = p_l;
                    }
                    else if (si > S_HL && si <= S_TL) {
                        const double temp_L = two_inv_gamma_plus_1 + gamma_minus_1 * (v_l - si) / (gamma_plus_1 * c_l);
                        rho_G = rho_l * pow(temp_L, two_over_gamma_minus_1);
                        v_G = 2.0 * (c_l + gamma_minus_1 * v_l * 0.5 + si) / gamma_plus_1;
                        p_G = p_l * pow(temp_L, two_gamma_over_gamma_minus_1);
                    }
                    else if (si > S_TL && si <= v_star) {
                        rho_G = rho_star_L;
                        v_G = v_star;
                        p_G = p_star;
                    }
                }
                
                // 上侧波系处理（对应y方向的"右"侧）
                if (p_star >= p_r) { // 上激波
                    const double part1_r = gamma_plus_1 * p_star + gamma_minus_1 * p_r;
                    const double part2_r = gamma_minus_1 * p_star + gamma_plus_1 * p_r;
                    const double rho_star_r = rho_r * part1_r / part2_r;
                    const double A_r = 2.0 / (gamma_plus_1 * rho_r);
                    const double B_r = gamma_minus_1 * p_r / gamma_plus_1;
                    const double S_r = v_r + (1.0 / rho_r) * sqrt((B_r + p_star) / A_r);
                    
                    if (si > v_star && si <= S_r) {
                        rho_G = rho_star_r;
                        v_G = v_star;
                        p_G = p_star;
                    }
                    else if (si > S_r) {
                        rho_G = rho_r;
                        v_G = v_r;
                        p_G = p_r;
                    }
                }
                else { // 上稀疏波
                    const double p_ratio_R = p_star / p_r;
                    const double rho_star_r = rho_r * pow(p_ratio_R, inv_gamma);
                    const double aR_star = c_r * pow(p_ratio_R, gamma_minus_1_over_2gamma);
                    const double S_HR = v_r + c_r;
                    const double S_TR = v_star + aR_star;
                    
                    if (si > v_star && si <= S_TR) {
                        rho_G = rho_star_r;
                        v_G = v_star;
                        p_G = p_star;
                    }
                    else if (si > S_TR && si <= S_HR) {
                        const double temp_R = two_inv_gamma_plus_1 - gamma_minus_1 * (v_r - si) / (gamma_plus_1 * c_r);
                        rho_G = rho_r * pow(temp_R, two_over_gamma_minus_1);
                        v_G = 2.0 * (-c_r + gamma_minus_1 * v_r * 0.5 + si) / gamma_plus_1;
                        p_G = p_r * pow(temp_R, two_gamma_over_gamma_minus_1);
                    }
                    else if (si > S_HR) {
                        rho_G = rho_r;
                        v_G = v_r;
                        p_G = p_r;
                    }
                }
                
                // 计算能量修正项（y方向）
                const double rhobar = sqrt(rho_l * rho_r);
                const double T_R = p_r / rho_r;
                const double T_L = p_l / rho_l;
                const double s_L_P = v_l - c_l;
                const double s_R_P = v_r + c_r;
                const double denom = s_R_P - s_L_P;
                const double energy_correction = (M_gamma * inv_gamma_minus_1) * rhobar * (T_R - T_L) * fabs(s_L_P * s_R_P / denom);
                
                // 存储y方向通量
                z[0][j][k] = rho_G * v_G;  // 质量通量
                z[1][j][k] = 0.0;  // x动量通量（y方向为0）
                z[2][j][k] = rho_G * v_G * v_G + p_G;  // y动量通量
                z[3][j][k] = ((0.5 * rho_G * v_G * v_G + p_G * inv_gamma_minus_1) + p_G) * v_G - energy_correction;  // 能量通量
            }
        }
    }
}*/

static inline void ER_Flux(int var, int rows, int GC, double (*x)[rows], double (*y)[rows] ,double (*z)[rows]) {
    
    for (int j = GC-1; j <= rows-GC; j++) {
        double rho_F = 0,u_F=0,p_F=0;
        double rho_l = x[0][j];
        double rho_r = y[0][j];
        double u_l = x[1][j];
        double u_r = y[1][j];
        double p_l = x[2][j];
        double p_r = y[2][j];
        double c_l = sqrt(M_gamma*p_l/rho_l);
        double c_r = sqrt(M_gamma*p_r/rho_r);
        double p_star = ExactRieamnna_pstar(rho_l, rho_r, u_l, u_r, p_l, p_r, M_gamma, M_gamma ,1e-10, 1e6);
		double u_star = ExactRieamnna_ustar(p_star, rho_l, rho_r, u_l,u_r, p_l,p_r, M_gamma,M_gamma);
        
         
        //分别考虑左激波，左稀疏波的情况
        if (p_star >=  p_l){
            double part1 = (M_gamma + 1)* p_star + (M_gamma - 1) * p_l ;
            double part2 = (M_gamma - 1) * p_star + (M_gamma + 1) * p_l;
            double rho_star_L = rho_l * part1 / part2;
            double A_L = 2.0 / ((M_gamma + 1) * rho_l);
            //计算激波的速度
            double B_L = (M_gamma - 1) * p_l  / (M_gamma + 1);
            double S_L = u_l - (1 / rho_l) * sqrt((B_L + p_star) / A_L); 
            double si = 0;
            if (si <= S_L){
                rho_F = rho_l;
                u_F = u_l;
                p_F = p_l;
            }
            else if (S_L < si && si <= u_star)
            {
                rho_F = rho_star_L;
                u_F = u_star;
                p_F = p_star;
            }
        }
        else {
            //计算接触间断左侧密度
            double part1 = p_star;
            double part2 = p_l;
            double rho_star_L = rho_l * pow(part1/part2, 1/M_gamma);
            //计算左稀疏波波头和波尾的速度
            double aL_star = c_l * pow((p_star / p_l), ((M_gamma - 1) / (2 *M_gamma)));
            double S_HL = u_l - c_l;
            double S_TL = u_star - aL_star;
            double si=0;
            if (si <= S_HL){
                rho_F = rho_l;
                u_F = u_l;
                p_F = p_l;
            }
            else if (si <= S_TL && si > S_HL)
            {
                rho_F = rho_l * pow(2/(M_gamma + 1) + (M_gamma - 1) * (u_l - si) / ((M_gamma + 1) * c_l), (2 / (M_gamma - 1))); 
                u_F = 2 * (c_l + (M_gamma - 1) * u_l/2 + si) / (M_gamma + 1);
                p_F = p_l * pow(2 / (M_gamma + 1) + (M_gamma - 1) * (u_l - si) / ((M_gamma + 1) * c_l),  2 * M_gamma / (M_gamma - 1));
            }
            else if (si <= u_star && si > S_TL)
            {
                rho_F = rho_star_L;
                u_F = u_star;
                p_F = p_star;
            }
        }

        //分别考虑右稀疏波和右激波的情况
        if ( p_star >=  p_r){
            double part1 =  (M_gamma + 1) * p_star + (M_gamma - 1) * p_r;
            double part2 =  (M_gamma - 1) * p_star + (M_gamma + 1) * p_r;
            double rho_star_r = rho_r * part1 / part2;
            double A_r = 2 / ((M_gamma + 1) * rho_r);
            //计算激波的速度
            double B_r = (M_gamma - 1) * p_r / (M_gamma + 1);
            double S_r = u_r + (1 / rho_r) * sqrt((B_r + p_star) / A_r); 
            double si = 0;
            if (si <= S_r &&  u_star < si)
            {
                rho_F = rho_star_r;
                u_F = u_star;
                p_F = p_star;
            }
            else if (S_r < 0){
                rho_F = rho_r;
                u_F = u_r;
                p_F = p_r;
            }
        }
        else{
            //计算接触间断左侧密度
            double part1 = p_star;
            double part2 = p_r;
            double rho_star_r = rho_r * pow(part1 / part2, 1/M_gamma);
            //计算左稀疏波波头和波尾的速度
            double aR_star = c_r * pow((p_star / p_r), ((M_gamma - 1) / (2 *M_gamma)));
            double S_HR = u_r + c_r;
            double S_TR = u_star + aR_star;   
            double si=0;
            if (si <= S_TR && si > u_star)
            {
                rho_F = rho_star_r;
                u_F = u_star;
                p_F = p_star;
            }    
            else if (si<= S_HR && si >S_TR)
            {
                rho_F = rho_r * pow(2/(M_gamma + 1) - (M_gamma - 1) * (u_r - si) / ((M_gamma + 1) * c_r), (2 / (M_gamma - 1))); 
                u_F = 2 * (-c_r + (M_gamma - 1) * u_r / 2 + si) / (M_gamma + 1);
                p_F = p_r * pow(2 / (M_gamma + 1) - (M_gamma - 1) * (u_r - si) / ((M_gamma + 1) * c_r),  2 * M_gamma / (M_gamma - 1));
            }
            else if (si > S_HR){
                rho_F = rho_r;
                u_F = u_r;
                p_F = p_r;
            }
        }

        
        double rhobar = sqrt(rho_l*rho_r);
        double T_R = p_r/rho_r;
        double T_L = p_l/rho_l;
        double s_L_P = u_l - c_l;
        double s_R_P = u_r + c_r;
        z[0][j] = rho_F * u_F; 
        z[1][j] = rho_F * u_F *u_F + p_F; 
        z[2][j] = ((0.5*rho_F * u_F *u_F + p_F/(M_gamma-1)) + p_F) * u_F - ((M_gamma)/(M_gamma-1)) * rhobar *  (T_R - T_L ) * (fabs(s_L_P * s_R_P/(s_R_P-s_L_P))); 

        if (j == 200) {
            //double p_diff = p_star - 2.926650;
            //p_star_accumulator += p_diff;
            //printf("p_Refect = %f, p_Refect - p_Exact = %f, 累计差值 = %f.\n",p_star, p_diff, p_star_accumulator);
            //int k = 0;
            //z[1][j] = p_star;

            Flux_test += rho_F * u_F;
            printf("累加通量 = %f.\n", Flux_test);
        }
    }   
}
//这里是计算精确Riemann解

double F_K_P(double p, double p_refer, double rho_refer, double gamma,double c_refer){
	//计算f_K(p) K=L/R
    //:param p: 变量 求解过程中的压强
    //:param _refer: 固定量要么是L 要么是R
    //:return: f_K 函数值
	if (p >=  p_refer)
	{
		double A_K = 2 / ((gamma + 1) * rho_refer);
        double B_K = (gamma - 1) * (p_refer) / (gamma + 1);
        return (p - p_refer) * sqrt(A_K / (B_K + p ));
	}
	else
		return (2 * c_refer / (gamma - 1)) * (pow(p/p_refer,(gamma - 1) / (2 *gamma)) - 1);
}

double dF_K_P(double p, double p_refer, double rho_refer, double gamma,double c_refer){
    //计算f'_K(p) K=L/R 它是 f_K(p) 关于 p 的导数
    //:param p: 变量 求解过程中的压强
    //:param _refer: 固定量要么是L 要么是R
	//:return:f'_K(p) 函数值
	if (p >=  p_refer)
	{
		double A_K = 2 / ((gamma + 1) * rho_refer);
        double B_K = (gamma - 1) * (p_refer) / (gamma + 1);
		double part1 = sqrt((A_K) / (B_K + p ));
        double part2 = 1 - (p - p_refer) / (2 * (B_K + p));
        return part1 * part2;

	}
	else{
		double part1 = pow(p/p_refer, -(gamma + 1) / (2 * gamma));
        double part2 = 1 / (rho_refer * c_refer);
        return part1 * part2;
	}
}


// 定义一个函数来计算f(x)
double ExactRieamnna_pstar(double rhol, double rhor, double ul,double ur, double pl, \
							double pr, double gammal,double gammar ,double tol, double maxit){
	
	//利用状态方方程计算声速度
	double cl = sqrt(gammal * pl/rhol);
	double cr = sqrt(gammar * pr/rhor);
	//计算迭代初始值
    double p = max_of_two(tol, 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rhor + rhol) * (cl+cr));
    double p_new = 0.0;
	double f_p = 0; 
	double df_p = 0;

	for ( int i = 0; i < maxit; i++)
	{
		//迭代初始值问题
		if (p < 0.0)
		{
			printf("negative pressure\n");
			exit(1); 
		}

		//计算Fpp的函数
        f_p = F_K_P(p, pl,rhol,gammal,cl) + F_K_P(p, pr,rhor,gammar,cr) + ur - ul;
        df_p = dF_K_P(p, pl,rhol,gammal,cl) + dF_K_P(p, pr,rhor,gammar,cr);
        p_new = p - f_p / df_p;

		if(p_new < 0.0 ){ 
            printf("negative pressure\n");
			exit(1); 
		}
		if (2 * fabs(p_new - p) / (p + p_new) < tol){
			p = p_new;
            break;
		}
		else{
			p = p_new;
		}
	}


	double p_star = p;
	return p_star;
}
double ExactRieamnna_ustar(double p_star, double rhol, double rhor, double ul,double ur, double pl, \
							double pr, double gammal,double gammar){
	
	//利用状态方方程计算声速度
	double cl = sqrt(gammal * pl/rhol);
	double cr = sqrt(gammar * pr/rhor);
	//计算迭代初始值
	double u_star = 0.5 * (ul + ur) + 0.5 * (F_K_P(p_star, pr,rhor,gammar,cr) - F_K_P(p_star, pl,rhol,gammal,cl));
    return u_star;
}


/*                               ************************************                               */
/*                               ************************************                               */
/*                              Riemann Solver：Exact and Approximate                               */
/*                               ************************************                               */
/*                               ************************************                               */

/*                                      ******************                                          */
/*                                  具有人工热传导的Riemann Solver                                   */
/*                                      ******************                                          */
/**
 * HLLHC (Harten-Lax-van Leer with Heat Capacity correction) Flux Calculator for Euler Equations (2D)
 * Computes numerical fluxes using HLL approximate Riemann solver with heat capacity correction term
 * HLLHC adds an energy correction term to improve accuracy for flows with significant heat transfer
 * 
 * @param dir Direction: 1 for x-direction fluxes, 2 for y-direction fluxes
 * @param var Number of variables (currently fixed at 4)
 * @param rows Number of grid points in y-direction
 * @param cols Number of grid points in x-direction
 * @param GC Number of ghost cells
 * @param L_state Left state conservative variables at cell interfaces
 * @param R_state Right state conservative variables at cell interfaces
 * @param flux Output flux arrays
 */
static inline void HLLHC_Flux(int dir, int var, int rows, int cols, int GC, 
                              double (*L_state)[rows][cols], 
                              double (*R_state)[rows][cols], 
                              double (*flux)[rows][cols]) {
    
    // Precompute frequently used constants
    const double gamma_minus_1 = M_gamma - 1.0;
    const double gamma_ratio = M_gamma / gamma_minus_1;      // γ/(γ-1) for enthalpy
    const double inv_gamma_minus_1 = 1.0 / gamma_minus_1;    // 1/(γ-1) for internal energy
    
    if (dir == 1) {
        // X-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j < rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables using division optimization
                // Compute reciprocal once to avoid multiple divisions
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;      // x-velocity left state
                const double u_R = rhou_R * inv_rho_R;      // x-velocity right state
                const double v_L = rhov_L * inv_rho_L;      // y-velocity left state
                const double v_R = rhov_R * inv_rho_R;      // y-velocity right state
                
                // Compute velocity squares for kinetic energy calculation
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law: p = (ρe - 0.5ρ(u²+v²)) * (γ-1)
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;
                
                // Compute speed of sound: a = sqrt(γ * p / ρ)
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);
                
                // Compute total enthalpy: H = 0.5(u²+v²) + γ/(γ-1) * p/ρ
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in x-direction:
                // F(ρ) = ρu, F(ρu) = ρu² + p, F(ρv) = ρuv, F(ρe) = ρHu
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;
                
                // Compute Roe-averaged variables for wave speed estimation
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double rhobar = sqrt_rho_L * sqrt_rho_R;          // Roe-averaged density
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;  // Roe-averaged x-velocity
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;  // Roe-averaged y-velocity
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;  // Roe-averaged enthalpy

                // Compute approximate wave speeds using Roe-averaged variables,
                // following Einfeldt's approach to Roe eigenvalue estimation
                // This provides more robust wave speed estimates, especially for strong shocks
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;

                // Einfeldt's correction term for better wave speed estimation
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R * inv_sum_sqrt * inv_sum_sqrt;
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (u_R - u_L) * (u_R - u_L));
                
                // Check for valid wave speeds (avoid numerical issues)
                /*if (isnan(sqrt_d_sq)) {
                    // Fall back to simple Roe average if Einfeldt's method fails
                    const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * (ubar_sq + vbar_sq)));
                    sqrt_d_sq = cbar;
                }*/
                
                // Wave speeds: s_left = ubar - â, s_right = ubar + â
                // where â is the modified sound speed from Einfeldt's method
                const double sleft = ubar - sqrt_d_sq;
                const double sright = ubar + sqrt_d_sq;

                // Compute denominator for heat capacity correction term
                //const double denom = rho_L * (sleft - u_L) - rho_R * (sright - u_R);
                const double rho_HLL = (sright*rho_R - sleft*rho_L + rho_FL - rho_FR)/(sright-sleft);
               // const double denom = 0.0;
                
                // Compute internal energy per unit mass: e = p/(ρ(γ-1))
                const double e_L = p_L * inv_rho_L * inv_gamma_minus_1;
                const double e_R = p_R * inv_rho_R * inv_gamma_minus_1;
                
                // Determine HLLHC numerical flux based on wave speeds
                // HLLHC adds a heat capacity correction term to the energy flux
                double rho_F, rhou_F, rhov_F, rhoe_F;
                
                if (sleft >= 0.0) {
                    // Entirely left-going waves: use left flux
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if (sright <= 0.0) {
                    // Entirely right-going waves: use right flux
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }
                else { 
                    // Intermediate region (sleft < 0 && sright > 0): use HLLHC flux formula
                    // HLLHC adds a correction term to energy flux for better heat transfer accuracy
                    const double inv_sdiff = 1.0 / (sright - sleft);
                    const double sleft_sright = sleft * sright;
                    
                    // Standard HLL flux for mass and momentum
                    rho_F = (sright * rho_FL - sleft * rho_FR + sleft_sright * (rho_R - rho_L)) * inv_sdiff;
                    rhou_F = (sright * rhou_FL - sleft * rhou_FR + sleft_sright * (rhou_R - rhou_L)) * inv_sdiff;
                    rhov_F = (sright * rhov_FL - sleft * rhov_FR + sleft_sright * (rhov_R - rhov_L)) * inv_sdiff;
                    
                    // HLLHC energy flux: adds heat capacity correction term (e_R - e_L) * denom
                    rhoe_F = (sright * rhoe_FL - sleft * rhoe_FR + sleft_sright * (rhoe_R - rhoe_L)) * inv_sdiff 
                             + (sleft_sright * inv_sdiff) * (e_R - e_L) * rho_HLL;
                }

                // Store computed fluxes in output array
                flux[0][j][k] = rho_F; 
                flux[1][j][k] = rhou_F; 
                flux[2][j][k] = rhov_F;
                flux[3][j][k] = rhoe_F; 
            }
        }   
    }
    else if (dir == 2) {
        // Y-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;

                // Compute speed of sound: a = sqrt(γ * p / ρ)
                const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                const double a_R = sqrt(M_gamma * p_R * inv_rho_R);

                // Compute total enthalpy
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in y-direction:
                // G(ρ) = ρv, G(ρu) = ρuv, G(ρv) = ρv² + p, G(ρe) = ρHv
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;
            
                // Compute Roe-averaged variables
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double rhobar = sqrt_rho_L * sqrt_rho_R;
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;

                // Compute approximate wave speeds using Roe-averaged variables,
                // following Einfeldt's approach for y-direction
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;

                // Einfeldt's correction term for y-direction wave speeds
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R * inv_sum_sqrt * inv_sum_sqrt;
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (v_R - v_L) * (v_R - v_L));
                
                // Check for valid wave speeds
                if (isnan(sqrt_d_sq)) {
                    // Fall back to simple Roe average
                    const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * (ubar_sq + vbar_sq)));
                    sqrt_d_sq = cbar;
                }
                
                // Wave speeds for y-direction: based on vbar instead of ubar
                const double sleft = vbar - sqrt_d_sq;
                const double sright = vbar + sqrt_d_sq;

                // Compute denominator for heat capacity correction term (y-direction)
                //const double denom = rho_L * (sleft - v_L) - rho_R * (sright - v_R);
                const double rho_HLL = (sright*rho_R - sleft*rho_L + rho_GL - rho_GR)/(sright-sleft);
               //const double denom = 0.0;
                
                // Compute internal energy per unit mass
                const double e_L = p_L * inv_rho_L * inv_gamma_minus_1;
                const double e_R = p_R * inv_rho_R * inv_gamma_minus_1;

                // Determine HLLHC numerical flux based on wave speeds
                double rho_G, rhou_G, rhov_G, rhoe_G;
                
                if (sleft >= 0.0) {
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if (sright <= 0.0) {
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }
                else { 
                    // Intermediate region: use HLLHC flux formula
                    const double inv_sdiff = 1.0 / (sright - sleft);
                    const double sleft_sright = sleft * sright;
                    
                    rho_G = (sright * rho_GL - sleft * rho_GR + sleft_sright * (rho_R - rho_L)) * inv_sdiff;
                    rhou_G = (sright * rhou_GL - sleft * rhou_GR + sleft_sright * (rhou_R - rhou_L)) * inv_sdiff;
                    rhov_G = (sright * rhov_GL - sleft * rhov_GR + sleft_sright * (rhov_R - rhov_L)) * inv_sdiff;
                    rhoe_G = (sright * rhoe_GL - sleft * rhoe_GR + sleft_sright * (rhoe_R - rhoe_L)) * inv_sdiff 
                             + (sleft_sright * inv_sdiff) * (e_R - e_L) * rho_HLL;
                }

                // Store computed fluxes
                flux[0][j][k] = rho_G; 
                flux[1][j][k] = rhou_G; 
                flux[2][j][k] = rhov_G;
                flux[3][j][k] = rhoe_G; 
            }
        }
    }
}


/**
 * HLLCHC (Harten-Lax-van Leer Contact with Heat Capacity correction) Flux Calculator for Euler Equations (2D)
 * Computes numerical fluxes using HLLC approximate Riemann solver with heat capacity correction term
 * Combines the contact discontinuity capturing ability of HLLC with improved energy accuracy
 * 
 * @param dir Direction: 1 for x-direction fluxes, 2 for y-direction fluxes
 * @param var Number of variables (currently fixed at 4)
 * @param rows Number of grid points in y-direction
 * @param cols Number of grid points in x-direction
 * @param GC Number of ghost cells
 * @param L_state Left state conservative variables at cell interfaces
 * @param R_state Right state conservative variables at cell interfaces
 * @param flux Output flux arrays
 */
static inline void HLLCHC_Flux(int dir, int var, int rows, int cols, int GC, 
                               double (*L_state)[rows][cols], 
                               double (*R_state)[rows][cols], 
                               double (*flux)[rows][cols]) {
    
    // Precompute frequently used constants
    const double gamma_minus_1 = M_gamma - 1.0;
    const double gamma_ratio = M_gamma / gamma_minus_1;      // γ/(γ-1) for enthalpy
    const double inv_gamma_minus_1 = 1.0 / gamma_minus_1;    // 1/(γ-1) for internal energy
    
    if (dir == 1) {
        // X-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j < rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables using division optimization
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;      // x-velocity left state
                const double u_R = rhou_R * inv_rho_R;      // x-velocity right state
                const double v_L = rhov_L * inv_rho_L;      // y-velocity left state
                const double v_R = rhov_R * inv_rho_R;      // y-velocity right state
                
                // Compute velocity squares for kinetic energy calculation
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law: p = (ρe - 0.5ρ(u²+v²)) * (γ-1)
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;
                
                // Compute speed of sound: a = sqrt(γ * p / ρ)
                //const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                //const double a_R = sqrt(M_gamma * p_R * inv_rho_R);
                
                // Compute total enthalpy: H = 0.5(u²+v²) + γ/(γ-1) * p/ρ
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                //Calculate The Left And Right Sound Speed Values
                const double a_L = sqrt((M_gamma - 1.0) * (H_L - 0.5 * (u_L_sq + v_L_sq)));
                const double a_R = sqrt((M_gamma - 1.0) * (H_R - 0.5 * (u_R_sq + v_R_sq)));
                
                // Compute Euler fluxes in x-direction:
                // F(ρ) = ρu, F(ρu) = ρu² + p, F(ρv) = ρuv, F(ρe) = ρHu
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;
                
                // Compute Roe-averaged variables for intermediate state estimation
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;  // Roe-averaged x-velocity
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;  // Roe-averaged y-velocity
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;  // Roe-averaged enthalpy
                
                // Compute approximate wave speeds using Roe-averaged variables,
                // following Einfeldt's approach to Roe eigenvalue estimation
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;

                // Einfeldt's correction term for better wave speed estimation
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R / (sum_sqrt * sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (u_R - u_L) * (u_R - u_L));
                
                // Check for valid wave speeds (avoid numerical issues)
                if (isnan(sqrt_d_sq)) {
                    // Fall back to simple Roe average if Einfeldt's method fails
                    const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * (ubar_sq + vbar_sq)));
                    sqrt_d_sq = cbar;
                }
                
                // Wave speeds: s_left = ubar - â, s_right = ubar + â
                // where â is the modified sound speed from Einfeldt's method
                const double sleft = ubar - sqrt_d_sq;
                const double sright = ubar + sqrt_d_sq;
            
                // Compute intermediate wave speed s_star and pressure p_star
                // s_star represents the speed of the contact discontinuity
                const double denom = rho_L * (sleft - u_L) - rho_R * (sright - u_R);
                const double s_star = (p_R - p_L + rho_L * u_L * (sleft - u_L) 
                                     - rho_R * u_R * (sright - u_R)) / denom;
                
                // Compute intermediate pressure (constant across contact discontinuity)
                const double p_star = p_L + rho_L * (sleft - u_L) * (s_star - u_L);
                
                // Compute intermediate state densities (star region)
                const double u_stat_L = rho_L * (sleft - u_L) / (sleft - s_star);
                const double u_stat_R = rho_R * (sright - u_R) / (sright - s_star);

                // Compute intermediate state energy fluxes
                const double fe_star_L = rhoe_L * inv_rho_L + (s_star - u_L) 
                                       * (s_star + p_L / (rho_L * (sleft - u_L)));
                const double fe_star_R = rhoe_R * inv_rho_R + (s_star - u_R) 
                                       * (s_star + p_R / (rho_R * (sright - u_R)));

                // Compute internal energy per unit mass: e = p/(ρ(γ-1))
                const double e_L = inv_rho_L * p_L * inv_gamma_minus_1;
                const double e_R = inv_rho_R * p_R * inv_gamma_minus_1;
                
                // Precompute terms for energy correction
                const double inv_sdiff = 1.0 / (sright - sleft);
                const double sleft_sright = sleft * sright;
                const double energy_corr_term = sleft_sright * inv_sdiff * (e_R - e_L);

                // Initialize flux variables
                double rho_F = 0.0, rhou_F = 0.0, rhov_F = 0.0, rhoe_F = 0.0;
                
                // Determine HLLCHC numerical flux based on wave speeds
                // HLLCHC adds a heat capacity correction term to HLLC's energy flux
                if (sleft >= 0.0) {
                    // Entirely left-going waves: use left flux
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if (sleft < 0.0 && s_star >= 0.0) {
                    // Left and star region: sleft < 0 ≤ s_star
                    rho_F = rho_FL + sleft * (u_stat_L - rho_L);
                    rhou_F = rhou_FL + sleft * (u_stat_L * s_star - rhou_L);
                    rhov_F = rhov_FL + sleft * (u_stat_L * v_L - rhov_L);
                    // Add energy correction term for improved heat transfer accuracy
                    rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L - rhoe_L) + energy_corr_term * u_stat_L;
                }
                else if (s_star < 0.0 && sright > 0.0) {
                    // Star and right region: s_star < 0 < sright
                    rho_F = rho_FR + sright * (u_stat_R - rho_R);
                    rhou_F = rhou_FR + sright * (u_stat_R * s_star - rhou_R);
                    rhov_F = rhov_FR + sright * (u_stat_R * v_R - rhov_R);
                    // Add energy correction term for improved heat transfer accuracy
                    rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R - rhoe_R) + energy_corr_term * u_stat_R;
                }
                else {  // sright <= 0.0
                    // Entirely right-going waves: use right flux
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }
                
                // Store computed fluxes in output array
                flux[0][j][k] = rho_F; 
                flux[1][j][k] = rhou_F; 
                flux[2][j][k] = rhov_F;
                flux[3][j][k] = rhoe_F; 
            }
        }   
    }
    else if (dir == 2) {
        // Y-direction flux calculation - fully parallelized
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;

                 // Compute speed of sound: a = sqrt(γ * p / ρ)
                //const double a_L = sqrt(M_gamma * p_L * inv_rho_L);
                //const double a_R = sqrt(M_gamma * p_R * inv_rho_R);

                // Compute total enthalpy
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                //Calculate The Left And Right Sound Speed Values
                const double a_L = sqrt((M_gamma - 1.0) * (H_L - 0.5 * (u_L_sq + v_L_sq)));
                const double a_R = sqrt((M_gamma - 1.0) * (H_R - 0.5 * (u_R_sq + v_R_sq)));

                // Compute Euler fluxes in y-direction:
                // G(ρ) = ρv, G(ρu) = ρuv, G(ρv) = ρv² + p, G(ρe) = ρHv
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;
            
                // Compute Roe-averaged variables
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;

                // Compute approximate wave speeds using Roe-averaged variables,
                // following Einfeldt's approach for y-direction
                const double ubar_sq = ubar * ubar;
                const double vbar_sq = vbar * vbar;

                // Einfeldt's correction term for y-direction wave speeds
                double eta_sq = 0.5 * sqrt_rho_L * sqrt_rho_R / (sum_sqrt * sum_sqrt);
                double sqrt_d_sq = sqrt((sqrt_rho_L * a_L * a_L + sqrt_rho_R * a_R * a_R) * inv_sum_sqrt 
                                      + eta_sq * (v_R - v_L) * (v_R - v_L));
                
                // Check for valid wave speeds
                /*if (isnan(sqrt_d_sq)) {
                    // Fall back to simple Roe average
                    const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * (ubar_sq + vbar_sq)));
                    sqrt_d_sq = cbar;
                }*/
                
                // Wave speeds for y-direction: based on vbar instead of ubar
                const double sleft = vbar - sqrt_d_sq;
                const double sright = vbar + sqrt_d_sq;
            
                // Compute intermediate wave speed and pressure for y-direction
                const double denom = rho_L * (sleft - v_L) - rho_R * (sright - v_R);
                const double s_star = (p_R - p_L + rho_L * v_L * (sleft - v_L) 
                                     - rho_R * v_R * (sright - v_R)) / denom;
                
                const double p_star = p_L + rho_L * (sleft - v_L) * (s_star - v_L);
                const double v_stat_L = rho_L * (sleft - v_L) / (sleft - s_star);
                const double v_stat_R = rho_R * (sright - v_R) / (sright - s_star);

                const double fe_star_L = rhoe_L * inv_rho_L + (s_star - v_L) 
                                       * (s_star + p_L / (rho_L * (sleft - v_L)));
                const double fe_star_R = rhoe_R * inv_rho_R + (s_star - v_R) 
                                       * (s_star + p_R / (rho_R * (sright - v_R)));

                // Compute internal energy per unit mass for y-direction
                const double e_L = inv_rho_L * p_L * inv_gamma_minus_1;
                const double e_R = inv_rho_R * p_R * inv_gamma_minus_1;
                
                // Precompute terms for energy correction (y-direction)
                const double inv_sdiff = 1.0 / (sright - sleft);
                const double sleft_sright = sleft * sright;
                const double energy_corr_term = sleft_sright * inv_sdiff * (e_R - e_L);

                // Initialize flux variables for y-direction
                double rho_G = 0.0, rhou_G = 0.0, rhov_G = 0.0, rhoe_G = 0.0;
                
                // Determine HLLCHC numerical flux based on wave speeds (y-direction)
                if (sleft >= 0.0) {
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if (sleft < 0.0 && s_star >= 0.0) {
                    rho_G = rho_GL + sleft * (v_stat_L - rho_L);
                    rhou_G = rhou_GL + sleft * (v_stat_L * u_L - rhou_L);
                    rhov_G = rhov_GL + sleft * (v_stat_L * s_star - rhov_L);
                    rhoe_G = rhoe_GL + sleft * (v_stat_L * fe_star_L - rhoe_L) + energy_corr_term * v_stat_L;
                }
                else if (s_star < 0.0 && sright > 0.0) {
                    rho_G = rho_GR + sright * (v_stat_R - rho_R);
                    rhou_G = rhou_GR + sright * (v_stat_R * u_R - rhou_R);
                    rhov_G = rhov_GR + sright * (v_stat_R * s_star - rhov_R);
                    rhoe_G = rhoe_GR + sright * (v_stat_R * fe_star_R - rhoe_R) + energy_corr_term * v_stat_R;
                }
                else {  // sright <= 0.0
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }
                
                // Store computed fluxes for y-direction
                flux[0][j][k] = rho_G; 
                flux[1][j][k] = rhou_G; 
                flux[2][j][k] = rhov_G;
                flux[3][j][k] = rhoe_G; 
            }
        }
    }    
}

/**
 * RoeHC (Roe with Heat Capacity correction) Flux Calculator for Euler Equations (2D)
 * Computes numerical fluxes using Roe's approximate Riemann solver with heat capacity correction
 * for improved energy conservation in flows with significant heat transfer
 * 
 * @param dir Direction: 1 for x-direction fluxes, 2 for y-direction fluxes
 * @param var Number of variables (currently fixed at 4)
 * @param rows Number of grid points in y-direction
 * @param cols Number of grid points in x-direction
 * @param GC Number of ghost cells
 * @param L_state Left state conservative variables at cell interfaces
 * @param R_state Right state conservative variables at cell interfaces
 * @param flux Output flux arrays
 */
static inline void RoeHC_Flux(int dir, int var, int rows, int cols, int GC, 
                              double (*L_state)[rows][cols], 
                              double (*R_state)[rows][cols], 
                              double (*flux)[rows][cols]) {
    
    // Parameters for entropy fix (Harten's entropy correction)
    const double epsilon = 1e-6;
    const double epsilon2 = epsilon * epsilon;
    const double two_epsilon = 2.0 * epsilon;
    
    // Small value to prevent division by zero
    const double eps = 1e-12;
    
    // Precompute frequently used constants
    const double gamma_minus_1 = M_gamma - 1.0;
    const double gamma_ratio = M_gamma / gamma_minus_1;      // γ/(γ-1) for enthalpy
    const double inv_gamma_minus_1 = 1.0 / gamma_minus_1;    // 1/(γ-1) for internal energy
    
    if (dir == 1) { // x-direction flux calculation
        #pragma omp parallel for collapse(2)
        for (int j = GC-1; j < rows-GC; j++) {
            for (int k = GC; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables using division optimization
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;      // x-velocity left state
                const double u_R = rhou_R * inv_rho_R;      // x-velocity right state
                const double v_L = rhov_L * inv_rho_L;      // y-velocity left state
                const double v_R = rhov_R * inv_rho_R;      // y-velocity right state
                
                // Compute velocity squares for kinetic energy calculation
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                // Compute pressure using ideal gas law: p = (ρe - 0.5ρ(u²+v²)) * (γ-1)
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;

                // Compute internal energy per unit mass: e = p/(ρ(γ-1))
                // Used for heat capacity correction term
                const double e_L = inv_rho_L * p_L * inv_gamma_minus_1;
                const double e_R = inv_rho_R * p_R * inv_gamma_minus_1;

                // Compute total enthalpy: H = 0.5(u²+v²) + γ/(γ-1) * p/ρ
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in x-direction:
                // F(ρ) = ρu, F(ρu) = ρu² + p, F(ρv) = ρuv, F(ρe) = ρHu
                const double rho_FL = rho_L * u_L;
                const double rho_FR = rho_R * u_R;
                const double rhou_FL = rho_L * u_L_sq + p_L;
                const double rhou_FR = rho_R * u_R_sq + p_R;
                const double rhov_FL = rho_L * v_L * u_L;
                const double rhov_FR = rho_R * v_R * u_R;
                const double rhoe_FL = rho_L * H_L * u_L;
                const double rhoe_FR = rho_R * H_R * u_R;

                // Compute Roe-averaged variables for intermediate state estimation
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double rhobar = sqrt_rho_L * sqrt_rho_R;          // Roe-averaged density
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;  // Roe-averaged x-velocity
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;  // Roe-averaged y-velocity
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;  // Roe-averaged enthalpy
                const double q2bar = ubar * ubar + vbar * vbar;         // Roe-averaged velocity magnitude squared
                
                // Compute Roe-averaged speed of sound: cbar = √[(γ-1)(Hbar - 0.5*q2bar)]
                const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * q2bar));

                // Compute eigenvalues using Roe-averaged variables
                // For x-direction: λ₁ = ubar-cbar, λ₂ = ubar, λ₃ = ubar, λ₄ = ubar+cbar
                double eigenvalues[4];
                eigenvalues[0] = ubar - cbar;  // Left acoustic wave
                eigenvalues[1] = ubar;         // Entropy wave
                eigenvalues[2] = ubar;         // Shear wave
                eigenvalues[3] = ubar + cbar;  // Right acoustic wave

                // Apply entropy fix (Harten's entropy correction)
                // Prevents expansion shocks and ensures entropy condition is satisfied
                for (int i = 0; i < 4; i++) {
                    const double abs_lambda = fabs(eigenvalues[i]);
                    if (abs_lambda < epsilon) {
                        // Smooth transition near zero eigenvalues
                        eigenvalues[i] = (eigenvalues[i] * eigenvalues[i] + epsilon2) / two_epsilon;
                    } else {
                        eigenvalues[i] = abs_lambda;
                    }
                }

                // Compute wave strengths (α coefficients) for characteristic decomposition
                const double drho = rho_R - rho_L;      // Density jump
                const double du = u_R - u_L;            // x-velocity jump
                const double dv = v_R - v_L;            // y-velocity jump
                const double dp = p_R - p_L;            // Pressure jump
                
                const double cbar_sq = cbar * cbar;
                const double inv_2cbar_sq = 1.0 / (2.0 * cbar_sq + eps);
                const double inv_cbar_sq = 1.0 / (cbar_sq + eps);

                // Wave strengths (α) from the jump in characteristic variables
                const double alpha1 = (dp - rhobar * cbar * du) * inv_2cbar_sq;  // Left acoustic wave
                const double alpha2 = drho - dp * inv_cbar_sq;                   // Entropy wave
                const double alpha3 = rhobar * dv;                               // Shear wave
                const double alpha4 = (dp + rhobar * cbar * du) * inv_2cbar_sq;  // Right acoustic wave

                // Right eigenvectors of the Roe matrix for x-direction
                const double eigenvectors[4][4] = {
                    {1.0, ubar - cbar, vbar, Hbar - ubar * cbar},  // Left acoustic wave
                    {1.0, ubar, vbar, 0.5 * q2bar},                // Entropy wave
                    {0.0, 0.0, 1.0, vbar},                         // Shear wave
                    {1.0, ubar + cbar, vbar, Hbar + ubar * cbar}   // Right acoustic wave
                };

                // Compute Roe numerical flux: F_Roe = 0.5*(F_L + F_R) - 0.5*Σ|λ|αR
                // Central part (average of left and right fluxes)
                double rho_F = 0.5 * (rho_FL + rho_FR);
                double rhou_F = 0.5 * (rhou_FL + rhou_FR);
                double rhov_F = 0.5 * (rhov_FL + rhov_FR);
                double rhoe_F = 0.5 * (rhoe_FL + rhoe_FR);

                // Add dissipation term (characteristic-based upwinding)
                const double wave_strengths[4] = {alpha1, alpha2, alpha3, alpha4};
                
                // Loop over all characteristic waves
                for (int i = 0; i < 4; i++) {
                    const double term = 0.5 * eigenvalues[i] * wave_strengths[i];
                    rho_F -= term * eigenvectors[i][0];
                    rhou_F -= term * eigenvectors[i][1];
                    rhov_F -= term * eigenvectors[i][2];
                    rhoe_F -= term * eigenvectors[i][3];
                }

                // Add heat capacity correction term to energy flux
                // This term improves energy conservation for flows with significant heat transfer
                const double energy_correction = 0.5 * eigenvalues[3] * (e_R - e_L) * rhobar;
                
                // Store computed fluxes with energy correction
                flux[0][j][k] = rho_F;
                flux[1][j][k] = rhou_F;
                flux[2][j][k] = rhov_F;
                flux[3][j][k] = rhoe_F - energy_correction;
            }
        }
    }
    else if (dir == 2) { // y-direction flux calculation
        #pragma omp parallel for collapse(2)
        for (int j = GC; j < rows-GC; j++) {
            for (int k = GC-1; k < cols-GC; k++) {
                // Read left and right conservative variables
                const double rho_L = L_state[0][j][k];
                const double rho_R = R_state[0][j][k];
                const double rhou_L = L_state[1][j][k];
                const double rhou_R = R_state[1][j][k];
                const double rhov_L = L_state[2][j][k];
                const double rhov_R = R_state[2][j][k];
                const double rhoe_L = L_state[3][j][k];
                const double rhoe_R = R_state[3][j][k];

                // Compute primitive variables
                const double inv_rho_L = 1.0 / rho_L;
                const double inv_rho_R = 1.0 / rho_R;
                
                const double u_L = rhou_L * inv_rho_L;
                const double u_R = rhou_R * inv_rho_R;
                const double v_L = rhov_L * inv_rho_L;
                const double v_R = rhov_R * inv_rho_R;
                
                const double u_L_sq = u_L * u_L;
                const double u_R_sq = u_R * u_R;
                const double v_L_sq = v_L * v_L;
                const double v_R_sq = v_R * v_R;
                
                const double p_L = (rhoe_L - 0.5 * rho_L * (u_L_sq + v_L_sq)) * gamma_minus_1;
                const double p_R = (rhoe_R - 0.5 * rho_R * (u_R_sq + v_R_sq)) * gamma_minus_1;
                
                // Compute internal energy per unit mass
                const double e_L = inv_rho_L * p_L * inv_gamma_minus_1;
                const double e_R = inv_rho_R * p_R * inv_gamma_minus_1;

                // Compute total enthalpy
                const double H_L = 0.5 * (u_L_sq + v_L_sq) + gamma_ratio * (p_L * inv_rho_L);
                const double H_R = 0.5 * (u_R_sq + v_R_sq) + gamma_ratio * (p_R * inv_rho_R);

                // Compute Euler fluxes in y-direction:
                // G(ρ) = ρv, G(ρu) = ρuv, G(ρv) = ρv² + p, G(ρe) = ρHv
                const double rho_GL = rho_L * v_L;
                const double rho_GR = rho_R * v_R;
                const double rhou_GL = rho_L * u_L * v_L;
                const double rhou_GR = rho_R * u_R * v_R;
                const double rhov_GL = rho_L * v_L_sq + p_L;
                const double rhov_GR = rho_R * v_R_sq + p_R;
                const double rhoe_GL = rho_L * H_L * v_L;
                const double rhoe_GR = rho_R * H_R * v_R;

                // Compute Roe-averaged variables for y-direction
                const double sqrt_rho_L = sqrt(rho_L);
                const double sqrt_rho_R = sqrt(rho_R);
                const double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                const double inv_sum_sqrt = 1.0 / sum_sqrt;
                
                const double rhobar = sqrt_rho_L * sqrt_rho_R;
                const double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) * inv_sum_sqrt;
                const double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) * inv_sum_sqrt;
                const double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) * inv_sum_sqrt;
                const double q2bar = ubar * ubar + vbar * vbar;
                
                // Compute Roe-averaged speed of sound for y-direction
                const double cbar = sqrt(gamma_minus_1 * (Hbar - 0.5 * q2bar));

                // Compute eigenvalues for y-direction
                // For y-direction: λ₁ = vbar-cbar, λ₂ = vbar, λ₃ = vbar, λ₄ = vbar+cbar
                double eigenvalues[4];
                eigenvalues[0] = vbar - cbar;
                eigenvalues[1] = vbar;
                eigenvalues[2] = vbar;
                eigenvalues[3] = vbar + cbar;

                // Apply entropy fix
                for (int i = 0; i < 4; i++) {
                    const double abs_lambda = fabs(eigenvalues[i]);
                    if (abs_lambda < epsilon) {
                        eigenvalues[i] = (eigenvalues[i] * eigenvalues[i] + epsilon2) / two_epsilon;
                    } else {
                        eigenvalues[i] = abs_lambda;
                    }
                }

                // Compute wave strengths for y-direction
                const double drho = rho_R - rho_L;
                const double du = u_R - u_L;
                const double dv = v_R - v_L;
                const double dp = p_R - p_L;
                
                const double cbar_sq = cbar * cbar;
                const double inv_2cbar_sq = 1.0 / (2.0 * cbar_sq + eps);
                const double inv_cbar_sq = 1.0 / (cbar_sq + eps);

                // Wave strengths for y-direction (note: du and dv swapped compared to x-direction)
                const double alpha1 = (dp - rhobar * cbar * dv) * inv_2cbar_sq;  // Left acoustic wave
                const double alpha2 = drho - dp * inv_cbar_sq;                   // Entropy wave
                const double alpha3 = rhobar * du;                               // Shear wave
                const double alpha4 = (dp + rhobar * cbar * dv) * inv_2cbar_sq;  // Right acoustic wave

                // Right eigenvectors of the Roe matrix for y-direction
                const double eigenvectors[4][4] = {
                    {1.0, ubar, vbar - cbar, Hbar - vbar * cbar},  // Left acoustic wave
                    {1.0, ubar, vbar, 0.5 * q2bar},                // Entropy wave
                    {0.0, 1.0, 0.0, ubar},                         // Shear wave
                    {1.0, ubar, vbar + cbar, Hbar + vbar * cbar}   // Right acoustic wave
                };

                // Compute Roe numerical flux for y-direction
                double rho_G = 0.5 * (rho_GL + rho_GR);
                double rhou_G = 0.5 * (rhou_GL + rhou_GR);
                double rhov_G = 0.5 * (rhov_GL + rhov_GR);
                double rhoe_G = 0.5 * (rhoe_GL + rhoe_GR);

                // Add dissipation term
                const double wave_strengths[4] = {alpha1, alpha2, alpha3, alpha4};
                
                for (int i = 0; i < 4; i++) {
                    const double term = 0.5 * eigenvalues[i] * wave_strengths[i];
                    rho_G -= term * eigenvectors[i][0];
                    rhou_G -= term * eigenvectors[i][1];
                    rhov_G -= term * eigenvectors[i][2];
                    rhoe_G -= term * eigenvectors[i][3];
                }

                // Add heat capacity correction term to energy flux (y-direction)
                const double energy_correction = 0.5 * eigenvalues[3] * (e_R - e_L) * rhobar;
                
                // Store computed fluxes with energy correction
                flux[0][j][k] = rho_G;
                flux[1][j][k] = rhou_G;
                flux[2][j][k] = rhov_G;
                flux[3][j][k] = rhoe_G - energy_correction;
            }
        }
    }
}


static inline void HLL_Flux_HeatConduction(int var, int rows, int GC, double (*x)[rows], double (*y)[rows] ,double (*z)[rows], double u_point) {
    
    for (int j = GC-1; j <= rows-GC; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        double a_L = sqrt(M_gamma * p_L / rho_L);
        double a_R = sqrt(M_gamma * p_R / rho_R);

        //计算需要使用的参数
        //计算声速

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

        //计算温度变量
        double rhobar = sqrt(rho_L*rho_R);
//        double rhobar = 0.5*(rho_L+rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;
        double e_L = T_L/(M_gamma-1);
        double e_R = T_R/(M_gamma-1);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(M_gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(M_gamma-1);
        double rhoe_LL = rhoe_L/rho_L;
        double rhoe_RR = rhoe_R/rho_R;
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double Splus = max_of_two(fabs(u_L)  + a_L, fabs(u_R) + a_R);

        double s_HLL = 0.5 * (sleft + sright);
        

        //确定HLL数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        if (sleft > 0){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft <= 0 && sright >=0){
            double rho_HLL = (sright*rho_R - sleft*rho_L + rho_FL - rho_FR)/(sright-sleft);
            rho_F = (sright*rho_FL - sleft*rho_FR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
            rhou_F = (sright*rhou_FL - sleft*rhou_FR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
            rhoe_F = (sright*rhoe_FL - sleft*rhoe_FR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft) +  (sleft*sright)/(sright-sleft) * (e_R-e_L) * rho_HLL;
            
        }
        else if (sright < 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLLC_Flux_HeatConduction(int var, int rows, int GC, double (*x)[rows], double (*y)[rows] ,double (*z)[rows], double u_point) {
    
    for (int j = GC-1; j <= rows-GC; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        double a_L = sqrt(M_gamma * p_L / rho_L);
        double a_R = sqrt(M_gamma * p_R / rho_R);

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

         //计算温度变量
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;

        double e_L = T_L / (M_gamma-1);
        double e_R = T_R / (M_gamma-1);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(M_gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(M_gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double s_star = (p_R - p_L + rho_L*u_L*(sleft - u_L) - rho_R*u_R*(sright - u_R))\
                                 /(rho_L*(sleft - u_L) - rho_R*(sright - u_R));
        double p_star = p_L + rho_L * (u_L - sleft) * (u_L -s_star);

        //HLLC通量需要的中间变量
        double rho_star_L = rho_L * (sleft-u_L)/(sleft-s_star);
        double rho_star_R = rho_R * (sright-u_R)/(sright-s_star);
        double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L)));
        double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R)));

    //    double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L))) + (p_R/rho_R - p_L/rho_L)/(M_gamma-1);
    //    double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R))) - (p_R/rho_R - p_L/rho_L)/(M_gamma-1);

        double T_star_L = p_star / rho_star_L;
        double T_star_R = p_star / rho_star_R;

        //热扩散相对的速度计算
        //double Splus = max_of_two(fabs(u_L) + a_L, fabs(u_R) + a_R);
        double s_L_P = sleft - u_point;
        double s_R_P = sright - u_point;

        //定义运动的波的相对波速：
        double s_Contact_L = u_L - u_point;
        double s_LW_L = u_L - a_L - u_point;
        double s_RW_L = u_L + a_L - u_point;
        double s_Contact_R = u_R- u_point;
        double s_LW_R = u_R - a_R- u_point;
        double s_RW_R = u_R + a_R- u_point;
        double ma = u_L/cbar;
        //Splus = fabs(ubar) + cbar;

        double rho_HLL = (sright*rho_R - sleft*rho_L + rho_FL - rho_FR)/(sright-sleft);
        double rhou_HLL = (sright*rhou_R - sleft*rhou_L + rhou_FL - rhou_FR)/(sright-sleft);
        double rhoe_HLL = (sright*rhoe_R - sleft*rhoe_L + rhoe_FL - rhoe_FR)/(sright-sleft);


        //HLLC数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        
        if (sleft >= 0){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && s_star >=0 ){
            rho_F = rho_FL + sleft * (rho_star_L - rho_L);
            rhou_F = rhou_FL + sleft *(rho_star_L * s_star  - rhou_L);
//            rhoe_F = rhoe_FL + sleft * (rho_star_L * fe_star_L -  rhoe_L) - (1.0/4.0) * rho_HLL * (sright-sleft) * (e_R-e_L);
            rhoe_F = rhoe_FL + sleft * (rho_star_L * fe_star_L -  rhoe_L) + (sleft*sright)/(sright-sleft) * (e_R-e_L) * rho_HLL;
        }
        else if(s_star < 0 && sright > 0){
            rho_F = rho_FR + sright * (rho_star_R - rho_R);
            rhou_F = rhou_FR + sright *(rho_star_R * s_star  - rhou_R);
//            rhoe_F = rhoe_FR + sright * (rho_star_R * fe_star_R -  rhoe_R) - (1.0/4.0) * rho_HLL * (sright-sleft) * (e_R-e_L);
//            rhoe_F = 0.5*(rhoe_FR + sright * (rho_star_R * fe_star_R -  rhoe_R) + rhoe_FL + sleft * (rho_star_L * fe_star_L -  rhoe_L));
            rhoe_F = rhoe_FR + sright * (rho_star_R * fe_star_R -  rhoe_R) + (sleft*sright)/(sright-sleft) * (e_R-e_L) * rho_HLL;
        }
        else if (sright  <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }


        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 

    }   
}



static inline void Roe_Flux_HeatConduction(int var, int rows, int GC, double (*x)[rows], double (*y)[rows] ,double (*z)[rows], double u_point) {
    double epsilon = 1e-6;
    
    for (int j = GC-1; j <= rows-GC; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        double a_L = sqrt(M_gamma * p_L / rho_L);
        double a_R = sqrt(M_gamma * p_R / rho_R);

        //计算需要使用的参数

        //计算热力学变量
        double H_L = 0.5 * pow(u_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;
        double e_L = T_L / (M_gamma-1);
        double e_R = T_R / (M_gamma-1);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均物理量
        double rhobar = sqrt(rho_L * rho_R); 
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));
        double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2)));

        double sleft = ubar - cbar;
        double sright = ubar + cbar;

         //热扩散相对的速度计算
        double Splus = max_of_two(fabs(u_L) + a_L, fabs(u_R) + a_R);
        double s_L_P = sleft - u_point;
        double s_R_P = sright - u_point;

        //定义运动的波的相对波速：
        double s_Contact_L = u_L - u_point;
        double s_LW_L = u_L - a_L - u_point;
        double s_RW_L = u_L + a_L - u_point;
        double s_Contact_R = u_R- u_point;
        double s_LW_R = u_R - a_R- u_point;
        double s_RW_R = u_R + a_R- u_point;

        //利用Roe平均计算稳定项
        double lambda1, lambda2, lambda3;
        double alpha1, alpha2, alpha3;
        lambda1 = fabs(ubar - cbar);
        if (lambda1 < epsilon){
            lambda1 = (fabs(ubar - cbar) + pow(epsilon,2)) / (2 * epsilon);
        }
        lambda2 = fabs(ubar);
        if (lambda2 < epsilon){
            lambda2 = (fabs(ubar) + pow(epsilon,2)) / (2 * epsilon);
        }
        lambda3 = fabs(ubar + cbar);
        if (lambda3 < epsilon){
            lambda3 = (fabs(ubar + cbar) + pow(epsilon,2)) / (2 * epsilon);
        }

        alpha1 = ((p_R - p_L) - rhobar * cbar * (u_R - u_L)) / (2 * pow(cbar,2));
        alpha2 = (rho_R - rho_L) - (p_R - p_L) / pow(cbar,2);
        alpha3 = ((p_R - p_L) + rhobar * cbar * (u_R - u_L)) / (2 * pow(cbar,2));


        //Roe数值通量
        double rho_HLL = (sright*rho_R - sleft*rho_L + rho_FL - rho_FR)/(sright-sleft);
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        rho_F = 0.5*(rho_FL + rho_FR) - 0.5 *(lambda1*alpha1 + lambda2*alpha2 + lambda3*alpha3);
        rhou_F = 0.5 * (rhou_FL + rhou_FR)\
                    - 0.5*(lambda1 * alpha1 * (ubar- cbar)\
                    + lambda2 * alpha2 * ubar\
                    + lambda3 * alpha3 * (ubar + cbar));
        rhoe_F = 0.5*(rhoe_FL + rhoe_FR)\
                    -0.5*(lambda1 * alpha1 * (Hbar - ubar*cbar)\
                    +lambda2 * alpha2 * 0.5 * pow(ubar,2)\
                    +lambda3 * alpha3 * (Hbar + ubar*cbar));

        if (s_L_P * s_R_P < 0)
            rhoe_F = rhoe_F + (sleft*sright)/(sright-sleft) * (e_R-e_L) * rho_HLL;

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}


static inline void ER_Flux_PlusHeat(int var, int rows, int GC, double (*x)[rows], double (*y)[rows] ,double (*z)[rows]) {
    
    for (int j = GC-1; j <= rows-GC; j++) {
        double rho_F = 0,u_F=0,p_F=0;
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        double c_L = sqrt(M_gamma*p_L/rho_L);
        double c_R = sqrt(M_gamma*p_R/rho_R);

        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;
        double Splus = max_of_two(fabs(u_L)  + c_L, fabs(u_R) + c_R);
        double rhoabr = sqrt(rho_L * rho_R);

        double H_L = 0.5 * pow(u_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(M_gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(M_gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;

        double p_star = ExactRieamnna_pstar(rho_L, rho_R, u_L, u_R, p_L, p_R, M_gamma, M_gamma ,1e-4, 1e3);
		double u_star = ExactRieamnna_ustar(p_star, rho_L, rho_R, u_L,u_R, p_L,p_R, M_gamma,M_gamma);
        
         
        //分别考虑左激波，左稀疏波的情况
        if (p_star >=  p_L){
            double part1 = (M_gamma + 1)* p_star + (M_gamma - 1) * p_L ;
            double part2 = (M_gamma - 1) * p_star + (M_gamma + 1) * p_L;
            double rho_star_L = rho_L * part1 / part2;
            double A_L = 2.0 / ((M_gamma + 1) * rho_L);
            //计算激波的速度
            double B_L = (M_gamma - 1) * p_L  / (M_gamma + 1);
            double S_L = u_L - (1 / rho_L) * sqrt((B_L + p_star) / A_L); 
            double si = 0;
            if (si <= S_L){
                rho_F = rho_L;
                u_F = u_L;
                p_F = p_L;
            }
            else if (S_L < si && si <= u_star)
            {
                rho_F = rho_star_L;
                u_F = u_star;
                p_F = p_star;
            }
        }
        else {
            //计算接触间断左侧密度
            double part1 = p_star;
            double part2 = p_L;
            double rho_star_L = rho_L * pow(part1/part2, 1/M_gamma);
            //计算左稀疏波波头和波尾的速度
            double aL_star = c_L * pow((p_star / p_L), ((M_gamma - 1) / (2 *M_gamma)));
            double S_HL = u_L - c_L;
            double S_TL = u_star - aL_star;
            double si=0;
            if (si <= S_HL){
                rho_F = rho_L;
                u_F = u_L;
                p_F = p_L;
            }
            else if (si <= S_TL && si > S_HL)
            {
                rho_F = rho_L * pow(2/(M_gamma + 1) + (M_gamma - 1) * (u_L - si) / ((M_gamma + 1) * c_L), (2 / (M_gamma - 1))); 
                u_F = 2 * (c_L + (M_gamma - 1) * u_L/2 + si) / (M_gamma + 1);
                p_F = p_L * pow(2 / (M_gamma + 1) + (M_gamma - 1) * (u_L - si) / ((M_gamma + 1) * c_L),  2 * M_gamma / (M_gamma - 1));
            }
            else if (si <= u_star && si > S_TL)
            {
                rho_F = rho_star_L;
                u_F = u_star;
                p_F = p_star;
            }
        }

        //分别考虑右稀疏波和右激波的情况
        if ( p_star >=  p_R){
            double part1 =  (M_gamma + 1) * p_star + (M_gamma - 1) * p_R;
            double part2 =  (M_gamma - 1) * p_star + (M_gamma + 1) * p_R;
            double rho_star_R = rho_R * part1 / part2;
            double A_R = 2 / ((M_gamma + 1) * rho_R);
            //计算激波的速度
            double B_R = (M_gamma - 1) * p_R / (M_gamma + 1);
            double S_R = u_R + (1 / rho_R) * sqrt((B_R + p_star) / A_R); 
            double si = 0;
            if (si <= S_R &&  u_star < si)
            {
                rho_F = rho_star_R;
                u_F = u_star;
                p_F = p_star;
            }
            else if (S_R < 0){
                rho_F = rho_R;
                u_F = u_R;
                p_F = p_R;
            }
        }
        else{
            //计算接触间断左侧密度
            double part1 = p_star;
            double part2 = p_R;
            double rho_star_R = rho_R * pow(part1 / part2, 1/M_gamma);
            //计算左稀疏波波头和波尾的速度
            double aR_star = c_R * pow((p_star / p_R), ((M_gamma - 1) / (2 *M_gamma)));
            double S_HR = u_R + c_R;
            double S_TR = u_star + aR_star;   
            double si=0;
            if (si <= S_TR && si > u_star)
            {
                rho_F = rho_star_R;
                u_F = u_star;
                p_F = p_star;
            }    
            else if (si<= S_HR && si >S_TR)
            {
                rho_F = rho_R * pow(2/(M_gamma + 1) - (M_gamma - 1) * (u_R - si) / ((M_gamma + 1) * c_R), (2 / (M_gamma - 1))); 
                u_F = 2 * (-c_R + (M_gamma - 1) * u_R / 2 + si) / (M_gamma + 1);
                p_F = p_R * pow(2 / (M_gamma + 1) - (M_gamma - 1) * (u_R - si) / ((M_gamma + 1) * c_R),  2 * M_gamma / (M_gamma - 1));
            }
            else if (si > S_HR){
                rho_F = rho_R;
                u_F = u_R;
                p_F = p_R;
            }
        }

        z[0][j] = rho_F * u_F; 
        z[1][j] = rho_F * u_F *u_F + p_F ; 
        z[2][j] = ((0.5*rho_F * u_F *u_F + p_F/(M_gamma-1)) + p_F) * u_F -  Splus * (rhoe_R - rhoe_L) * rhoabr; 
    }   
}

/*……………………………………………………………………………………………………*/
//“The algorithmic description of Marquina’s flux formula is as follows:” ([Donat 和 Marquina, 1996, p. 44]
static inline void RS_Marquina(int var, int rows,double (*y)[rows],double (*z)[rows], double dt,double dx) {
    double conserl[3][rows],conserr[3][rows];
    double fluxl[3][rows],fluxr[3][rows];
    double prileft[3][rows], priright[3][rows];
    double pri[var][rows];
    double slope[3][rows],a[3][rows],b[3][rows];
    
    double eigen_l[3][3][rows],eigen_r[3][3][rows];
    double w_l[var][rows], w_r[var][rows];
    double phi_fl[var][rows], phi_fr[var][rows];
    double phi_fp[var][rows], phi_fm[var][rows];
    double flux[3][rows];
    double lamda[3][rows];
    double alpha[var][rows];

    //初始化数组
     for (int k = 0; k < var; k++){
        for (int j = 0; j < rows; j++){
            w_l[k][j] = 0.0;
            w_r[k][j] = 0.0;
            phi_fl[k][j] = 0.0;
            phi_fr[k][j] = 0.0;
            phi_fp[k][j] = 0.0;
            phi_fm[k][j] = 0.0;
            flux[k][j] = 0.0;
        }
        
    }

//    Con_to_Pri_1D(3,rows,pri,y);
    //计算特征矩阵
	for(int j=0;  j < rows; j++) {
			
        double  q2, c2, b1, b2;
        double _u, _H, _c;
        //preparing some interval value
        _u = pri[1][j];
        _H = 0.5 * pow(pri[1][j],2) + M_gamma * pri[2][j] /((M_gamma-1) * pri[0][j]);
        q2 = _u*_u ;
        c2 = (M_gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
        _c = sqrt(c2);
        b1 = (M_gamma - 1.0)/c2;
        b2 = 1.0 + b1*q2 - b1*_H;
        // left eigen vectors 
        eigen_l[0][0][j] = 0.5*(b2 + _u/_c);
        eigen_l[0][1][j] = -0.5*(b1*_u + 1/_c);
        eigen_l[0][2][j] = 0.5*b1;
            
        eigen_l[1][0][j] = -q2 + _H;
        eigen_l[1][1][j] = _u;;
        eigen_l[1][2][j] = -1.0;

        eigen_l[2][0][j] = 0.5*(b2 - _u/_c);
        eigen_l[2][1][j] = 0.5*(-b1*_u + 1/_c);
        eigen_l[2][2][j] = 0.5*b1;

        //right eigen vectors
        eigen_r[0][0][j] = 1.0;
        eigen_r[0][1][j] = b1;
        eigen_r[0][2][j] = 1.0;
            
        eigen_r[1][0][j] = _u - _c;
        eigen_r[1][1][j] = _u*b1;
        eigen_r[1][2][j] = _u + _c;

        eigen_r[2][0][j] = _H - _u*_c;
        eigen_r[2][1][j] = _H*b1 - 1.0;
        eigen_r[2][2][j] = _H + _u*_c;


        //计算对应特征值
        lamda[0][j] = _u - _c ;
        lamda[1][j] = _u;
        lamda[2][j] = _u + _c;
    }


    switch (Recon_Accur){
        case 0:
//            Reconstruction_Godunov(var,rows,3,y,conserl,conserr,dx);
            break;
        case 1:
//            TVD_Reconstruction(var,rows,3,y,conserl,conserr,dx);
            break;
        case 3:
//            WENO3_Reconstruction(var,rows,3,y,conserl,conserr);
            break;
        case 5:
//            WENO5_Reconstruction_C(var,rows,3,y,conserl,conserr);
            break;
        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }
                              

    
    //执行计算Marquina flux
//    Con_to_Pri_1D(3,rows,prileft,conserl);
//    Con_to_Pri_1D(3,rows,priright,conserr);

//    initEulerflux1D(var, rows, conserl, fluxl);
//    initEulerflux1D(var, rows, conserr, fluxr);

    // 投影到特征空间
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < rows-1; j++){
            for (int m = 0; m < 3; m++){
                w_l[k][j] += conserl[m][j]*eigen_l[k][m][j];
                phi_fl[k][j] += fluxl[m][j]*eigen_l[k][m][j];
            }
        }
        
    }
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < rows-1; j++){
            for (int m = 0; m < 3; m++){
                w_r[k][j] += conserr[m][j]*eigen_l[k][m][j+1];
                phi_fr[k][j] += fluxr[m][j]*eigen_l[k][m][j+1];
            }
        }
        
    }

    //计算通量分量
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < rows-1; j++){
            if (lamda[k][j]  * lamda[k][j+1] > 0){
                if (lamda[k][j] > 0){
                    phi_fp[k][j] = phi_fl[k][j];
                    phi_fm[k][j] = 0.0;
                }
                else{
                    phi_fp[k][j] = 0.0;
                    phi_fm[k][j] = phi_fr[k][j];
                }
            }
            else{
                alpha[k][j] = max_of_two(fabs(lamda[k][j]),fabs(lamda[k][j+1]));
                phi_fp[k][j] = 0.5*(phi_fl[k][j] + alpha[k][j] * w_l[k][j]);
                phi_fm[k][j] = 0.5*(phi_fr[k][j] - alpha[k][j] * w_r[k][j]);
            }
        }
    }

    //投影到物理空间计算Flux
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < rows-1; j++){
            for (int m = 0; m < 3; m++){
               flux[k][j] += phi_fp[m][j]*eigen_r[k][m][j] + phi_fm[m][j] * eigen_r[k][m][j+1];
            }
        }  
    }

    for (int i = 0; i < var; i++){
        for (int j = 1; j < rows-1; j++){
            z[i][j] = flux[i][j];
        }
    }

}

static inline void Source_Gravity(int var, int rows, int cols, int GC, 
                                 double (*x)[rows][cols], double (*source)[rows][cols]) {
    
    // Initialize source term to zero
    #pragma omp parallel for collapse(2)
    for (int i = GC; i <= rows-GC-1; i++) {
        for (int j = GC; j <= cols-GC-1; j++) {
            for (int k = 0; k < var; k++) {
                source[k][i][j] = 0.0;
            }
        }
    }
    
    // Set gravity source term for Euler equations
    #pragma omp parallel for collapse(2)
    for (int i = GC; i <= rows-GC-1; i++) {
        for (int j = GC; j <= cols-GC-1; j++) {
            double rho = x[0][i][j];          // Density
            double rhov = x[2][i][j];         // y-momentum
            double v = rhov / rho;            // y-velocity (avoid division by zero)
            
            // Gravity source terms for Euler equations
            source[0][i][j] = 0.0;           // Mass equation: no source term
            source[1][i][j] = 0.0;           // x-momentum equation: no source term  
            source[2][i][j] = rho * Gravity; // y-momentum equation: gravity force
            source[3][i][j] = rho * Gravity * v; // Energy equation: work done by gravity
        }
    }
}

/*
 * Compute spatial discretization term for Euler equations
 * L = -∂f/∂x - ∂g/∂y + source
 * Using finite volume method with first-order upwind scheme
 */
static inline void Space_Discrete_Item(int var, int rows, int cols, int GC, double delta_x, double delta_y,
                                        double (*f)[rows][cols], double (*g)[rows][cols], double (*s)[rows][cols],
                                        double (*L)[rows][cols]) {
    
    const double inv_dx = 1.0 / delta_x;  // Inverse of x-spacing for efficiency
    const double inv_dy = 1.0 / delta_y;  // Inverse of y-spacing for efficiency
    
    // Compute spatial derivative for interior cells only
    #pragma omp parallel for collapse(3)
    for (int i = GC; i <= rows - GC - 1; i++) {
        for (int j = GC; j <= cols - GC - 1; j++) {
            for (int k = 0; k < var; k++) {
                // Finite volume discretization: ∂f/∂x ≈ (f_{i+1/2} - f_{i-1/2})/Δx
                L[k][i][j] = -(f[k][i][j] - f[k][i-1][j]) * inv_dx 
                             -(g[k][i][j] - g[k][i][j-1]) * inv_dy 
                             + s[k][i][j];
            }
        }
    }
}

#endif  