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
static inline double WENO3_L();
static inline double WENO3_R();
static inline double WENO5_L();
static inline double WENO5_R();


/*                               ************************************                               */
/*                               ************************************                               */
/*                                    Reconstruction Scheme                                         */
/*                               ************************************                               */
/*                               ************************************                               */


/*                                      ******************                                          */
/*                                            Godunov                                               */
/*                                      ******************                                          */
static inline void Reconstruction_Godunov(int dir, int var, int rows, int cols, int GC, 
                                          double (*y)[rows][cols],
                                          double (*conserl)[rows][cols], 
                                          double (*conserr)[rows][cols], 
                                          double delta_x) {
    int i, j, k;
    double epsilo = 1e-6;
    
    
    if (Characteriz) {
        // 动态分配内存
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Chara_Var)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*W_L)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*W_R)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        
        if (!Pri || !Chara_Var || !Eigen_L || !Eigen_R || !W_L || !W_R) {
            fprintf(stderr, "Memory allocation failed in Reconstruction_Godunov\n");
            free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
            return;
        }

        // 转换到原始变量
        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            Compute_Eigen_2D(1.0, 0.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // 并行计算特征变量
            #pragma omp parallel for
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            // 并行赋值 W_L 和 W_R
            #pragma omp parallel for
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        W_L[i][j][k] = Chara_Var[i][j][k];
                        W_R[i][j][k] = Chara_Var[i][j+1][k];
                    }
                }
            }

            // 并行计算守恒量
            #pragma omp parallel for
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        else if (dir == 2) {
            Compute_Eigen_2D(0.0, 1.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // 并行计算特征变量
            #pragma omp parallel for
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            // 并行赋值 W_L 和 W_R
            #pragma omp parallel for
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        W_L[i][j][k] = Chara_Var[i][j][k];
                        W_R[i][j][k] = Chara_Var[i][j][k+1];
                    }
                }
            }

            // 并行计算守恒量
            #pragma omp parallel for
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        
        free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
    } 
    else {
        if (dir == 1) {
            #pragma omp parallel for
            for (int j = GC - 1; j < rows - GC; j++) {
                for (int k = GC; k < cols - GC; k++) {
                    for (int i = 0; i < var; i++) {
                        conserl[i][j][k] = y[i][j][k];
                        conserr[i][j][k] = y[i][j + 1][k];
                    }
                }
            }
        } else if (dir == 2) {
            #pragma omp parallel for
            for (int j = GC; j < rows - GC; j++) {
                for (int k = GC - 1; k < cols - GC; k++) {
                    for (int i = 0; i < var; i++) {
                        conserl[i][j][k] = y[i][j][k];
                        conserr[i][j][k] = y[i][j][k + 1];
                    }
                }
            }
        }
    }
}
/*                                      ******************                                          */
/*                                            TVD                                                   */
/*                                      ******************                                          */
static inline void TVD_Reconstruction(int dir, int var, int rows, int cols, int GC, 
                                      double (*y)[rows][cols],
                                      double (*conserl)[rows][cols], 
                                      double (*conserr)[rows][cols], 
                                      double delta_x, double delta_y) {
    int i, j, k;
    
    if (Characteriz) {
        // 动态分配内存
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Chara_Var)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*W_L)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*W_R)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        
        if (!Pri || !Chara_Var || !Eigen_L || !Eigen_R || !W_L || !W_R) {
            fprintf(stderr, "Memory allocation failed in TVD_Reconstruction\n");
            free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
            return;
        }

        // ========== 并行初始化 ==========
        #pragma omp parallel for 
        for (i = 0; i < var; i++) {
            for (j = 0; j < rows; j++) {
                for (k = 0; k < cols; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }
        
        #pragma omp parallel for
        for (i = 0; i < var; i++) {
            for (int ii = 0; ii < var; ii++) {
                for (j = 0; j < rows; j++) {
                    for (k = 0; k < cols; k++) {
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                }
            }
        }

        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            Compute_Eigen_2D(1.0, 0.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            // ========== 计算特征变量 ==========
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            // ========== TVD重构 ==========
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = Chara_Var[i][j-2+nn][k];
                        }
                        W_L[i][j][k] = TVD_minmod_L(&fu[2], delta_x);
                        W_R[i][j][k] = TVD_minmod_R(&fu[2], delta_x);
                    }
                }
            }

            // ========== 转换回守恒量 ==========
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        else if (dir == 2) {
            Compute_Eigen_2D(0.0, 1.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j][k-2+nn];
                        }
                        W_L[i][j][k] = TVD_minmod_L(&fu[2], delta_y);
                        W_R[i][j][k] = TVD_minmod_R(&fu[2], delta_y);
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        
        free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
    } 
    else {
        if (dir == 1) {
            #pragma omp parallel for collapse(3)
            for (int j = GC - 1; j < rows - GC; j++) {
                for (int k = GC; k < cols - GC; k++) {
                    for (int i = 0; i < var; i++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j - 2 + nn][k];
                        }
                        conserl[i][j][k] = TVD_vanleer_L(&fu[2], delta_x);
                        conserr[i][j][k] = TVD_vanleer_R(&fu[2], delta_x);
                    }
                }
            }
        } else if (dir == 2) {
            #pragma omp parallel for collapse(3)
            for (int j = GC; j < rows - GC; j++) {
                for (int k = GC - 1; k < cols - GC; k++) {
                    for (int i = 0; i < var; i++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j][k - 2 + nn];
                        }
                        conserl[i][j][k] = TVD_vanleer_L(&fu[2], delta_y);
                        conserr[i][j][k] = TVD_vanleer_R(&fu[2], delta_y);
                    }
                }
            }
        }
    }
    
}

/*                                      ******************                                          */
/*                                               WENO                                               */
/*                                      ******************                                          */
// 三阶WENO重构
static inline void WENO3_Reconstruction(int dir, int var, int rows, int cols, int GC,
                                        double (*y)[rows][cols],
                                        double (*conserl)[rows][cols], 
                                        double (*conserr)[rows][cols]) {
    int i, j, k;
    
    if (Characteriz) {
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Chara_Var)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*W_L)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*W_R)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        
        if (!Pri || !Chara_Var || !Eigen_L || !Eigen_R || !W_L || !W_R) {
            fprintf(stderr, "Memory allocation failed in WENO3_Reconstruction\n");
            free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
            return;
        }

        // ========== 并行初始化 ==========
        #pragma omp parallel for collapse(3)
        for (i = 0; i < var; i++) {
            for (j = 0; j < rows; j++) {
                for (k = 0; k < cols; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }
        
        #pragma omp parallel for collapse(4)
        for (i = 0; i < var; i++) {
            for (int ii = 0; ii < var; ii++) {
                for (j = 0; j < rows; j++) {
                    for (k = 0; k < cols; k++) {
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                }
            }
        }

        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            Compute_Eigen_2D(1.0, 0.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = Chara_Var[i][j-2+nn][k];
                        }
                        W_L[i][j][k] = WENO3_L(&fu[2]);
                        W_R[i][j][k] = WENO3_R(&fu[2]);
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        else if (dir == 2) {
            Compute_Eigen_2D(0.0, 1.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = Chara_Var[i][j][k-2+nn];
                        }
                        W_L[i][j][k] = WENO3_L(&fu[2]);
                        W_R[i][j][k] = WENO3_R(&fu[2]);
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j <= rows-GC; j++) {
                    for (k = GC-1; k <= cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        
        free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
    } 
    else {
        if (dir == 1) {
            #pragma omp parallel for collapse(3)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j-2+nn][k];
                        }
                        conserl[i][j][k] = WENO3_L(&fu[2]);
                        conserr[i][j][k] = WENO3_R(&fu[2]);
                    }
                }
            }
        }
        else if (dir == 2) {
            #pragma omp parallel for collapse(3)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j][k-2+nn];
                        }
                        conserl[i][j][k] = WENO3_L(&fu[2]);
                        conserr[i][j][k] = WENO3_R(&fu[2]);
                    }
                }
            }
        }
    }
    
}



// 五阶WENO重构
static inline void WENO5_Reconstruction(int dir, int var, int rows, int cols, int GC,
                                        double (*y)[rows][cols],
                                        double (*conserl)[rows][cols], 
                                        double (*conserr)[rows][cols]) {
    int i, j, k;
    

    if (Characteriz) {
        double (*Pri)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Chara_Var)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*Eigen_L)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*Eigen_R)[var][rows][cols] = malloc(var * sizeof(double[var][rows][cols]));
        double (*W_L)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        double (*W_R)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
        
        if (!Pri || !Chara_Var || !Eigen_L || !Eigen_R || !W_L || !W_R) {
            fprintf(stderr, "Memory allocation failed in WENO5_Reconstruction\n");
            free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
            return;
        }

        // ========== 并行初始化 ==========
        #pragma omp parallel for collapse(3)
        for (i = 0; i < var; i++) {
            for (j = 0; j < rows; j++) {
                for (k = 0; k < cols; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }
        
        #pragma omp parallel for collapse(4)
        for (i = 0; i < var; i++) {
            for (int ii = 0; ii < var; ii++) {
                for (j = 0; j < rows; j++) {
                    for (k = 0; k < cols; k++) {
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                }
            }
        }

        Con_to_Pri_2D(var, rows, cols, Pri, y);
        
        if (dir == 1) {
            Compute_Eigen_2D(1.0, 0.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = Chara_Var[i][j-2+nn][k];
                        }
                        W_L[i][j][k] = WENO5_L(&fu[2]);
                        W_R[i][j][k] = WENO5_R(&fu[2]);
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        else if (dir == 2) {
            Compute_Eigen_2D(0.0, 1.0, var, rows, cols, Pri, Eigen_L, Eigen_R);
            
            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double sum = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum += y[ii][j][k] * Eigen_L[i][ii][j][k];
                        }
                        Chara_Var[i][j][k] = sum;
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = Chara_Var[i][j][k-2+nn];
                        }
                        W_L[i][j][k] = WENO5_L(&fu[2]);
                        W_R[i][j][k] = WENO5_R(&fu[2]);
                    }
                }
            }

            #pragma omp parallel for collapse(3) private(i, j, k)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double sum_l = 0.0, sum_r = 0.0;
                        for (int ii = 0; ii < var; ii++) {
                            sum_l += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            sum_r += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                        conserl[i][j][k] = sum_l;
                        conserr[i][j][k] = sum_r;
                    }
                }
            }
        }
        
        free(Pri); free(Chara_Var); free(Eigen_L); free(Eigen_R); free(W_L); free(W_R);
    } 
    else {
        if (dir == 1) {
            #pragma omp parallel for collapse(3)
            for (i = 0; i < var; i++) {
                for (j = GC-1; j < rows-GC; j++) {
                    for (k = GC; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j-2+nn][k];
                        }
                        conserl[i][j][k] = WENO5_L(&fu[2]);
                        conserr[i][j][k] = WENO5_R(&fu[2]);
                    }
                }
            }
        }
        else if (dir == 2) {
            #pragma omp parallel for collapse(3)
            for (i = 0; i < var; i++) {
                for (j = GC; j < rows-GC; j++) {
                    for (k = GC-1; k < cols-GC; k++) {
                        double fu[6];
                        for (int nn = 0; nn < 6; nn++) {
                            fu[nn] = y[i][j][k-2+nn];
                        }
                        conserl[i][j][k] = WENO5_L(&fu[2]);
                        conserr[i][j][k] = WENO5_R(&fu[2]);
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
    a1 = 0.1 / pow(epsilon + s1, 2);
    a2 = 0.6 / pow(epsilon + s2, 2);
    a3 = 0.3 / pow(epsilon + s3, 2);

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
    a1 = 0.1 / pow(epsilon + s1, 2);
    a2 = 0.6 / pow(epsilon + s2, 2);
    a3 = 0.3 / pow(epsilon + s3, 2);

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