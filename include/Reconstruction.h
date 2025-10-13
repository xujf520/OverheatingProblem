#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <stdio.h>
#include <math.h>
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
/*                                         重构格式计算                                              */
/*                                      ******************                                          */



/*                                      ******************                                          */
/*                                                TVD                                               */
/*                                      ******************                                          */



static inline void Reconstruction_Godunov(int dir, int rows, int cols, int depth, int GC, double (*y)[cols][depth], \
                                            double (*conserl)[cols][depth], double (*conserr)[cols][depth], double delta_x) {
    int i, j, k;
    double epsilo = 1e-6;
    
    if (Characteriz) {
        // 动态分配内存
        double (*Pri)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Chara_Var)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Eigen_L)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*Eigen_R)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*W_L)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*W_R)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        // 检查内存分配是否成功
        if (Pri == NULL || Chara_Var == NULL || W_L == NULL || W_R == NULL) {
            fprintf(stderr, "Memory allocation failed in Reconstruction_Godunov\n");
            // 释放已分配的内存
            free(Pri);
            free(Chara_Var);
            free(Eigen_L);
            free(Eigen_R);
            free(W_L);
            free(W_R);
            return;
        }

        // 初始化数组
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                for (k = 0; k < depth; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    for (int ii = 0; ii < rows; ii++){
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }

        // 特征重构代码（暂时注释掉）
        Con_to_Pri_2D(rows,cols,depth,Pri,y);
        if (dir == 1){
            Compute_Eigen_2D(1.0, 0.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];
        


            for (i = 0; i < rows; i++) 
                 for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++){
                        W_L[i][j][k] = Chara_Var[i][j][k];
                        W_R[i][j][k] = Chara_Var[i][j+1][k];
                    }


            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                            
        }

        else if (dir == 2){
            Compute_Eigen_2D(0.0, 1.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];
        

            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++){
                        W_L[i][j][k] = Chara_Var[i][j][k];
                        W_R[i][j][k] = Chara_Var[i][j][k+1];
                    }


            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                            

        }
        // 释放动态分配的内存
        free(Pri);
        free(Chara_Var);
        free(Eigen_L);
        free(Eigen_R);
        free(W_L);
        free(W_R);
    } 
    else{
        if (dir == 1 )
            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++)
                    for ( k = GC; k < depth-GC; k++){
                        conserl[i][j][k] = y[i][j][k];
                        conserr[i][j][k] = y[i][j+1][k];
                    }
                  
        else if (dir == 2)
            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        conserl[i][j][k] = y[i][j][k];
                        conserr[i][j][k] = y[i][j][k+1];
                    }
                 
    }

}



static inline void TVD_Reconstruction(int dir, int rows, int cols, int depth, int GC, double (*y)[cols][depth], \
                                            double (*conserl)[cols][depth], double (*conserr)[cols][depth], double delta_x) {
    int i, j, k;
 

    if (Characteriz) {
        // 动态分配内存
        double (*Pri)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Chara_Var)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Eigen_L)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*Eigen_R)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*W_L)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*W_R)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        // 检查内存分配是否成功
        if (Pri == NULL || Chara_Var == NULL || W_L == NULL || W_R == NULL || Eigen_L==NULL || Eigen_R==NULL) {
            fprintf(stderr, "Memory allocation failed in Reconstruction_Godunov\n");
            // 释放已分配的内存
            free(Pri);
            free(Chara_Var);
            free(Eigen_L);
            free(Eigen_R);
            free(W_L);
            free(W_R);
            return;
        }

        // 初始化数组
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                for (k = 0; k < depth; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    for (int ii = 0; ii < rows; ii++){
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }

        // 特征重构代码（暂时注释掉）
        Con_to_Pri_2D(rows,cols,depth,Pri,y);
        if (dir == 1){
            Compute_Eigen_2D(1.0, 0.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];



            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++)
                    for ( k = GC; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = Chara_Var[i][j-2+nn][k];

                        W_L[i][j][k] = TVD_minmod_L(&fu[2],delta_x);
                        W_R[i][j][k] = TVD_minmod_R(&fu[2],delta_x);
                    }
        
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                            
        }

        else if (dir == 2){
            Compute_Eigen_2D(0.0, 1.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];


            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j][k-2+nn];

                        W_L[i][j][k] = TVD_minmod_L(&fu[2],delta_x);
                        W_R[i][j][k] = TVD_minmod_R(&fu[2],delta_x);
                    }

            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                            

        }
        // 释放动态分配的内存
        free(Pri);
        free(Chara_Var);
        free(Eigen_L);
        free(Eigen_R);
        free(W_L);
        free(W_R);
    } 

    else{
        if (dir == 1){
            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++)
                    for ( k = GC; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j-2+nn][k];

                        conserl[i][j][k] = TVD_minmod_L(&fu[2],delta_x);
                        conserr[i][j][k] = TVD_minmod_R(&fu[2],delta_x);
                    }
        }
        else if (dir == 2){
            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j][k-2+nn];

                        conserl[i][j][k] = TVD_minmod_L(&fu[2],delta_x);
                        conserr[i][j][k] = TVD_minmod_R(&fu[2],delta_x);
                    }
        }        

        
    }
}
/*                                      ******************                                          */
/*                                               WENO                                               */
/*                                      ******************                                          */


// 三阶WENO重构
static inline void WENO3_Reconstruction(int dir, int rows, int cols, int depth, int GC, double (*y)[cols][depth], \
                                            double (*conserl)[cols][depth], double (*conserr)[cols][depth]) {
    int i, j, k;

   
    if (Characteriz) {
         // 动态分配内存
        double (*Pri)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Chara_Var)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Eigen_L)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*Eigen_R)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*W_L)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*W_R)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        // 检查内存分配是否成功
        if (Pri == NULL || Chara_Var == NULL || W_L == NULL || W_R == NULL || Eigen_L==NULL || Eigen_R==NULL) {
            fprintf(stderr, "Memory allocation failed in Reconstruction_Godunov\n");
            // 释放已分配的内存
            free(Pri);
            free(Chara_Var);
            free(Eigen_L);
            free(Eigen_R);
            free(W_L);
            free(W_R);
            return;
        }

        // 初始化数组
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                for (k = 0; k < depth; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    for (int ii = 0; ii < rows; ii++){
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }

        // 特征重构代码（暂时注释掉）
        Con_to_Pri_2D(rows,cols,depth,Pri,y);
        if (dir == 1){
            Compute_Eigen_2D(1.0, 0.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];



            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++)
                    for ( k = GC; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = Chara_Var[i][j-2+nn][k];

                        W_L[i][j][k] = WENO3_L(&fu[2]);
                        W_R[i][j][k] = WENO3_R(&fu[2]);
                    }
        
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                            
        }

        else if (dir == 2){
            Compute_Eigen_2D(0.0, 1.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];


            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j][k-2+nn];

                        W_L[i][j][k] = WENO3_L(&fu[2]);
                        W_R[i][j][k] = WENO3_R(&fu[2]);
                    }

            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j <= cols-GC; j++) 
                    for ( k = GC-1; k <= depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                            

        }
        // 释放动态分配的内存
        free(Pri);
        free(Chara_Var);
        free(Eigen_L);
        free(Eigen_R);
        free(W_L);
        free(W_R);
    } 
    else{
        if (dir == 1){
            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++)
                    for ( k = GC; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j-2+nn][k];

                        conserl[i][j][k] = WENO3_L(&fu[2]);
                        conserr[i][j][k] = WENO3_R(&fu[2]);
                    }
        }
        else if (dir == 2){
            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j][k-2+nn];

                        conserl[i][j][k] = WENO3_L(&fu[2]);
                        conserr[i][j][k] = WENO3_R(&fu[2]);
                    }
                }        
    }

}


// 五阶WENO重构
static inline void WENO5_Reconstruction(int dir, int rows, int cols, int depth, int GC, double (*y)[cols][depth], \
                                            double (*conserl)[cols][depth], double (*conserr)[cols][depth]) {
    int i, j, k;

   
    if (Characteriz) {
         // 动态分配内存
        double (*Pri)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Chara_Var)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*Eigen_L)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*Eigen_R)[rows][cols][depth] = malloc(rows * sizeof(double[rows][cols][depth]));
        double (*W_L)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        double (*W_R)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
        // 检查内存分配是否成功
        if (Pri == NULL || Chara_Var == NULL || W_L == NULL || W_R == NULL || Eigen_L==NULL || Eigen_R==NULL) {
            fprintf(stderr, "Memory allocation failed in Reconstruction_Godunov\n");
            // 释放已分配的内存
            free(Pri);
            free(Chara_Var);
            free(Eigen_L);
            free(Eigen_R);
            free(W_L);
            free(W_R);
            return;
        }

        // 初始化数组
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                for (k = 0; k < depth; k++) {
                    Pri[i][j][k] = 0.0;
                    Chara_Var[i][j][k] = 0.0;
                    for (int ii = 0; ii < rows; ii++){
                        Eigen_L[i][ii][j][k] = 0.0;
                        Eigen_R[i][ii][j][k] = 0.0;
                    }
                    W_L[i][j][k] = 0.0;
                    W_R[i][j][k] = 0.0;
                    conserl[i][j][k] = 0.0;
                    conserr[i][j][k] = 0.0;
                }
            }
        }

        // 特征重构代码（暂时注释掉）
        Con_to_Pri_2D(rows,cols,depth,Pri,y);

        if (dir == 1){
            Compute_Eigen_2D(1.0, 0.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j < cols-GC; j++) 
                    for ( k = GC; k < depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];



            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++) 
                    for ( k = GC; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = Chara_Var[i][j-2+nn][k];

                        W_L[i][j][k] = WENO5_L(&fu[2]);
                        W_R[i][j][k] = WENO5_R(&fu[2]);
                    }
        
            for (i = 0; i < rows; i++) 
                for ( j = GC-1; j < cols-GC; j++) 
                    for ( k = GC; k < depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j+1][k];
                        }
                            
        }

        else if (dir == 2){
            Compute_Eigen_2D(0.0, 1.0, rows, cols,depth ,Pri, Eigen_L, Eigen_R);
            
            for (i = 0; i < rows; i++) 
                for ( j = GC; j < cols-GC; j++) 
                    for ( k = GC-1; k < depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++)
                            Chara_Var[i][j][k] += y[ii][j][k] * Eigen_L[i][ii][j][k];


            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = Chara_Var[i][j][k-2+nn];

                        W_L[i][j][k] = WENO5_L(&fu[2]);
                        W_R[i][j][k] = WENO5_R(&fu[2]);
                    }

            for (i = 0; i < rows; i++) 
                for ( j = GC; j < cols-GC; j++) 
                    for ( k = GC-1; k < depth-GC; k++)
                        for (int ii = 0; ii < rows; ii++){
                            conserl[i][j][k] += W_L[ii][j][k] * Eigen_R[i][ii][j][k];
                            conserr[i][j][k] += W_R[ii][j][k] * Eigen_R[i][ii][j][k+1];
                        }
                            

        }
        // 释放动态分配的内存
        free(Pri);
        free(Chara_Var);
        free(Eigen_L);
        free(Eigen_R);
        free(W_L);
        free(W_R);
    } 
    else{
        if (dir == 1){
            for ( i = 0; i < rows; i++)
                for ( j = GC-1; j < cols-GC; j++)
                    for ( k = GC; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j-2+nn][k];

                        conserl[i][j][k] = WENO5_L(&fu[2]);
                        conserr[i][j][k] = WENO5_R(&fu[2]);
                    }
        }
        else if (dir == 2){
            for ( i = 0; i < rows; i++)
                for ( j = GC; j < cols-GC; j++)
                    for ( k = GC-1; k < depth-GC; k++){
                        double fu[10] = {0.0};
                        for (int nn = 0; nn < 6; nn++)
                            fu[nn] = y[i][j][k-2+nn];

                        conserl[i][j][k] = WENO5_L(&fu[2]);
                        conserr[i][j][k] = WENO5_R(&fu[2]);
                    }
                }        
    }

}


/*                                      ******************                                          */
/*                                            重构函数                                              */
/*                                      ******************                                          */


//-------------------------------------------------------------------------------------------------
//       the 2th TVD Scheme
//  reference Appendix of paper: 
//A non-osillatory Eulerian Approch to interfaces in Multimaterial flows
//(Ghost Fluid Method) JCP vol 152, p457-492
//-------------------------------------------------------------------------------------------------

double TVD_minmod_L(double *f, double delta)
{

	int k;
	double v1, v2, v3;
	double slope1,slope2,slope;

	//assign value to v1, v2,...
	k = 0;
	v1 = *(f + k - 1);
	v2 = *(f + k);
	v3 = *(f + k + 1); 
	slope1 = (v2 - v1) / delta;
	slope2 = (v3 - v2) / delta;
	slope = min_mod(slope1,slope2);
	
	//reconstruction cell slope
	return  v2 + 0.5 * slope * delta;
}

double TVD_minmod_R(double *f, double delta)
{

	int k;
	double v1, v2, v3;
	double slope1,slope2,slope;

	//assign value to v1, v2,...
	k = 1;
	v1 = *(f + k - 1);
	v2 = *(f + k);
	v3 = *(f + k + 1); 
	slope1 = (v2 - v1) / delta;
	slope2 = (v3 - v2) / delta;
	slope = min_mod(slope1,slope2);

	return  v2 - 0.5 * slope * delta;
}

double TVD_vanleer_L(double *f, double delta)
{

	int k;
	double v1, v2, v3;
	double slope1,slope2,slope;

	//assign value to v1, v2,...
	k = 0;
	v1 = *(f + k - 1);
	v2 = *(f + k);
	v3 = *(f + k + 1); 
	slope1 = (v2 - v1) / delta;
	slope2 = (v3 - v2) / delta;
	slope = van_leer(slope1,slope2);

	//reconstruction cell slope
	return  v2 + 0.5 * slope * delta;
}

double TVD_vanleer_R(double *f, double delta)
{

	int k;
	double v1, v2, v3;
	double slope1,slope2,slope;

	//assign value to v1, v2,...
	k = 1;
	v1 = *(f + k - 1);
	v2 = *(f + k);
	v3 = *(f + k + 1); 
	slope1 = (v2 - v1) / delta;
	slope2 = (v3 - v2) / delta;
	slope = van_leer(slope1,slope2);

	return  v2 - 0.5 * slope * delta;
}



//-------------------------------------------------------------------------------------------------
//  the  WENO Scheme
//  reference Appendix of paper: 
//A non-osillatory Eulerian Approch to interfaces in Multimaterial flows
//(Ghost Fluid Method) JCP vol 152, p457-492
//-------------------------------------------------------------------------------------------------

static inline double WENO3_L(double *f)
{

	int k;
	double v1, v2, v3;
	double s1, s2;
	double a1, a2, a3, w1, w2, w3;
    double epsilo = 1e-6;

	//assign value to v1, v2,...
	k = 0;
	v1 = *(f + k - 1);
	v2 = *(f + k);
	v3 = *(f + k + 1);

	//smoothness indicator
	s1 = (v2-v1) * (v2-v1);
	s2 = (v3-v2) * (v3-v2);

	//weights
	a1 = (1./3)/pow(epsilo + s1, 2);
	a2 = (2./3)/pow(epsilo + s2, 2);

	w1 = a1/(a1 + a2);
	w2 = a2/(a1 + a2);

    if ( w1 < 0.0 || w2 < 0.0)
		printf("Negative weights appear in WENO!!!\n");
	
	//return weighted average
	return  w1 * (-0.5*v1 + 1.5*v2) + w2 * (0.5*v2 + 0.5*v3);
}



static inline double WENO3_R(double *f)
{

	int k;
	double v1, v2, v3;
	double s1, s2;
	double a1, a2,a3, w1, w2, w3;
    double epsilon = 1e-6;

	//assign value to v1, v2,...
	k = 1;
	v1 = *(f + k + 1);
	v2 = *(f + k );
	v3 = *(f + k - 1);

	//smoothness indicator
	s1 = (v2-v1) * (v2-v1);
	s2 = (v3-v2) * (v3-v2);

	//weights
	a1 = (1./3)/pow(epsilon + s1, 2);
	a2 = (2./3)/pow(epsilon + s2, 2);

	w1 = a1/(a1 + a2);
	w2 = a2/(a1 + a2);

	if ( w1 < 0.0 || w2 < 0.0)
		printf("Negative weights appear in WENO!!!\n");
	//return weighted average
	return  w1 * (-0.5*v1 + 1.5*v2) + w2 * (0.5*v2 + 0.5*v3);

}




static inline double WENO5_L(double *f){

	int k;
	double v1, v2, v3, v4, v5;
	double s1, s2, s3;
	double a1, a2, a3, w1, w2, w3;
    double epsilon = 1.0e-6;

	//assign value to v1, v2,...
	k = 0;
	v1 = *(f + k - 2);
	v2 = *(f + k - 1);
	v3 = *(f + k);
	v4 = *(f + k + 1); 
	v5 = *(f + k + 2);

	//smoothness indicator
	s1 = 13.0/12.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) 
	   + 0.25*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	s2 = 13.0/12.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) 
	   + 0.25*(v2 - v4)*(v2 - v4);
	s3 = 13.0/12.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) 
	   + 0.25*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);

	//weights
	a1 = 0.1/pow(epsilon + s1, 2);
	a2 = 0.6/pow(epsilon + s2, 2);
	a3 = 0.3/pow(epsilon + s3, 2);

	w1 = a1/(a1 + a2 + a3);
	w2 = a2/(a1 + a2 + a3);
	w3 = a3/(a1 + a2 + a3);
    if ( w1 < 0.0 || w2 < 0.0 || w3 < 0.0 )
		printf("Negative weights appear in WENO!!!\n");
	

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3)/6.0
		  + w2*(-v2 + 5.0*v3 + 2.0*v4)/6.0
		  + w3*(2.0*v3 + 5.0*v4 - v5)/6.0;

}

static inline double WENO5_R(double *f)
{

	int k;
	double v1, v2, v3, v4, v5;
	double s1, s2,s3;
	double a1, a2,a3, w1, w2, w3;
    double epsilon = 1.0e-6;

	//assign value to v1, v2,...
	k = 1;
	v1 = *(f + k + 2);
	v2 = *(f + k + 1);
	v3 = *(f + k);
	v4 = *(f + k - 1); 
	v5 = *(f + k - 2);

	//smoothness indicator
	s1 = 13.0/12.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) 
	   + 0.25*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	s2 = 13.0/12.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) 
	   + 0.25*(v2 - v4)*(v2 - v4);
	s3 = 13.0/12.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) 
	   + 0.25*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);

	//weights
	a1 = 0.1/pow(epsilon + s1, 2);
	a2 = 0.6/pow(epsilon + s2, 2);
	a3 = 0.3/pow(epsilon + s3, 2);

	w1 = a1/(a1 + a2 +a3);
	w2 = a2/(a1 + a2 +a3);
	w3 = a3/(a1 + a2 +a3);
	if ( w1 < 0.0 || w2 < 0.0 || w3 < 0.0 )
		printf("Negative weights appear in WENO!!!\n");
	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3)/6.0
		  + w2*(-v2 + 5.0*v3 + 2.0*v4)/6.0
		  + w3*(2.0*v3 + 5.0*v4 - v5)/6.0;

}



#endif  