#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <stdio.h>
#include <math.h>
#include "function.h"
#include "Characteriz.h"
#include "Golbal.h"



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



static inline void Reconstruction_Godunov(int rows, int cols, int GC, double (*y)[cols], double (*conserl)[cols], double (*conserr)[cols], double delta_x) {
    int i, j;
    double epsilo = 1e-6;

    double slope[rows][cols];
    double Pri[rows][cols],Chara_Var[rows][cols];
    double W_L[rows][cols],W_R[rows][cols];
    double eigen_l[rows][rows][cols],eigen_r[rows][rows][cols];
 
    //初始化数组
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            Chara_Var[i][j] = 0.0;
            W_L[i][j] = 0.0;
            W_R[i][j] = 0.0;
            conserl[i][j] = 0.0;
            conserr[i][j] = 0.0;
        }
    }

    if (Characteriz){
        //计算特征矩阵
        Con_to_Pri_1D(3,cols,Pri,y);
        Compute_Eigen_Matrix(3, rows, cols, Pri, eigen_l, eigen_r);
        //投影到特征空间，计算特征变量：
        for (int k = 0; k < rows; k++){
            for (int j = 0; j < cols; j++){
                for (int m = 0; m < rows; m++){
                    Chara_Var[k][j] += y[m][j] * eigen_l[k][m][j];
                }
            }   
        }

        for ( i = 0; i < rows; i++){
            for ( j = GC-1; j <= cols-GC; j++){
                W_L[i][j] = Chara_Var[i][j];
                W_R[i][j] = Chara_Var[i][j+1];
            }
        }

        for (int k = 0; k < rows; k++){
            for (int j = GC-1; j <= cols-GC; j++){
                for (int m = 0; m < rows; m++){
                    conserl[k][j] += W_L[m][j] * eigen_r[k][m][j];
                    conserr[k][j] += W_R[m][j] * eigen_r[k][m][j+1];
                }
            }  
        }
    }
    else{
        for (int i = 0; i < rows; i++){
            for (int j = GC-1; j <= cols-GC-1; j++){
                    conserl[i][j] = y[i][j];
                    conserr[i][j] = y[i][j+1];
            }  
        }
    }

}

static inline void TVD_Reconstruction(int rows, int cols, int GC,double (*y)[cols], double (*conserl)[cols], double (*conserr)[cols], double delta_x) {
    int i, j;
    double epsilo = 1e-6;

    double slope[rows][cols],slope_a[rows][cols],slope_b[rows][cols];
    double Pri[rows][cols],Chara_Var[rows][cols];
    double W_L[rows][cols],W_R[rows][cols];
    double eigen_l[rows][rows][cols],eigen_r[rows][rows][cols];
 
    //初始化数组
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            slope[i][j] = 0.0;
            slope_a[i][j] = 0.0;
            slope_b[i][j] = 0.0;
            Chara_Var[i][j] = 0.0;
            W_L[i][j] = 0.0;
            W_R[i][j] = 0.0;
            conserl[i][j] = 0.0;
            conserr[i][j] = 0.0;
        }
    }

    if (Characteriz){
        Con_to_Pri_1D(3,cols,Pri,y);
        //计算特征矩阵
        Compute_Eigen_Matrix(3, rows, cols, Pri, eigen_l, eigen_r);
        //投影到特征空间，计算特征变量：
        for (int k = 0; k < rows; k++){
            for (int j = 0; j < cols; j++){
                for (int m = 0; m < rows; m++){
                    Chara_Var[k][j] += y[m][j] * eigen_l[k][m][j];
                }
            }   
        }

        //重构特征变量
        for ( i = 0; i < rows; i++){
            for ( j = GC-1; j < cols-GC; j++){
                slope_a[i][j] = (Chara_Var[i][j] - Chara_Var[i][j-1])/delta_x;
                slope_b[i][j] = (Chara_Var[i][j+1] - Chara_Var[i][j])/delta_x;
            }
        }

        //重构
        for (int i = 0; i < rows; i++){   
            for (int j = GC-1; j < cols-GC; j++)
            {
                //vanleer
                //slope[i][j]= ((sgn(a[i][j])+sgn(b[i][j]))*a[i][j]*b[i][j])/(fabs(a[i][j])+fabs(b[i][j])+1e-15);
                //minibee
                slope[i][j] = 0.5 * (sgn(slope_a[i][j])+sgn(slope_b[i][j])) * min_of_two(fabs(slope_a[i][j]), fabs(slope_b[i][j]));
                //vanalbada
                //slope[i][j]=(fmax(a[i][j]*b[i][j],0) * (a[i][j]+b[i][j]))/(pow(a[i][j],2)+pow(b[i][j],2)+10e-6);
            }
        }


        for ( i = 0; i < rows; i++){
            for ( j = GC-1; j <= cols-GC; j++){
                W_L[i][j] = Chara_Var[i][j] + 0.5*slope[i][j] * delta_x;
                W_R[i][j] = Chara_Var[i][j+1] - 0.5*slope[i][j+1] * delta_x;
            }
        }

        for (int k = 0; k < rows; k++){
            for (int j = GC-1; j <= cols-GC-1; j++){
                for (int m = 0; m < rows; m++){
                    conserl[k][j] += W_L[m][j] * eigen_r[k][m][j];
                    conserr[k][j] += W_R[m][j] * eigen_r[k][m][j+1];
                }
            }  
        }
    }
    else{
        //重构变量
        for ( i = 0; i < rows; i++){
            for ( j = GC-1; j < cols-GC; j++){
                slope_a[i][j] = (y[i][j] - y[i][j-1])/delta_x;
                slope_b[i][j] = (y[i][j+1] - y[i][j])/delta_x;
            }
        }

        //重构
        for (int i = 0; i < rows; i++){   
            for (int j = GC-1; j < cols-GC; j++)
            {
                //vanleer
                slope[i][j]= ((sgn(slope_a[i][j])+sgn(slope_b[i][j]))*slope_a[i][j]*slope_b[i][j])/(fabs(slope_a[i][j])+fabs(slope_b[i][j])+1e-15);
                //minibee
                //slope[i][j] = 0.5 * (sgn(slope_a[i][j])+sgn(slope_b[i][j])) * min_of_two(fabs(slope_a[i][j]), fabs(slope_b[i][j]));
                //vanalbada
                //slope[i][j]=(fmax(a[i][j]*b[i][j],0) * (a[i][j]+b[i][j]))/(pow(a[i][j],2)+pow(b[i][j],2)+10e-6);
            }
        }

        for ( i = 0; i < rows; i++){
            for ( j = GC-1; j <= cols-GC; j++){
                conserl[i][j] = y[i][j] + 0.5*slope[i][j] * delta_x;
                conserr[i][j] = y[i][j+1] - 0.5*slope[i][j+1] * delta_x;
            }
        }

    }
    
    

}


/*                                      ******************                                          */
/*                                               WENO                                               */
/*                                      ******************                                          */


// 三阶WENO重构
static inline void WENO3_Reconstruction(int rows, int cols, int GC,double (*y)[cols], double (*conserl)[cols], double (*conserr)[cols]) {
    int i, j;
    double epsilo = 1e-8;
    
    // 声明局部变量
    double beta_0[rows][cols], beta_1[rows][cols];
    double alphal_0[rows][cols], alphal_1[rows][cols];
    double alphar_0[rows][cols], alphar_1[rows][cols];
    double weightl_0[rows][cols], weightl_1[rows][cols];
    double weightr_0[rows][cols], weightr_1[rows][cols];

    double Pri[rows][cols],Chara_Var[rows][cols];
    double W_L[rows][cols],W_R[rows][cols];
    double eigen_l[rows][rows][cols],eigen_r[rows][rows][cols];
 
    //初始化数组
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            Chara_Var[i][j] = 0.0;
            W_L[i][j] = W_R[i][j] = 0.0;
            beta_0[i][j] = beta_1[i][j] = 0.0;
            alphal_0[i][j] = alphal_1[i][j] = 0.0;
            alphar_0[i][j] = alphar_1[i][j] = 0.0;
            weightl_0[i][j] = weightl_1[i][j] = 0.0;
            weightr_0[i][j] = weightr_1[i][j] = 0.0;
            conserl[i][j] = conserr[i][j] = 0.0;
        }
    }

    if (Characteriz){
        Con_to_Pri_1D(rows,cols,Pri,y);
        //计算特征矩阵
        Compute_Eigen_Matrix(rows, rows, cols, Pri, eigen_l, eigen_r);
        //投影到特征空间，计算特征变量：
        for (int k = 0; k < rows; k++){
            for (int j = 0; j < cols; j++){
                for (int m = 0; m < rows; m++){
                    Chara_Var[k][j] += y[m][j] * eigen_l[k][m][j];
                }
            }   
        }


        for (i = 0; i < rows; i++) {
            for (j = GC-1; j <= cols-GC; j++) {

                //计算间断因子
                beta_0[i][j] = pow(Chara_Var[i][j+1] - Chara_Var[i][j], 2);
                beta_1[i][j] = pow(Chara_Var[i][j] - Chara_Var[i][j-1], 2);
                
                //计算非线性的系数
                alphal_0[i][j] = (2.0 / 3.0) * (1.0 / pow(epsilo + beta_0[i][j], 2));
                alphal_1[i][j] = (1.0 / 3.0) * (1.0 / pow(epsilo + beta_1[i][j], 2));
                alphar_0[i][j] = (1.0 / 3.0) * (1.0 / pow(epsilo + beta_0[i][j], 2));
                alphar_1[i][j] = (2.0 / 3.0) * (1.0 / pow(epsilo + beta_1[i][j], 2));
            
                //计算非线性权
                weightl_0[i][j] = alphal_0[i][j] / (alphal_0[i][j] + alphal_1[i][j]);
                weightl_1[i][j] = alphal_1[i][j] / (alphal_0[i][j] + alphal_1[i][j]);
                weightr_0[i][j] = alphar_0[i][j] / (alphar_0[i][j] + alphar_1[i][j]);
                weightr_1[i][j] = alphar_1[i][j] / (alphar_0[i][j] + alphar_1[i][j]);
            }
        
            for (j = GC-1; j <= cols-GC; j++) {
                W_L[i][j] = weightl_0[i][j] * (0.5 * Chara_Var[i][j] + 0.5 * Chara_Var[i][j+1]) \
                            + weightl_1[i][j] * (-0.5 * Chara_Var[i][j-1] + 1.5 * Chara_Var[i][j]);
                W_R[i][j] = weightr_0[i][j+1] * (1.5 * Chara_Var[i][j+1] - 0.5 * Chara_Var[i][j+2]) \
                            + weightr_1[i][j+1] * (0.5 * Chara_Var[i][j] + 0.5 * Chara_Var[i][j+1]);
            }
        }

        for (int k = 0; k < rows; k++){
            for (int j = GC-1; j <= cols-GC-1; j++){
                for (int m = 0; m < rows; m++){
                    conserl[k][j] += W_L[m][j] * eigen_r[k][m][j];
                    conserr[k][j] += W_R[m][j] * eigen_r[k][m][j+1];
                }
            }  
        }
    }

    else{

        for (i = 0; i < rows; i++) {
            for (j = GC-1; j <= cols-GC; j++) {
                //计算间断因子
                beta_0[i][j] = pow(y[i][j+1] - y[i][j], 2.0);
                beta_1[i][j] = pow(y[i][j] - y[i][j-1], 2.0);
            }   
        }

        for (i = 0; i < rows; i++) {
            for (j = GC-1; j <= cols-GC; j++) {                
                //计算非线性的系数
                alphal_0[i][j] = (2.0 / 3.0) * (1.0 / pow(epsilo + beta_0[i][j], 2));
                alphal_1[i][j] = (1.0 / 3.0) * (1.0 / pow(epsilo + beta_1[i][j], 2));
                alphar_0[i][j] = (1.0 / 3.0) * (1.0 / pow(epsilo + beta_0[i][j], 2));
                alphar_1[i][j] = (2.0 / 3.0) * (1.0 / pow(epsilo + beta_1[i][j], 2));
            }
        }

        for (i = 0; i < rows; i++) {
            for (j = GC-1; j <= cols-GC; j++) {
                //计算非线性权
                weightl_0[i][j] = alphal_0[i][j] / (alphal_0[i][j] + alphal_1[i][j]);
                weightl_1[i][j] = alphal_1[i][j] / (alphal_0[i][j] + alphal_1[i][j]);
                weightr_0[i][j] = alphar_0[i][j] / (alphar_0[i][j] + alphar_1[i][j]);
                weightr_1[i][j] = alphar_1[i][j] / (alphar_0[i][j] + alphar_1[i][j]);
            }
        
            
        }

        for (i = 0; i < rows; i++){
           for (j = GC-1; j <= cols-GC; j++) {
                conserl[i][j] = weightl_0[i][j] * (0.5 * y[i][j] + 0.5 * y[i][j+1])\
                                + weightl_1[i][j] * (-0.5 * y[i][j-1] + 1.5 * y[i][j]);
                conserr[i][j] = weightr_0[i][j+1] * (1.5 * y[i][j+1] - 0.5 * y[i][j+2])\
                                + weightr_1[i][j+1] * (0.5 * y[i][j] + 0.5 * y[i][j+1]);
            }
        }

    }
    
}

// 五阶WENO重构
static inline void WENO5_Reconstruction(int rows, int cols, int GC, double (*y)[cols], double (*conserl)[cols], double (*conserr)[cols]) {
    int i, j;
    double epsilo = 1e-6;
    
    // 声明局部变量
    double beta_0[rows][cols], beta_1[rows][cols], beta_2[rows][cols];
    double alphal_0[rows][cols], alphal_1[rows][cols], alphal_2[rows][cols];
    double alphar_0[rows][cols], alphar_1[rows][cols], alphar_2[rows][cols];
    double weightl_0[rows][cols], weightl_1[rows][cols], weightl_2[rows][cols];
    double weightr_0[rows][cols], weightr_1[rows][cols], weightr_2[rows][cols];

    // 初始化数组
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            beta_0[i][j] = beta_1[i][j] = beta_2[i][j] = 0.0;
            alphal_0[i][j] = alphal_1[i][j] = alphal_2[i][j] = 0.0;
            alphar_0[i][j] = alphar_1[i][j] = alphar_2[i][j] = 0.0;
            weightl_0[i][j] = weightl_1[i][j] = weightl_2[i][j] = 0.0;
            weightr_0[i][j] = weightr_1[i][j] = weightr_2[i][j] = 0.0;
            conserl[i][j] = conserr[i][j] = 0.0;
        }
    }


    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            beta_0[i][j] = (13.0 / 12.0) * pow(y[i][j] - 2.0 * y[i][j+1] + y[i][j+2], 2)
                          + (1.0 / 4.0) * pow(3.0 * y[i][j] - 4.0 * y[i][j+1] + y[i][j+2], 2);
            beta_1[i][j] = (13.0 / 12.0) * pow(y[i][j-1] - 2.0 * y[i][j] + y[i][j+1], 2)
                          + (1.0 / 4.0) * pow(y[i][j-1] - y[i][j+1], 2);
            beta_2[i][j] = (13.0 / 12.0) * pow(y[i][j-2] - 2.0 * y[i][j-1] + y[i][j], 2)
                          + (1.0 / 4.0) * pow(y[i][j-2] - 4.0 * y[i][j-1] + 3.0 * y[i][j], 2);
        }
    }
    
    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            alphal_0[i][j] = 0.3 * (1.0 / pow(epsilo + beta_0[i][j], 2));
            alphal_1[i][j] = 0.6 * (1.0 / pow(epsilo + beta_1[i][j], 2));
            alphal_2[i][j] = 0.1 * (1.0 / pow(epsilo + beta_2[i][j], 2));
            alphar_0[i][j] = 0.1 * (1.0 / pow(epsilo + beta_0[i][j], 2));
            alphar_1[i][j] = 0.6 * (1.0 / pow(epsilo + beta_1[i][j], 2));
            alphar_2[i][j] = 0.3 * (1.0 / pow(epsilo + beta_2[i][j], 2));
        }
    }
    
    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            weightl_0[i][j] = alphal_0[i][j] / (alphal_0[i][j] + alphal_1[i][j] + alphal_2[i][j]);
            weightl_1[i][j] = alphal_1[i][j] / (alphal_0[i][j] + alphal_1[i][j] + alphal_2[i][j]);
            weightl_2[i][j] = alphal_2[i][j] / (alphal_0[i][j] + alphal_1[i][j] + alphal_2[i][j]);
            weightr_0[i][j] = alphar_0[i][j] / (alphar_0[i][j] + alphar_1[i][j] + alphar_2[i][j]);
            weightr_1[i][j] = alphar_1[i][j] / (alphar_0[i][j] + alphar_1[i][j] + alphar_2[i][j]);
            weightr_2[i][j] = alphar_2[i][j] / (alphar_0[i][j] + alphar_1[i][j] + alphar_2[i][j]);
        }
    }
    
    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC-1; j++) {
            conserl[i][j] = weightl_0[i][j] * ((1.0/3.0)*y[i][j] + (5.0/6.0)*y[i][j+1] - (1.0/6.0)*y[i][j+2]) 
                          + weightl_1[i][j] * (-(1.0/6.0)*y[i][j-1] + (5.0/6.0)*y[i][j] + (1.0/3.0)*y[i][j+1])
                          + weightl_2[i][j] * ((1.0/3.0)*y[i][j-2] - (7.0/6.0)*y[i][j-1] + (11.0/6.0)*y[i][j]);

            conserr[i][j] = weightr_0[i][j+1] * ((11.0/6.0)*y[i][j+1] - (7.0/6.0)*y[i][j+2] + (1.0/3.0)*y[i][j+3]) 
                            + weightr_1[i][j+1] * ((1.0/3.0)*y[i][j] + (5.0/6.0)*y[i][j+1] - (1.0/6.0)*y[i][j+2])
                            + weightr_2[i][j+1] * (-(1.0/6.0)*y[i][j-1] + (5.0/6.0)*y[i][j] + (1.0/3.0)*y[i][j+1]);
        }
    }
}


// 五阶WENO重构
static inline void WENO5_Reconstruction_C(int rows, int cols, int GC, double (*y)[cols], double (*conserl)[cols], double (*conserr)[cols]){
    int i, j;
    double epsilo = 1e-6;
    
    // 声明局部变量
   double beta_0[rows][cols], beta_1[rows][cols], beta_2[rows][cols];
    double alphal_0[rows][cols], alphal_1[rows][cols], alphal_2[rows][cols];
    double alphar_0[rows][cols], alphar_1[rows][cols], alphar_2[rows][cols];
    double weightl_0[rows][cols], weightl_1[rows][cols], weightl_2[rows][cols];
    double weightr_0[rows][cols], weightr_1[rows][cols], weightr_2[rows][cols];
    double Pri[rows][cols],Chara_Var[rows][cols];
    double W_L[rows][cols],W_R[rows][cols];
    double eigen_l[rows][rows][cols],eigen_r[rows][rows][cols];
 

    // 初始化数组
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            Chara_Var[i][j] = 0.0;
            W_L[i][j] = W_R[i][j] = 0.0;
            beta_0[i][j] = beta_1[i][j] = beta_2[i][j] = 0.0;
            alphal_0[i][j] = alphal_1[i][j] = alphal_2[i][j] = 0.0;
            alphar_0[i][j] = alphar_1[i][j] = alphar_2[i][j] = 0.0;
            weightl_0[i][j] = weightl_1[i][j] = weightl_2[i][j] = 0.0;
            weightr_0[i][j] = weightr_1[i][j] = weightr_2[i][j] = 0.0;
            conserl[i][j] = conserr[i][j] = 0.0;
        }
    }

    Con_to_Pri_1D(3,cols,Pri,y);
   //计算特征矩阵
    if (Characteriz)
        Compute_Eigen_Matrix(3, rows, cols, Pri, eigen_l, eigen_r);
    //投影到特征空间，计算特征变量：
    for (int k = 0; k < rows; k++){
        for (int j = 0; j < cols; j++){
            for (int m = 0; m < rows; m++){
                Chara_Var[k][j] += y[m][j] * eigen_l[k][m][j];
            }
        }   
    }

    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            beta_0[i][j] = (13.0 / 12.0) * pow(Chara_Var[i][j] - 2.0 * Chara_Var[i][j+1] + Chara_Var[i][j+2], 2)
                          + (1.0 / 4.0) * pow(3.0 * Chara_Var[i][j] - 4.0 * Chara_Var[i][j+1] + Chara_Var[i][j+2], 2);
            beta_1[i][j] = (13.0 / 12.0) * pow(Chara_Var[i][j-1] - 2.0 * Chara_Var[i][j] + Chara_Var[i][j+1], 2)
                          + (1.0 / 4.0) * pow(Chara_Var[i][j-1] - Chara_Var[i][j+1], 2);
            beta_2[i][j] = (13.0 / 12.0) * pow(Chara_Var[i][j-2] - 2.0 * Chara_Var[i][j-1] + Chara_Var[i][j], 2)
                          + (1.0 / 4.0) * pow(Chara_Var[i][j-2] - 4.0 * Chara_Var[i][j-1] + 3.0 * Chara_Var[i][j], 2);
        }
    }
    
    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            alphal_0[i][j] = 0.3 * (1.0 / pow(epsilo + beta_0[i][j], 2));
            alphal_1[i][j] = 0.6 * (1.0 / pow(epsilo + beta_1[i][j], 2));
            alphal_2[i][j] = 0.1 * (1.0 / pow(epsilo + beta_2[i][j], 2));
            alphar_0[i][j] = 0.1 * (1.0 / pow(epsilo + beta_0[i][j], 2));
            alphar_1[i][j] = 0.6 * (1.0 / pow(epsilo + beta_1[i][j], 2));
            alphar_2[i][j] = 0.3 * (1.0 / pow(epsilo + beta_2[i][j], 2));
        }
    }
    
    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            weightl_0[i][j] = alphal_0[i][j] / (alphal_0[i][j] + alphal_1[i][j] + alphal_2[i][j]);
            weightl_1[i][j] = alphal_1[i][j] / (alphal_0[i][j] + alphal_1[i][j] + alphal_2[i][j]);
            weightl_2[i][j] = alphal_2[i][j] / (alphal_0[i][j] + alphal_1[i][j] + alphal_2[i][j]);
            weightr_0[i][j] = alphar_0[i][j] / (alphar_0[i][j] + alphar_1[i][j] + alphar_2[i][j]);
            weightr_1[i][j] = alphar_1[i][j] / (alphar_0[i][j] + alphar_1[i][j] + alphar_2[i][j]);
            weightr_2[i][j] = alphar_2[i][j] / (alphar_0[i][j] + alphar_1[i][j] + alphar_2[i][j]);
        }
    }
    
    for (i = 0; i < rows; i++) {
        for (j = GC-1; j <= cols-GC; j++) {
            W_L[i][j] = weightl_0[i][j] * ((1.0/3.0)*Chara_Var[i][j] + (5.0/6.0)*Chara_Var[i][j+1] - (1.0/6.0)*Chara_Var[i][j+2]) 
                          + weightl_1[i][j] * (-(1.0/6.0)*Chara_Var[i][j-1] + (5.0/6.0)*Chara_Var[i][j] + (1.0/3.0)*Chara_Var[i][j+1])
                          + weightl_2[i][j] * ((1.0/3.0)*Chara_Var[i][j-2] - (7.0/6.0)*Chara_Var[i][j-1] + (11.0/6.0)*Chara_Var[i][j]);

            W_R[i][j] = weightr_0[i][j+1] * ((11.0/6.0)*Chara_Var[i][j+1] - (7.0/6.0)*Chara_Var[i][j+2] + (1.0/3.0)*Chara_Var[i][j+3]) 
                            + weightr_1[i][j+1] * ((1.0/3.0)*Chara_Var[i][j] + (5.0/6.0)*Chara_Var[i][j+1] - (1.0/6.0)*Chara_Var[i][j+2])
                            + weightr_2[i][j+1] * (-(1.0/6.0)*Chara_Var[i][j-1] + (5.0/6.0)*Chara_Var[i][j] + (1.0/3.0)*Chara_Var[i][j+1]);
        }
    }

    for (int k = 0; k < rows; k++){
        for (int j = GC-1; j <= cols-GC-1; j++){
            for (int m = 0; m < rows; m++){
                conserl[k][j] += W_L[m][j] * eigen_r[k][m][j];
                conserr[k][j] += W_R[m][j] * eigen_r[k][m][j+1];
            }
        }  
    }
}








//计算有限体积格式误差
static inline void Scheme_Error(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    
    int l =  cols- 2*GC;

    for (int i = 0; i < rows; i++)
    {
        Error[i] = 0.0;
    }
    

    WENO3_Reconstruction(rows,cols,GC,x,Error_L,Error_R);

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


#endif  