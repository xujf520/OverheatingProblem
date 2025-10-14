#ifndef FLUX_H
#define FLUX_H

#include <stdio.h>
#include <math.h>
#include "Characteriz.h"
#include "Golbal.h"

// 定义全局变量来累加 p_star - 2.926650
static int Ite = 0;
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
//Flux计算方法


//Flux计算方法
static inline void HLL_Flux(int dir, int rows, int cols, int depth, int GC, double (*x)[cols][depth], double (*y)[cols][depth] ,double (*z)[cols][depth]) {


    if (dir == 1){
       for (int j = GC-1; j < cols-GC; j++) {
            for (int k = GC; k < depth-GC; k++){
                //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);
                
                //计算声速
                double a_L = sqrt(M_gamma * p_L / rho_L);
                double a_R = sqrt(M_gamma * p_R / rho_R);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量

                double rho_FL = rho_L * u_L;
                double rho_FR = rho_R * u_R;
                double rhou_FL = rho_L * u_L * u_L + p_L;
                double rhou_FR = rho_R * u_R * u_R + p_R;
                double rhov_FL = rho_L * v_L * u_L;
                double rhov_FR = rho_R * v_R * u_R;
                double rhoe_FL = rho_L * H_L * u_L;
                double rhoe_FR = rho_R * H_R * u_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = ubar - cbar;
                double sright = ubar + cbar;

                //确定HLL数值通量
                double rho_F = 0, rhou_F = 0, rhov_F = 0,rhoe_F = 0;
                if (sleft >= 0 ){
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if(sleft < 0 && sright > 0){
                    rho_F = (sright*rho_FL - sleft*rho_FR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
                    rhou_F = (sright*rhou_FL - sleft*rhou_FR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
                    rhov_F = (sright*rhov_FL - sleft*rhov_FR + sleft*sright * (rhov_R - rhov_L))/(sright-sleft);
                    rhoe_F = (sright*rhoe_FL - sleft*rhoe_FR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
                }
                else if (sright <= 0){
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }

                z[0][j][k] = rho_F; 
                z[1][j][k] = rhou_F; 
                z[2][j][k] = rhov_F;
                z[3][j][k] = rhoe_F; 
            }
        }   
    }

    else if (dir == 2){
        for (int j = GC; j < cols-GC; j++) {
            for (int k = GC-1; k < depth-GC; k++){
               //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量

                double rho_GL = rho_L * v_L;
                double rho_GR = rho_R * v_R;
                double rhou_GL = rho_L * u_L * v_L;
                double rhou_GR = rho_R * u_R * v_R;
                double rhov_GL = rho_L * v_L*v_L + p_L;
                double rhov_GR = rho_R * v_R*v_R + p_R;
                double rhoe_GL = rho_L * H_L * v_L;
                double rhoe_GR = rho_R * H_R * v_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = vbar - cbar;
                double sright = vbar + cbar;

                //确定HLL数值通量
                double rho_G = 0, rhou_G = 0, rhov_G = 0,rhoe_G = 0;
                if (sleft >= 0 ){
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if(sleft < 0 && sright >0){
                    rho_G = (sright*rho_GL - sleft*rho_GR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
                    rhou_G = (sright*rhou_GL - sleft*rhou_GR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
                    rhov_G = (sright*rhov_GL - sleft*rhov_GR + sleft*sright * (rhov_R - rhov_L))/(sright-sleft);
                    rhoe_G = (sright*rhoe_GL - sleft*rhoe_GR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
                }

                else if (sright <= 0){
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }

                z[0][j][k] = rho_G; 
                z[1][j][k] = rhou_G; 
                z[2][j][k] = rhov_G;
                z[3][j][k] = rhoe_G; 


            }
        }
    }
    
}



static inline void HLLC_Flux(int dir, int rows, int cols, int depth, int GC, double (*x)[cols][depth], double (*y)[cols][depth] ,double (*z)[cols][depth]) {


    if (dir == 1){
       for (int j = GC-1; j < cols-GC; j++) {
            for (int k = GC; k < depth-GC; k++){
                //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);
                
                //计算声速
                double a_L = sqrt(M_gamma * p_L / rho_L);
                double a_R = sqrt(M_gamma * p_R / rho_R);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量

                double rho_FL = rho_L * u_L;
                double rho_FR = rho_R * u_R;
                double rhou_FL = rho_L * u_L * u_L + p_L;
                double rhou_FR = rho_R * u_R * u_R + p_R;
                double rhov_FL = rho_L * v_L * u_L;
                double rhov_FR = rho_R * v_R * u_R;
                double rhoe_FL = rho_L * H_L * u_L;
                double rhoe_FR = rho_R * H_R * u_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = ubar - cbar;
                double sright = ubar + cbar;
            
                double s_star = (p_R - p_L + rho_L*u_L*(sleft - u_L) - rho_R*u_R*(sright - u_R))\
                                        /(rho_L*(sleft - u_L) - rho_R*(sright - u_R));
                                        
                double p_star = p_L + rho_L * (sleft-u_L) * (s_star-u_L);
                double u_stat_L = rho_L * (sleft-u_L)/(sleft-s_star);
                double u_stat_R = rho_R * (sright-u_R)/(sright-s_star);

                double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L)));
                double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R)));
                //确定HLLC数值通量
                double rho_F = 0, rhou_F = 0, rhov_F = 0,rhoe_F = 0;
                if (sleft >= 0 ){
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if(sleft < 0 && s_star >=0 ){
                    rho_F = rho_FL + sleft * (u_stat_L - rho_L);
                    rhou_F = rhou_FL + sleft * (u_stat_L * s_star  - rhou_L);
                    rhov_F = rhov_FL + sleft * (u_stat_L * v_L  - rhov_L);
                    rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L -  rhoe_L);
                }
                else if(s_star < 0 && sright > 0){
                    rho_F = rho_FR + sright * (u_stat_R - rho_R);
                    rhou_F = rhou_FR + sright * (u_stat_R * s_star  - rhou_R);
                    rhov_F = rhov_FR + sright * (u_stat_R * v_R  - rhov_R);
                    rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R -  rhoe_R);
    
                }
                else if (sright <= 0){
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }
                z[0][j][k] = rho_F; 
                z[1][j][k] = rhou_F; 
                z[2][j][k] = rhov_F;
                z[3][j][k] = rhoe_F; 
            }
        }   
    }

    else if (dir == 2){
        for (int j = GC; j < cols-GC; j++) {
            for (int k = GC-1; k < depth-GC; k++){
               //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量
                double rho_GL = rho_L * v_L;
                double rho_GR = rho_R * v_R;
                double rhou_GL = rho_L * u_L * v_L;
                double rhou_GR = rho_R * u_R * v_R;
                double rhov_GL = rho_L * v_L*v_L + p_L;
                double rhov_GR = rho_R * v_R*v_R + p_R;
                double rhoe_GL = rho_L * H_L * v_L;
                double rhoe_GR = rho_R * H_R * v_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = vbar - cbar;
                double sright = vbar + cbar;
            
                double s_star = (p_R - p_L + rho_L*v_L*(sleft - v_L) - rho_R*v_R*(sright - v_R))\
                                        /(rho_L*(sleft - v_L) - rho_R*(sright - v_R));
                                        
                double p_star = p_L + rho_L * (sleft-v_L) * (s_star-v_L);
                double v_stat_L = rho_L * (sleft-v_L)/(sleft-s_star);
                double v_stat_R = rho_R * (sright-v_R)/(sright-s_star);

                double fe_star_L = rhoe_L/rho_L + (s_star-v_L)*(s_star + p_L/(rho_L *(sleft-v_L)));
                double fe_star_R = rhoe_R/rho_R + (s_star-v_R)*(s_star + p_R/(rho_R *(sright-v_R)));
                //确定HLLC数值通量
                double rho_G = 0, rhou_G = 0, rhov_G = 0,rhoe_G = 0;
                if (sleft >= 0 ){
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if(sleft < 0 && s_star >=0 ){
                    rho_G = rho_GL + sleft * (v_stat_L - rho_L);
                    rhou_G = rhou_GL + sleft * (v_stat_L * u_L  - rhou_L);
                    rhov_G = rhov_GL + sleft * (v_stat_L * s_star   - rhov_L);
                    rhoe_G = rhoe_GL + sleft * (v_stat_L * fe_star_L -  rhoe_L);
                }
                else if(s_star < 0 && sright > 0){
                    rho_G = rho_GR + sright * (v_stat_R - rho_R);
                    rhou_G = rhou_GR + sright * (v_stat_R * u_R  - rhou_R);
                    rhov_G = rhov_GR + sright * (v_stat_R * s_star  - rhov_R);
                    rhoe_G = rhoe_GR + sright * (v_stat_R * fe_star_R -  rhoe_R);
    
                }
                else if (sright <= 0){
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }

                z[0][j][k] = rho_G; 
                z[1][j][k] = rhou_G; 
                z[2][j][k] = rhov_G;
                z[3][j][k] = rhoe_G; 

            }
        }
    }    
}

static inline void Roe_Flux(int dir, int rows, int cols, int depth, int GC, double (*x)[cols][depth], double (*y)[cols][depth], double (*z)[cols][depth]) {
    double epsilon = 1e-6;
    
    if (dir == 1) { // x方向通量
        for (int j = GC-1; j < cols-GC; j++) {
            for (int k = GC; k < depth-GC; k++) {
                // 读取左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];

                // 计算原始变量
                double u_L = rhou_L / rho_L;
                double u_R = rhou_R / rho_R;
                double v_L = rhov_L / rho_L;
                double v_R = rhov_R / rho_R;
                double p_L = (rhoe_L - 0.5 * rho_L * (u_L*u_L + v_L*v_L)) * (M_gamma-1);
                double p_R = (rhoe_R - 0.5 * rho_R * (u_R*u_R + v_R*v_R)) * (M_gamma-1);

                // 计算总焓H
                double H_L = 0.5 * (u_L*u_L + v_L*v_L) + (M_gamma/(M_gamma-1)) * (p_L/rho_L);
                double H_R = 0.5 * (u_R*u_R + v_R*v_R) + (M_gamma/(M_gamma-1)) * (p_R/rho_R);

                // 计算左右通量 (x方向)
                double rho_FL = rho_L * u_L;
                double rho_FR = rho_R * u_R;
                double rhou_FL = rho_L * u_L*u_L + p_L;
                double rhou_FR = rho_R * u_R*u_R + p_R;
                double rhov_FL = rho_L * u_L*v_L;
                double rhov_FR = rho_R * u_R*v_R;
                double rhoe_FL = rho_L * H_L * u_L;
                double rhoe_FR = rho_R * H_R * u_R;

                // 计算Roe平均
                double sqrt_rho_L = sqrt(rho_L);
                double sqrt_rho_R = sqrt(rho_R);
                double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                
                double rhobar = sqrt_rho_L * sqrt_rho_R;
                double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / sum_sqrt;
                double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) / sum_sqrt;
                double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / sum_sqrt;
                double q2bar = ubar*ubar + vbar*vbar;
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * q2bar));

                // 计算特征值
                double lambda[4];
                lambda[0] = ubar - cbar;
                lambda[1] = ubar;
                lambda[2] = ubar;
                lambda[3] = ubar + cbar;

                // 熵修正
                for (int i = 0; i < 4; i++) {
                    if (fabs(lambda[i]) < epsilon) {
                        lambda[i] = (lambda[i]*lambda[i] + epsilon*epsilon) / (2*epsilon);
                    } else {
                        lambda[i] = fabs(lambda[i]);
                    }
                }

                // 计算波强
                double drho = rho_R - rho_L;
                double du = u_R - u_L;
                double dv = v_R - v_L;
                double dp = p_R - p_L;

                double alpha1 = (dp - rhobar * cbar * du) / (2 * cbar*cbar);
                double alpha2 = drho - dp / (cbar*cbar);
                double alpha3 = rhobar * dv;
                double alpha4 = (dp + rhobar * cbar * du) / (2 * cbar*cbar);

                // 特征向量 (x方向)
                double K1[4] = {1, ubar-cbar, vbar, Hbar-ubar*cbar};
                double K2[4] = {1, ubar, vbar, 0.5*q2bar};
                double K3[4] = {0, 0, 1, vbar};
                double K4[4] = {1, ubar+cbar, vbar, Hbar+ubar*cbar};

                // Roe数值通量
                double rho_F = 0.5 * (rho_FL + rho_FR);
                double rhou_F = 0.5 * (rhou_FL + rhou_FR);
                double rhov_F = 0.5 * (rhov_FL + rhov_FR);
                double rhoe_F = 0.5 * (rhoe_FL + rhoe_FR);

                // 添加耗散项
                for (int i = 0; i < 4; i++) {
                    rho_F -= 0.5 * lambda[i] * alpha1 * K1[0];
                    rhou_F -= 0.5 * lambda[i] * alpha1 * K1[1];
                    rhov_F -= 0.5 * lambda[i] * alpha1 * K1[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha1 * K1[3];
                    
                    rho_F -= 0.5 * lambda[i] * alpha2 * K2[0];
                    rhou_F -= 0.5 * lambda[i] * alpha2 * K2[1];
                    rhov_F -= 0.5 * lambda[i] * alpha2 * K2[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha2 * K2[3];
                    
                    rho_F -= 0.5 * lambda[i] * alpha3 * K3[0];
                    rhou_F -= 0.5 * lambda[i] * alpha3 * K3[1];
                    rhov_F -= 0.5 * lambda[i] * alpha3 * K3[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha3 * K3[3];
                    
                    rho_F -= 0.5 * lambda[i] * alpha4 * K4[0];
                    rhou_F -= 0.5 * lambda[i] * alpha4 * K4[1];
                    rhov_F -= 0.5 * lambda[i] * alpha4 * K4[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha4 * K4[3];
                }

                z[0][j][k] = rho_F;
                z[1][j][k] = rhou_F;
                z[2][j][k] = rhov_F;
                z[3][j][k] = rhoe_F;
            }
        }
    }
    else if (dir == 2) { // y方向通量
        for (int j = GC; j < cols-GC; j++) {
            for (int k = GC-1; k < depth-GC; k++) {
                // 读取上下原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];

                // 计算原始变量
                double u_L = rhou_L / rho_L;
                double u_R = rhou_R / rho_R;
                double v_L = rhov_L / rho_L;
                double v_R = rhov_R / rho_R;
                double p_L = (rhoe_L - 0.5 * rho_L * (u_L*u_L + v_L*v_L)) * (M_gamma-1);
                double p_R = (rhoe_R - 0.5 * rho_R * (u_R*u_R + v_R*v_R)) * (M_gamma-1);

                // 计算总焓H
                double H_L = 0.5 * (u_L*u_L + v_L*v_L) + (M_gamma/(M_gamma-1)) * (p_L/rho_L);
                double H_R = 0.5 * (u_R*u_R + v_R*v_R) + (M_gamma/(M_gamma-1)) * (p_R/rho_R);

                // 计算左右通量 (y方向)
                double rho_GL = rho_L * v_L;
                double rho_GR = rho_R * v_R;
                double rhou_GL = rho_L * u_L*v_L;
                double rhou_GR = rho_R * u_R*v_R;
                double rhov_GL = rho_L * v_L*v_L + p_L;
                double rhov_GR = rho_R * v_R*v_R + p_R;
                double rhoe_GL = rho_L * H_L * v_L;
                double rhoe_GR = rho_R * H_R * v_R;

                // 计算Roe平均
                double sqrt_rho_L = sqrt(rho_L);
                double sqrt_rho_R = sqrt(rho_R);
                double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                
                double rhobar = sqrt_rho_L * sqrt_rho_R;
                double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / sum_sqrt;
                double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) / sum_sqrt;
                double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / sum_sqrt;
                double q2bar = ubar*ubar + vbar*vbar;
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * q2bar));

                // 计算特征值 (y方向)
                double lambda[4];
                lambda[0] = vbar - cbar;
                lambda[1] = vbar;
                lambda[2] = vbar;
                lambda[3] = vbar + cbar;

                // 熵修正
                for (int i = 0; i < 4; i++) {
                    if (fabs(lambda[i]) < epsilon) {
                        lambda[i] = (lambda[i]*lambda[i] + epsilon*epsilon) / (2*epsilon);
                    } else {
                        lambda[i] = fabs(lambda[i]);
                    }
                }

                // 计算波强
                double drho = rho_R - rho_L;
                double du = u_R - u_L;
                double dv = v_R - v_L;
                double dp = p_R - p_L;

                double alpha1 = (dp - rhobar * cbar * dv) / (2 * cbar*cbar);
                double alpha2 = drho - dp / (cbar*cbar);
                double alpha3 = rhobar * du;
                double alpha4 = (dp + rhobar * cbar * dv) / (2 * cbar*cbar);

                // 特征向量 (y方向)
                double K1[4] = {1, ubar, vbar-cbar, Hbar-vbar*cbar};
                double K2[4] = {1, ubar, vbar, 0.5*q2bar};
                double K3[4] = {0, 1, 0, ubar};
                double K4[4] = {1, ubar, vbar+cbar, Hbar+vbar*cbar};

                // Roe数值通量
                double rho_G= 0.5 * (rho_GL + rho_GR);
                double rhou_G = 0.5 * (rhou_GL + rhou_GR);
                double rhov_G = 0.5 * (rhov_GL + rhov_GR);
                double rhoe_G = 0.5 * (rhoe_GL + rhoe_GR);

                // 添加耗散项
                for (int i = 0; i < 4; i++) {
                    rho_G-= 0.5 * lambda[i] * alpha1 * K1[0];
                    rhou_G -= 0.5 * lambda[i] * alpha1 * K1[1];
                    rhov_G -= 0.5 * lambda[i] * alpha1 * K1[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha1 * K1[3];
                    
                    rho_G-= 0.5 * lambda[i] * alpha2 * K2[0];
                    rhou_G -= 0.5 * lambda[i] * alpha2 * K2[1];
                    rhov_G -= 0.5 * lambda[i] * alpha2 * K2[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha2 * K2[3];
                    
                    rho_G-= 0.5 * lambda[i] * alpha3 * K3[0];
                    rhou_G -= 0.5 * lambda[i] * alpha3 * K3[1];
                    rhov_G -= 0.5 * lambda[i] * alpha3 * K3[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha3 * K3[3];
                    
                    rho_G-= 0.5 * lambda[i] * alpha4 * K4[0];
                    rhou_G -= 0.5 * lambda[i] * alpha4 * K4[1];
                    rhov_G -= 0.5 * lambda[i] * alpha4 * K4[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha4 * K4[3];
                }

                z[0][j][k] = rho_G;
                z[1][j][k] = rhou_G;
                z[2][j][k] = rhov_G;
                z[3][j][k] = rhoe_G;
            }
        }
    }
}


static inline void ER_Flux(int rows, int cols, int GC, double (*x)[cols], double (*y)[cols] ,double (*z)[cols]) {
    
    for (int j = GC-1; j <= cols-GC; j++) {
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
//Flux计算方法


static inline void HLLHC_Flux(int dir, int rows, int cols, int depth, int GC, double (*x)[cols][depth], double (*y)[cols][depth] ,double (*z)[cols][depth]) {


    if (dir == 1){
       for (int j = GC-1; j < cols-GC; j++) {
            for (int k = GC; k < depth-GC; k++){
                //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);
                
                //计算声速
                double a_L = sqrt(M_gamma * p_L / rho_L);
                double a_R = sqrt(M_gamma * p_R / rho_R);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量

                double rho_FL = rho_L * u_L;
                double rho_FR = rho_R * u_R;
                double rhou_FL = rho_L * u_L * u_L + p_L;
                double rhou_FR = rho_R * u_R * u_R + p_R;
                double rhov_FL = rho_L * v_L * u_L;
                double rhov_FR = rho_R * v_R * u_R;
                double rhoe_FL = rho_L * H_L * u_L;
                double rhoe_FR = rho_R * H_R * u_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = ubar - cbar;
                double sright = ubar + cbar;

                //确定HLL数值通量
                double rho_F = 0, rhou_F = 0, rhov_F = 0,rhoe_F = 0;
                if (sleft >= 0 ){
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if(sleft < 0 && sright > 0){
                    rho_F = (sright*rho_FL - sleft*rho_FR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
                    rhou_F = (sright*rhou_FL - sleft*rhou_FR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
                    rhov_F = (sright*rhov_FL - sleft*rhov_FR + sleft*sright * (rhov_R - rhov_L))/(sright-sleft);
                    rhoe_F = (sright*rhoe_FL - sleft*rhoe_FR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
                }
                else if (sright <= 0){
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }

                z[0][j][k] = rho_F; 
                z[1][j][k] = rhou_F; 
                z[2][j][k] = rhov_F;
                z[3][j][k] = rhoe_F; 
            }
        }   
    }

    else if (dir == 2){
        for (int j = GC; j < cols-GC; j++) {
            for (int k = GC-1; k < depth-GC; k++){
               //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量

                double rho_GL = rho_L * v_L;
                double rho_GR = rho_R * v_R;
                double rhou_GL = rho_L * u_L * v_L;
                double rhou_GR = rho_R * u_R * v_R;
                double rhov_GL = rho_L * v_L*v_L + p_L;
                double rhov_GR = rho_R * v_R*v_R + p_R;
                double rhoe_GL = rho_L * H_L * v_L;
                double rhoe_GR = rho_R * H_R * v_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = vbar - cbar;
                double sright = vbar + cbar;

                //确定HLL数值通量
                double rho_G = 0, rhou_G = 0, rhov_G = 0,rhoe_G = 0;
                if (sleft >= 0 ){
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if(sleft < 0 && sright >0){
                    rho_G = (sright*rho_GL - sleft*rho_GR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
                    rhou_G = (sright*rhou_GL - sleft*rhou_GR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
                    rhov_G = (sright*rhov_GL - sleft*rhov_GR + sleft*sright * (rhov_R - rhov_L))/(sright-sleft);
                    rhoe_G = (sright*rhoe_GL - sleft*rhoe_GR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
                }

                else if (sright <= 0){
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }

                z[0][j][k] = rho_G; 
                z[1][j][k] = rhou_G; 
                z[2][j][k] = rhov_G;
                z[3][j][k] = rhoe_G; 


            }
        }
    }
    
}



static inline void HLLCHC_Flux(int dir, int rows, int cols, int depth, int GC, double (*x)[cols][depth], double (*y)[cols][depth] ,double (*z)[cols][depth]) {


    if (dir == 1){
       for (int j = GC-1; j < cols-GC; j++) {
            for (int k = GC; k < depth-GC; k++){
                //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);
                
                //计算声速
                double a_L = sqrt(M_gamma * p_L / rho_L);
                double a_R = sqrt(M_gamma * p_R / rho_R);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量

                double rho_FL = rho_L * u_L;
                double rho_FR = rho_R * u_R;
                double rhou_FL = rho_L * u_L * u_L + p_L;
                double rhou_FR = rho_R * u_R * u_R + p_R;
                double rhov_FL = rho_L * v_L * u_L;
                double rhov_FR = rho_R * v_R * u_R;
                double rhoe_FL = rho_L * H_L * u_L;
                double rhoe_FR = rho_R * H_R * u_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = ubar - cbar;
                double sright = ubar + cbar;
            
                double s_star = (p_R - p_L + rho_L*u_L*(sleft - u_L) - rho_R*u_R*(sright - u_R))\
                                        /(rho_L*(sleft - u_L) - rho_R*(sright - u_R));
                                        
                double p_star = p_L + rho_L * (sleft-u_L) * (s_star-u_L);
                double u_stat_L = rho_L * (sleft-u_L)/(sleft-s_star);
                double u_stat_R = rho_R * (sright-u_R)/(sright-s_star);

                double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L)));
                double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R)));
                //确定HLLC数值通量
                double rho_F = 0, rhou_F = 0, rhov_F = 0,rhoe_F = 0;
                if (sleft >= 0 ){
                    rho_F = rho_FL;
                    rhou_F = rhou_FL;
                    rhov_F = rhov_FL;
                    rhoe_F = rhoe_FL;
                }
                else if(sleft < 0 && s_star >=0 ){
                    rho_F = rho_FL + sleft * (u_stat_L - rho_L);
                    rhou_F = rhou_FL + sleft * (u_stat_L * s_star  - rhou_L);
                    rhov_F = rhov_FL + sleft * (u_stat_L * v_L  - rhov_L);
                    rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L -  rhoe_L);
                }
                else if(s_star < 0 && sright > 0){
                    rho_F = rho_FR + sright * (u_stat_R - rho_R);
                    rhou_F = rhou_FR + sright * (u_stat_R * s_star  - rhou_R);
                    rhov_F = rhov_FR + sright * (u_stat_R * v_R  - rhov_R);
                    rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R -  rhoe_R);
    
                }
                else if (sright <= 0){
                    rho_F = rho_FR;
                    rhou_F = rhou_FR;
                    rhov_F = rhov_FR;
                    rhoe_F = rhoe_FR;
                }
                z[0][j][k] = rho_F; 
                z[1][j][k] = rhou_F; 
                z[2][j][k] = rhov_F;
                z[3][j][k] = rhoe_F; 
            }
        }   
    }

    else if (dir == 2){
        for (int j = GC; j < cols-GC; j++) {
            for (int k = GC-1; k < depth-GC; k++){
               //读取已知的左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];


                double u_L = x[1][j][k]/x[0][j][k];
                double u_R = y[1][j][k]/y[0][j][k];
                double v_L = x[2][j][k]/x[0][j][k];
                double v_R = y[2][j][k]/y[0][j][k];
                double p_L = (x[3][j][k] - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2)))*(M_gamma-1);
                double p_R = (y[3][j][k] - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2)))*(M_gamma-1);

                //计算总焓H
                double H_L = 0.5 * pow(u_L,2) +  0.5 * pow(v_L,2) + (M_gamma/(M_gamma-1) ) * (p_L/rho_L);
                double H_R = 0.5 * pow(u_R,2) +  0.5 * pow(v_R,2) + (M_gamma/(M_gamma-1) ) * (p_R/rho_R);

                //计算左右守恒变量和通量
                double rho_GL = rho_L * v_L;
                double rho_GR = rho_R * v_R;
                double rhou_GL = rho_L * u_L * v_L;
                double rhou_GR = rho_R * u_R * v_R;
                double rhov_GL = rho_L * v_L*v_L + p_L;
                double rhov_GR = rho_R * v_R*v_R + p_R;
                double rhoe_GL = rho_L * H_L * v_L;
                double rhoe_GR = rho_R * H_R * v_R;
            
                //计算Roe平均
                double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
                double vbar = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R)/(sqrt(rho_L)+sqrt(rho_R));
                double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

                //利用Roe平均的变量计算近似波速
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * pow(ubar,2) - 0.5 * pow(vbar,2)));
                double sleft = vbar - cbar;
                double sright = vbar + cbar;
            
                double s_star = (p_R - p_L + rho_L*v_L*(sleft - v_L) - rho_R*v_R*(sright - v_R))\
                                        /(rho_L*(sleft - v_L) - rho_R*(sright - v_R));
                                        
                double p_star = p_L + rho_L * (sleft-v_L) * (s_star-v_L);
                double v_stat_L = rho_L * (sleft-v_L)/(sleft-s_star);
                double v_stat_R = rho_R * (sright-v_R)/(sright-s_star);

                double fe_star_L = rhoe_L/rho_L + (s_star-v_L)*(s_star + p_L/(rho_L *(sleft-v_L)));
                double fe_star_R = rhoe_R/rho_R + (s_star-v_R)*(s_star + p_R/(rho_R *(sright-v_R)));
                //确定HLLC数值通量
                double rho_G = 0, rhou_G = 0, rhov_G = 0,rhoe_G = 0;
                if (sleft >= 0 ){
                    rho_G = rho_GL;
                    rhou_G = rhou_GL;
                    rhov_G = rhov_GL;
                    rhoe_G = rhoe_GL;
                }
                else if(sleft < 0 && s_star >=0 ){
                    rho_G = rho_GL + sleft * (v_stat_L - rho_L);
                    rhou_G = rhou_GL + sleft * (v_stat_L * u_L  - rhou_L);
                    rhov_G = rhov_GL + sleft * (v_stat_L * s_star   - rhov_L);
                    rhoe_G = rhoe_GL + sleft * (v_stat_L * fe_star_L -  rhoe_L);
                }
                else if(s_star < 0 && sright > 0){
                    rho_G = rho_GR + sright * (v_stat_R - rho_R);
                    rhou_G = rhou_GR + sright * (v_stat_R * u_R  - rhou_R);
                    rhov_G = rhov_GR + sright * (v_stat_R * s_star  - rhov_R);
                    rhoe_G = rhoe_GR + sright * (v_stat_R * fe_star_R -  rhoe_R);
    
                }
                else if (sright <= 0){
                    rho_G = rho_GR;
                    rhou_G = rhou_GR;
                    rhov_G = rhov_GR;
                    rhoe_G = rhoe_GR;
                }

                z[0][j][k] = rho_G; 
                z[1][j][k] = rhou_G; 
                z[2][j][k] = rhov_G;
                z[3][j][k] = rhoe_G; 

            }
        }
    }    
}

static inline void RoeHC_Flux(int dir, int rows, int cols, int depth, int GC, double (*x)[cols][depth], double (*y)[cols][depth], double (*z)[cols][depth]) {
    double epsilon = 1e-6;
    
    if (dir == 1) { // x方向通量
        for (int j = GC-1; j < cols-GC; j++) {
            for (int k = GC; k < depth-GC; k++) {
                // 读取左右原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];

                // 计算原始变量
                double u_L = rhou_L / rho_L;
                double u_R = rhou_R / rho_R;
                double v_L = rhov_L / rho_L;
                double v_R = rhov_R / rho_R;
                double p_L = (rhoe_L - 0.5 * rho_L * (u_L*u_L + v_L*v_L)) * (M_gamma-1);
                double p_R = (rhoe_R - 0.5 * rho_R * (u_R*u_R + v_R*v_R)) * (M_gamma-1);

                // 计算总焓H
                double H_L = 0.5 * (u_L*u_L + v_L*v_L) + (M_gamma/(M_gamma-1)) * (p_L/rho_L);
                double H_R = 0.5 * (u_R*u_R + v_R*v_R) + (M_gamma/(M_gamma-1)) * (p_R/rho_R);

                // 计算左右通量 (x方向)
                double rho_FL = rho_L * u_L;
                double rho_FR = rho_R * u_R;
                double rhou_FL = rho_L * u_L*u_L + p_L;
                double rhou_FR = rho_R * u_R*u_R + p_R;
                double rhov_FL = rho_L * u_L*v_L;
                double rhov_FR = rho_R * u_R*v_R;
                double rhoe_FL = rho_L * H_L * u_L;
                double rhoe_FR = rho_R * H_R * u_R;

                // 计算Roe平均
                double sqrt_rho_L = sqrt(rho_L);
                double sqrt_rho_R = sqrt(rho_R);
                double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                
                double rhobar = sqrt_rho_L * sqrt_rho_R;
                double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / sum_sqrt;
                double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) / sum_sqrt;
                double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / sum_sqrt;
                double q2bar = ubar*ubar + vbar*vbar;
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * q2bar));

                // 计算特征值
                double lambda[4];
                lambda[0] = ubar - cbar;
                lambda[1] = ubar;
                lambda[2] = ubar;
                lambda[3] = ubar + cbar;

                // 熵修正
                for (int i = 0; i < 4; i++) {
                    if (fabs(lambda[i]) < epsilon) {
                        lambda[i] = (lambda[i]*lambda[i] + epsilon*epsilon) / (2*epsilon);
                    } else {
                        lambda[i] = fabs(lambda[i]);
                    }
                }

                // 计算波强
                double drho = rho_R - rho_L;
                double du = u_R - u_L;
                double dv = v_R - v_L;
                double dp = p_R - p_L;

                double alpha1 = (dp - rhobar * cbar * du) / (2 * cbar*cbar);
                double alpha2 = drho - dp / (cbar*cbar);
                double alpha3 = rhobar * dv;
                double alpha4 = (dp + rhobar * cbar * du) / (2 * cbar*cbar);

                // 特征向量 (x方向)
                double K1[4] = {1, ubar-cbar, vbar, Hbar-ubar*cbar};
                double K2[4] = {1, ubar, vbar, 0.5*q2bar};
                double K3[4] = {0, 0, 1, vbar};
                double K4[4] = {1, ubar+cbar, vbar, Hbar+ubar*cbar};

                // Roe数值通量
                double rho_F = 0.5 * (rho_FL + rho_FR);
                double rhou_F = 0.5 * (rhou_FL + rhou_FR);
                double rhov_F = 0.5 * (rhov_FL + rhov_FR);
                double rhoe_F = 0.5 * (rhoe_FL + rhoe_FR);

                // 添加耗散项
                for (int i = 0; i < 4; i++) {
                    rho_F -= 0.5 * lambda[i] * alpha1 * K1[0];
                    rhou_F -= 0.5 * lambda[i] * alpha1 * K1[1];
                    rhov_F -= 0.5 * lambda[i] * alpha1 * K1[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha1 * K1[3];
                    
                    rho_F -= 0.5 * lambda[i] * alpha2 * K2[0];
                    rhou_F -= 0.5 * lambda[i] * alpha2 * K2[1];
                    rhov_F -= 0.5 * lambda[i] * alpha2 * K2[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha2 * K2[3];
                    
                    rho_F -= 0.5 * lambda[i] * alpha3 * K3[0];
                    rhou_F -= 0.5 * lambda[i] * alpha3 * K3[1];
                    rhov_F -= 0.5 * lambda[i] * alpha3 * K3[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha3 * K3[3];
                    
                    rho_F -= 0.5 * lambda[i] * alpha4 * K4[0];
                    rhou_F -= 0.5 * lambda[i] * alpha4 * K4[1];
                    rhov_F -= 0.5 * lambda[i] * alpha4 * K4[2];
                    rhoe_F -= 0.5 * lambda[i] * alpha4 * K4[3];
                }

                z[0][j][k] = rho_F;
                z[1][j][k] = rhou_F;
                z[2][j][k] = rhov_F;
                z[3][j][k] = rhoe_F;
            }
        }
    }
    else if (dir == 2) { // y方向通量
        for (int j = GC; j < cols-GC; j++) {
            for (int k = GC-1; k < depth-GC; k++) {
                // 读取上下原始变量
                double rho_L = x[0][j][k];
                double rho_R = y[0][j][k];
                double rhou_L = x[1][j][k];
                double rhou_R = y[1][j][k];
                double rhov_L = x[2][j][k];
                double rhov_R = y[2][j][k];
                double rhoe_L = x[3][j][k];
                double rhoe_R = y[3][j][k];

                // 计算原始变量
                double u_L = rhou_L / rho_L;
                double u_R = rhou_R / rho_R;
                double v_L = rhov_L / rho_L;
                double v_R = rhov_R / rho_R;
                double p_L = (rhoe_L - 0.5 * rho_L * (u_L*u_L + v_L*v_L)) * (M_gamma-1);
                double p_R = (rhoe_R - 0.5 * rho_R * (u_R*u_R + v_R*v_R)) * (M_gamma-1);

                // 计算总焓H
                double H_L = 0.5 * (u_L*u_L + v_L*v_L) + (M_gamma/(M_gamma-1)) * (p_L/rho_L);
                double H_R = 0.5 * (u_R*u_R + v_R*v_R) + (M_gamma/(M_gamma-1)) * (p_R/rho_R);

                // 计算左右通量 (y方向)
                double rho_GL = rho_L * v_L;
                double rho_GR = rho_R * v_R;
                double rhou_GL = rho_L * u_L*v_L;
                double rhou_GR = rho_R * u_R*v_R;
                double rhov_GL = rho_L * v_L*v_L + p_L;
                double rhov_GR = rho_R * v_R*v_R + p_R;
                double rhoe_GL = rho_L * H_L * v_L;
                double rhoe_GR = rho_R * H_R * v_R;

                // 计算Roe平均
                double sqrt_rho_L = sqrt(rho_L);
                double sqrt_rho_R = sqrt(rho_R);
                double sum_sqrt = sqrt_rho_L + sqrt_rho_R;
                
                double rhobar = sqrt_rho_L * sqrt_rho_R;
                double ubar = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / sum_sqrt;
                double vbar = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) / sum_sqrt;
                double Hbar = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / sum_sqrt;
                double q2bar = ubar*ubar + vbar*vbar;
                double cbar = sqrt((M_gamma-1) * (Hbar - 0.5 * q2bar));

                // 计算特征值 (y方向)
                double lambda[4];
                lambda[0] = vbar - cbar;
                lambda[1] = vbar;
                lambda[2] = vbar;
                lambda[3] = vbar + cbar;

                // 熵修正
                for (int i = 0; i < 4; i++) {
                    if (fabs(lambda[i]) < epsilon) {
                        lambda[i] = (lambda[i]*lambda[i] + epsilon*epsilon) / (2*epsilon);
                    } else {
                        lambda[i] = fabs(lambda[i]);
                    }
                }

                // 计算波强
                double drho = rho_R - rho_L;
                double du = u_R - u_L;
                double dv = v_R - v_L;
                double dp = p_R - p_L;

                double alpha1 = (dp - rhobar * cbar * dv) / (2 * cbar*cbar);
                double alpha2 = drho - dp / (cbar*cbar);
                double alpha3 = rhobar * du;
                double alpha4 = (dp + rhobar * cbar * dv) / (2 * cbar*cbar);

                // 特征向量 (y方向)
                double K1[4] = {1, ubar, vbar-cbar, Hbar-vbar*cbar};
                double K2[4] = {1, ubar, vbar, 0.5*q2bar};
                double K3[4] = {0, 1, 0, ubar};
                double K4[4] = {1, ubar, vbar+cbar, Hbar+vbar*cbar};

                // Roe数值通量
                double rho_G= 0.5 * (rho_GL + rho_GR);
                double rhou_G = 0.5 * (rhou_GL + rhou_GR);
                double rhov_G = 0.5 * (rhov_GL + rhov_GR);
                double rhoe_G = 0.5 * (rhoe_GL + rhoe_GR);

                // 添加耗散项
                for (int i = 0; i < 4; i++) {
                    rho_G-= 0.5 * lambda[i] * alpha1 * K1[0];
                    rhou_G -= 0.5 * lambda[i] * alpha1 * K1[1];
                    rhov_G -= 0.5 * lambda[i] * alpha1 * K1[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha1 * K1[3];
                    
                    rho_G-= 0.5 * lambda[i] * alpha2 * K2[0];
                    rhou_G -= 0.5 * lambda[i] * alpha2 * K2[1];
                    rhov_G -= 0.5 * lambda[i] * alpha2 * K2[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha2 * K2[3];
                    
                    rho_G-= 0.5 * lambda[i] * alpha3 * K3[0];
                    rhou_G -= 0.5 * lambda[i] * alpha3 * K3[1];
                    rhov_G -= 0.5 * lambda[i] * alpha3 * K3[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha3 * K3[3];
                    
                    rho_G-= 0.5 * lambda[i] * alpha4 * K4[0];
                    rhou_G -= 0.5 * lambda[i] * alpha4 * K4[1];
                    rhov_G -= 0.5 * lambda[i] * alpha4 * K4[2];
                    rhoe_G -= 0.5 * lambda[i] * alpha4 * K4[3];
                }

                z[0][j][k] = rho_G;
                z[1][j][k] = rhou_G;
                z[2][j][k] = rhov_G;
                z[3][j][k] = rhoe_G;
            }
        }
    }
}


static inline void HLL_Flux_HeatConduction(int rows, int cols, int GC, double (*x)[cols], double (*y)[cols] ,double (*z)[cols], double u_point) {
    
    for (int j = GC-1; j <= cols-GC; j++) {
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



static inline void HLLC_Flux_HeatConduction(int rows, int cols, int GC, double (*x)[cols], double (*y)[cols] ,double (*z)[cols], double u_point) {
    
    for (int j = GC-1; j <= cols-GC; j++) {
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
        double ma = u_L/cbar;
        Splus = fabs(ubar) + cbar;

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



static inline void Roe_Flux_HeatConduction(int rows, int cols, int GC, double (*x)[cols], double (*y)[cols] ,double (*z)[cols], double u_point) {
    double epsilon = 1e-6;
    
    for (int j = GC-1; j <= cols-GC; j++) {
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


static inline void ER_Flux_PlusHeat(int rows, int cols, int GC, double (*x)[cols], double (*y)[cols] ,double (*z)[cols]) {
    
    for (int j = GC-1; j <= cols-GC; j++) {
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
static inline void RS_Marquina(int rows, int cols,double (*y)[cols],double (*z)[cols], double dt,double dx) {
    double conserl[3][cols],conserr[3][cols];
    double fluxl[3][cols],fluxr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double pri[rows][cols];
    double slope[3][cols],a[3][cols],b[3][cols];
    
    double eigen_l[3][3][cols],eigen_r[3][3][cols];
    double w_l[rows][cols], w_r[rows][cols];
    double phi_fl[rows][cols], phi_fr[rows][cols];
    double phi_fp[rows][cols], phi_fm[rows][cols];
    double flux[3][cols];
    double lamda[3][cols];
    double alpha[rows][cols];

    //初始化数组
     for (int k = 0; k < rows; k++){
        for (int j = 0; j < cols; j++){
            w_l[k][j] = 0.0;
            w_r[k][j] = 0.0;
            phi_fl[k][j] = 0.0;
            phi_fr[k][j] = 0.0;
            phi_fp[k][j] = 0.0;
            phi_fm[k][j] = 0.0;
            flux[k][j] = 0.0;
        }
        
    }

    Con_to_Pri_1D(3,cols,pri,y);
    //计算特征矩阵
	for(int j=0;  j < cols; j++) {
			
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
//            Reconstruction_Godunov(rows,cols,3,y,conserl,conserr,dx);
            break;
        case 1:
//            TVD_Reconstruction(rows,cols,3,y,conserl,conserr,dx);
            break;
        case 3:
//            WENO3_Reconstruction(rows,cols,3,y,conserl,conserr);
            break;
        case 5:
//            WENO5_Reconstruction_C(rows,cols,3,y,conserl,conserr);
            break;
        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }
                              

    
    //执行计算Marquina flux
    Con_to_Pri_1D(3,cols,prileft,conserl);
    Con_to_Pri_1D(3,cols,priright,conserr);

    initEulerflux1D(rows, cols, conserl, fluxl);
    initEulerflux1D(rows, cols, conserr, fluxr);

    // 投影到特征空间
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_l[k][j] += conserl[m][j]*eigen_l[k][m][j];
                phi_fl[k][j] += fluxl[m][j]*eigen_l[k][m][j];
            }
        }
        
    }
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_r[k][j] += conserr[m][j]*eigen_l[k][m][j+1];
                phi_fr[k][j] += fluxr[m][j]*eigen_l[k][m][j+1];
            }
        }
        
    }

    //计算通量分量
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
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
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
               flux[k][j] += phi_fp[m][j]*eigen_r[k][m][j] + phi_fm[m][j] * eigen_r[k][m][j+1];
            }
        }  
    }

    for (int i = 0; i < rows; i++){
        for (int j = 1; j < cols-1; j++){
            z[i][j] = flux[i][j];
        }
    }

}


#endif  