#ifndef WENO_FLUX_H
#define WENO_FLUX_H

#include <stdio.h>
#include <math.h>

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


static inline void Con_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols],double gamma){
    int j;
    for (j = 0; j <= cols-1; j++) {
        x[0][j] = y[0][j];
        x[1][j] = y[1][j]/y[0][j];
        x[2][j] = (gamma-1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j]);
    }
}

static inline void Pri_to_Con_1D(int rows, int cols, double (*x)[cols], double (*y)[cols],double gamma){
    int j;
    for (j = 0; j <= cols-1; j++) {
        y[0][j] = x[0][j];
        y[1][j] = x[0][j] * x[1][j];
        y[2][j] = 0.5 * x[0][j] * pow(x[1][j],2) + x[2][j]/(gamma-1);
    }
}

static inline void Pri_to_S_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j < cols-1; j++) {
        y[0][j] = x[2][j] / pow(x[0][j], 1.4);
        y[1][j] = x[1][j];
        y[2][j] = x[2][j];
    }
}

static inline void S_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j < cols-1; j++) {
        x[0][j] = pow(y[2][j]/y[0][j],1/1.4);
        x[1][j] = y[1][j];
        x[2][j] = y[2][j];
    }
}

static inline void ConS_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j < cols-1; j++) {
        x[0][j] = y[0][j];
        x[1][j] = y[1][j] / y[0][j];
        x[2][j] = y[2][j] * pow(y[0][j],0.4);
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
static inline void Rusanov_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {


        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数
        //计算声速
        double a_L = sqrt(gamma*p_L/rho_L);
        double a_R = sqrt(gamma*p_R/rho_R);

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;

        //计算一个最大波速
        double Splus = max_of_two(fabs(u_L)  + a_L, fabs(u_R) + a_R);
    
        //Rusanov数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        rho_F = 0.5*(rho_FL + rho_FR) - 0.5*Splus*(rho_R - rho_L);
        rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*Splus*(rhou_R - rhou_L);
        rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*Splus*(rhoe_R - rhoe_L);
        
        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}


static inline void Lax_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    double S_plus = 0;

    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        //计算声速
        double a_L = sqrt(gamma*p_L/rho_L);
        double a_R = sqrt(gamma*p_R/rho_R);
        //计算全局最大波速
        S_plus = max_of_three(fabs(u_L)+a_L, fabs(u_R)+a_R,S_plus);
    
    }   
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
     

        //计算需要使用的参数
        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;

        //Rusanov数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        rho_F = 0.5*(rho_FL + rho_FR) - 0.5*S_plus*(rho_R - rho_L);
        rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*S_plus*(rhou_R - rhou_L);
        rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*S_plus*(rhoe_R - rhoe_L);
        
        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLL_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        double a_L = sqrt(gamma * p_L / rho_L);
        double a_R = sqrt(gamma * p_R / rho_R);

        //计算需要使用的参数
        //计算声速

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);


        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double Splus = max_of_two(fabs(u_L)  + a_L, fabs(u_R) + a_R);

        //确定HLL数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        if (sleft >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && sright >0){
            rho_F = (sright*rho_FL - sleft*rho_FR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
            rhou_F = (sright*rhou_FL - sleft*rhou_FR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
            rhoe_F = (sright*rhoe_FL - sleft*rhoe_FR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
        }
        else if (sright <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLLC_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j <= cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double s_star = (p_R - p_L + rho_L*u_L*(sleft - u_L) - rho_R*u_R*(sright - u_R))\
                                 /(rho_L*(sleft - u_L) - rho_R*(sright - u_R));
                                
        double p_star = p_L + rho_L * (sleft-u_L) * (s_star-u_L);
        double u_stat_L = rho_L * (sleft-u_L)/(sleft-s_star);
        double u_stat_R = rho_R * (sright-u_R)/(sright-s_star);
        double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L)));
        double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R)));

     
        //HLLC数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        if (sleft >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && s_star >=0 ){
            rho_F = rho_FL + sleft * (u_stat_L - rho_L);
            rhou_F = rhou_FL + sleft *(u_stat_L * s_star  - rhou_L);
            rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L -  rhoe_L);

        }
        else if(s_star < 0 && sright > 0){
            rho_F = rho_FR + sright * (u_stat_R - rho_R);
            rhou_F = rhou_FR + sright *(u_stat_R * s_star  - rhou_R);
            rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R -  rhoe_R);
        }
        else if (sright <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
        
    }   
}



static inline void Roe_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    double epsilon = 1e-6;
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

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
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));


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

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 

    }   
}



static inline void StegerWarming_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_P = x[0][j];
        double rho_N = y[0][j];
        double u_P = x[1][j];
        double u_N = y[1][j];
        double p_P = x[2][j];
        double p_N = y[2][j];
        //计算声速
        double a_P = sqrt(gamma * p_P / rho_P);
        double a_N = sqrt(gamma * p_N / rho_N);
        //计算总焓H
        double H_P = 0.5 * pow(u_P,2) + (gamma/(gamma-1)) * (p_P/rho_P);
        double H_N = 0.5 * pow(u_N,2) + (gamma/(gamma-1)) * (p_N/rho_N);

        //计算Roe平均
        double rhobar = sqrt(rho_P * rho_N);
        double ubar = (sqrt(rho_P) * u_P + sqrt(rho_N) * u_N)/(sqrt(rho_P)+sqrt(rho_N));
        double Hbar = (sqrt(rho_P) * H_P + sqrt(rho_N) * H_N)/(sqrt(rho_P)+sqrt(rho_N));
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));

        //利用Roe平均计算近似波速&&&&同时分裂为正负波速度
        double lambdap1 = u_P - a_P;
        double lambdap2 = u_P;
        double lambdap3 = u_P + a_P;
        double lambdan1 = u_N - a_N;
        double lambdan2 = u_N;
        double lambdan3 = u_N + a_N;
        double lambda1_p = 0.5*(lambdap1  + fabs(lambdap1));
        double lambda2_p = 0.5*(lambdap2  + fabs(lambdap2));
        double lambda3_p = 0.5*(lambdap3  + fabs(lambdap3));
        double lambda1_n = 0.5*(lambdan1  - fabs(lambdan1));
        double lambda2_n = 0.5*(lambdan2  - fabs(lambdan2));
        double lambda3_n = 0.5*(lambdan3  - fabs(lambdan3));
        /*double lambda1 = ubar - cbar;
        double lambda2 = ubar;
        double lambda3 = ubar + cbar;
        double lambda1_p = 0.5*(lambda1  + fabs(lambda1));
        double lambda2_p = 0.5*(lambda2  + fabs(lambda2));
        double lambda3_p = 0.5*(lambda3  + fabs(lambda3));
        double lambda1_n = 0.5*(lambda1  - fabs(lambda1));
        double lambda2_n = 0.5*(lambda2  - fabs(lambda2));
        double lambda3_n = 0.5*(lambda3  - fabs(lambda3));*/

        //计算中间变量
        double F_rhou_PP =  u_P - a_P;
        double F_rhou_PN =  u_P + a_P;
        double F_rhou_NP =  u_N - a_N;
        double F_rhou_NN =  u_N + a_N;
        double F_rhoe_PP =  H_P - u_P *a_P;
        double F_rhoe_PN =  H_P + u_P *a_P;
        double F_rhoe_NP =  H_N - u_N *a_N;
        double F_rhoe_NN =  H_N + u_N *a_N;

        //计算分裂之后的正负通量
        double rho_FP = (0.5*rho_P/gamma) * (lambda1_p + 2.0*(gamma-1)*lambda2_p + lambda3_p);
        double rho_FN = (0.5*rho_N/gamma) * (lambda1_n + 2.0*(gamma-1)*lambda2_n + lambda3_n);
        double rhou_FP = (0.5*rho_P/gamma) * (F_rhou_PP*lambda1_p + 2.0*(gamma-1)*ubar*lambda2_p + F_rhou_PN* lambda3_p);
        double rhou_FN = (0.5*rho_N/gamma) * (F_rhou_NP*lambda1_n + 2.0*(gamma-1)*ubar*lambda2_n + F_rhou_NN* lambda3_n);
        double rhoe_FP = (0.5*rho_P/gamma) * (F_rhoe_PP*lambda1_p + (gamma-1)*pow(u_P,2)*lambda2_p + F_rhoe_PN * lambda3_p);
        double rhoe_FN = (0.5*rho_N/gamma) * (F_rhoe_NP*lambda1_n + (gamma-1)*pow(u_N,2)*lambda2_n + F_rhoe_NN * lambda3_n);

        z[0][j] = rho_FP + rho_FN; 
        z[1][j] = rhou_FP + rhou_FN; 
        z[2][j] = rhoe_FP + rhoe_FN; 
    }   
}


static inline void VanLeer_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_P = x[0][j];
        double rho_N = y[0][j];
        double u_P = x[1][j];
        double u_N = y[1][j];
        double p_P = x[2][j];
        double p_N = y[2][j];
        //计算声速
        double a_P = sqrt(gamma * p_P / rho_P);
        double a_N = sqrt(gamma * p_N / rho_N);

        //计算当地Ma
        double Ma_P = u_P / a_P;
        double Ma_N = u_N / a_N;
        
        //计算总焓H
        double H_P = 0.5 * pow(u_P,2) + (gamma/(gamma-1) ) * (p_P/rho_P);
        double H_N = 0.5 * pow(u_N,2) + (gamma/(gamma-1) ) * (p_N/rho_N);

        //计算相关系数：
        double gamma_1 = (gamma-1)/2;
        double gamma_2 = pow(gamma,2) - 1;

        double bP = (1.0/4.0) * rho_P * a_P;
        double bN =  -(1.0/4.0) * rho_N * a_N;
        double dP = pow(1+Ma_P, 2);
        double dN = pow(1-Ma_N, 2);
        double uP = 2.0 * a_P / gamma;
        double uN = 2.0 * a_N / gamma;
        double eP = 2.0 * pow (a_P, 2) /gamma_2;
        double eN = 2.0 * pow (a_N, 2) /gamma_2;


        //计算分裂之后的正负通量
        double rho_FP = bP * dP * 1.0;
        double rho_FN = bN * dN * 1.0;
        double rhou_FP = bP * dP * uP * (gamma_1*Ma_P + 1);
        double rhou_FN = bN * dN * uN * (gamma_1*Ma_N - 1);
        double rhoe_FP = bP * dP * eP * pow(gamma_1*Ma_P + 1, 2);
        double rhoe_FN = bN * dN * eN * pow(gamma_1*Ma_N - 1, 2);

        z[0][j] = rho_FP + rho_FN; 
        z[1][j] = rhou_FP + rhou_FN; 
        z[2][j] = rhoe_FP + rhoe_FN; 
    }   
}


static inline void LiouSteffen_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_P = x[0][j];
        double rho_N = y[0][j];
        double u_P = x[1][j];
        double u_N = y[1][j];
        double p_P = x[2][j];
        double p_N = y[2][j];
        //计算声速
        double a_P = sqrt(gamma * p_P / rho_P);
        double a_N = sqrt(gamma * p_N / rho_N);
        //计算当地Ma
        double Ma_P = u_P / a_P;
        double Ma_N = u_N / a_N;
        //计算总焓H
        double H_P = 0.5 * pow(u_P,2) + (gamma/(gamma-1) ) * (p_P/rho_P);
        double H_N = 0.5 * pow(u_N,2) + (gamma/(gamma-1) ) * (p_N/rho_N);

        //计算Roe平均
        double Ma_PP = (1.0/4.0) * pow(Ma_P+1, 2);
        double Fp_P = 0.5 * p_P * (1+Ma_P);
        if(fabs(Ma_P) > 1.0){
            Ma_PP = 0.5 * (Ma_P + fabs(Ma_P));
            Fp_P = 0.5 * p_P * (Ma_P + fabs(Ma_P)) /Ma_P;
        }

        double Ma_NN = -(1.0/4.0) * pow(Ma_N-1, 2);
        double Fp_N = 0.5 * p_N * (1-Ma_N);
        if(fabs(Ma_N) > 1.0){
            Ma_NN = 0.5 * (Ma_N - fabs(Ma_N));
            Fp_N = 0.5 * p_N * (Ma_N - fabs(Ma_N)) /Ma_N;
        }


        //计算左右守恒变量和通量
        double rho_FL_C = rho_P * a_P;
        double rho_FR_C = rho_N * a_N;
        double rhou_FL_C = rho_P * a_P * u_P;
        double rhou_FR_C = rho_N * a_N * u_N;
        double rhoe_FL_C = rho_P * H_P * a_P;
        double rhoe_FR_C = rho_N * H_N * a_N;

    
        //计算分裂之后的正负通量
        double rho_FP = Ma_PP * rho_FL_C;
        double rho_FN = Ma_NN * rho_FR_C;
        double rhou_FP = Ma_PP * rhou_FL_C + Fp_P;
        double rhou_FN = Ma_NN * rhou_FR_C + Fp_N;
        double rhoe_FP = Ma_PP * rhoe_FL_C;
        double rhoe_FN = Ma_NN * rhoe_FR_C;

        z[0][j] = rho_FP + rho_FN; 
        z[1][j] = rhou_FP + rhou_FN; 
        z[2][j] = rhoe_FP + rhoe_FN; 
    }   
}



static inline void XJF_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_P = x[0][j];
        double rho_N = y[0][j];
        double u_P = x[1][j];
        double u_N = y[1][j];
        double p_P = x[2][j];
        double p_N = y[2][j];
        //计算声速
        double a_P = sqrt(gamma * p_P / rho_P);
        double a_N = sqrt(gamma * p_N / rho_N);
        //计算总焓H
        double H_P = 0.5 * pow(u_P,2) + (gamma/(gamma-1)) * (p_P/rho_P);
        double H_N = 0.5 * pow(u_N,2) + (gamma/(gamma-1)) * (p_N/rho_N);

        //计算Roe平均
        double rhobar = sqrt(rho_P * rho_N);
        double ubar = (sqrt(rho_P) * u_P + sqrt(rho_N) * u_N)/(sqrt(rho_P)+sqrt(rho_N));
        double Hbar = (sqrt(rho_P) * H_P + sqrt(rho_N) * H_N)/(sqrt(rho_P)+sqrt(rho_N));
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));

        //利用Roe平均计算近似波速&&&&同时分裂为正负波速度
        double lambdap1 = u_P - a_P;
        double lambdap2 = u_P;
        double lambdap3 = u_P + a_P;
        double lambdan1 = u_N - a_N;
        double lambdan2 = u_N;
        double lambdan3 = u_N + a_N;
        double lambda1_p = 0.5*(lambdap1  + fabs(lambdap1));
        double lambda2_p = 0.5*(lambdap2  + fabs(lambdap2));
        double lambda3_p = 0.5*(lambdap3  + fabs(lambdap3));
        double lambda1_n = 0.5*(lambdan1  - fabs(lambdan1));
        double lambda2_n = 0.5*(lambdan2  - fabs(lambdan2));
        double lambda3_n = 0.5*(lambdan3  - fabs(lambdan3));
        /*double lambda1 = ubar - cbar;
        double lambda2 = ubar;
        double lambda3 = ubar + cbar;
        double lambda1_p = 0.5*(lambda1  + fabs(lambda1));
        double lambda2_p = 0.5*(lambda2  + fabs(lambda2));
        double lambda3_p = 0.5*(lambda3  + fabs(lambda3));
        double lambda1_n = 0.5*(lambda1  - fabs(lambda1));
        double lambda2_n = 0.5*(lambda2  - fabs(lambda2));
        double lambda3_n = 0.5*(lambda3  - fabs(lambda3));*/

        //计算中间变量
        double F_rhou_PP =  u_P - a_P;
        double F_rhou_PN =  u_P + a_P;
        double F_rhou_NP =  u_N - a_N;
        double F_rhou_NN =  u_N + a_N;
        double F_rhoe_PP =  H_P - u_P *a_P;
        double F_rhoe_PN =  H_P + u_P *a_P;
        double F_rhoe_NP =  H_N - u_N *a_N;
        double F_rhoe_NN =  H_N + u_N *a_N;

        //计算分裂之后的正负通量
        double rho_FP = (0.5*rho_P/gamma) * (lambda1_p + 2.0*(gamma-1)*lambda2_p + lambda3_p);
        double rho_FN = (0.5*rho_N/gamma) * (lambda1_n + 2.0*(gamma-1)*lambda2_n + lambda3_n);
        double rhou_FP = (0.5*rho_P/gamma) * (F_rhou_PP*lambda1_p + 2.0*(gamma-1)*ubar*lambda2_p + F_rhou_PN* lambda3_p);
        double rhou_FN = (0.5*rho_N/gamma) * (F_rhou_NP*lambda1_n + 2.0*(gamma-1)*ubar*lambda2_n + F_rhou_NN* lambda3_n);
        double rhoe_FP = (0.5*rho_P/gamma) * (F_rhoe_PP*lambda1_p + (gamma-1)*pow(u_P,2)*lambda2_p + F_rhoe_PN * lambda3_p);
        double rhoe_FN = (0.5*rho_N/gamma) * (F_rhoe_NP*lambda1_n + (gamma-1)*pow(u_N,2)*lambda2_n + F_rhoe_NN * lambda3_n);


        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        double a_L = sqrt(gamma * p_L / rho_L);
        double a_R = sqrt(gamma * p_R / rho_R);

        //计算需要使用的参数
        //计算声速

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);


        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        //double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
       //double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        //double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double Splus = max_of_two(fabs(u_L)  + a_L, fabs(u_R) + a_R);

        //确定HLL数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        if (sleft >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && sright >0){
            rho_F = rho_FP + rho_FN;
            rhou_F = rhou_FP + rhou_FN;
            rhoe_F = rhoe_FP + rhoe_FN;
        }
        else if (sright <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}


static inline void ER_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        double rho_F = 0,u_F=0,p_F=0;
        double rho_l = x[0][j];
        double rho_r = y[0][j];
        double u_l = x[1][j];
        double u_r = y[1][j];
        double p_l = x[2][j];
        double p_r = y[2][j];
        double c_l = sqrt(gamma*p_l/rho_l);
        double c_r = sqrt(gamma*p_r/rho_r);
        double p_star = ExactRieamnna_pstar(rho_l, rho_r, u_l, u_r, p_l, p_r, gamma, gamma ,1e-10, 1e6);
		double u_star = ExactRieamnna_ustar(p_star, rho_l, rho_r, u_l,u_r, p_l,p_r, gamma,gamma);
        
         
        //分别考虑左激波，左稀疏波的情况
        if (p_star >=  p_l){
            double part1 = (gamma + 1)* p_star + (gamma - 1) * p_l ;
            double part2 = (gamma - 1) * p_star + (gamma + 1) * p_l;
            double rho_star_L = rho_l * part1 / part2;
            double A_L = 2.0 / ((gamma + 1) * rho_l);
            //计算激波的速度
            double B_L = (gamma - 1) * p_l  / (gamma + 1);
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
            double rho_star_L = rho_l * pow(part1/part2, 1/gamma);
            //计算左稀疏波波头和波尾的速度
            double aL_star = c_l * pow((p_star / p_l), ((gamma - 1) / (2 *gamma)));
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
                rho_F = rho_l * pow(2/(gamma + 1) + (gamma - 1) * (u_l - si) / ((gamma + 1) * c_l), (2 / (gamma - 1))); 
                u_F = 2 * (c_l + (gamma - 1) * u_l/2 + si) / (gamma + 1);
                p_F = p_l * pow(2 / (gamma + 1) + (gamma - 1) * (u_l - si) / ((gamma + 1) * c_l),  2 * gamma / (gamma - 1));
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
            double part1 =  (gamma + 1) * p_star + (gamma - 1) * p_r;
            double part2 =  (gamma - 1) * p_star + (gamma + 1) * p_r;
            double rho_star_r = rho_r * part1 / part2;
            double A_r = 2 / ((gamma + 1) * rho_r);
            //计算激波的速度
            double B_r = (gamma - 1) * p_r / (gamma + 1);
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
            double rho_star_r = rho_r * pow(part1 / part2, 1/gamma);
            //计算左稀疏波波头和波尾的速度
            double aR_star = c_r * pow((p_star / p_r), ((gamma - 1) / (2 *gamma)));
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
                rho_F = rho_r * pow(2/(gamma + 1) - (gamma - 1) * (u_r - si) / ((gamma + 1) * c_r), (2 / (gamma - 1))); 
                u_F = 2 * (-c_r + (gamma - 1) * u_r / 2 + si) / (gamma + 1);
                p_F = p_r * pow(2 / (gamma + 1) - (gamma - 1) * (u_r - si) / ((gamma + 1) * c_r),  2 * gamma / (gamma - 1));
            }
            else if (si > S_HR){
                rho_F = rho_r;
                u_F = u_r;
                p_F = p_r;
            }
        }


        z[0][j] = rho_F * u_F; 
        z[1][j] = rho_F * u_F *u_F + p_F; 
        z[2][j] = ((0.5*rho_F * u_F *u_F + p_F/(gamma-1)) + p_F) * u_F; 

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
static inline void Rusanov_Flux_HeatConduction(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {


        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数
        //计算声速
        double a_L = sqrt(gamma*p_L/rho_L);
        double a_R = sqrt(gamma*p_R/rho_R);

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

            //计算温度变量
        double rhobar = sqrt(rho_L*rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;

        //计算一个最大波速
        double Splus = max_of_two(fabs(u_L)  + a_L, fabs(u_R) + a_R);
    
        //Rusanov数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        rho_F = 0.5*(rho_FL + rho_FR) - 0.5*Splus*(rho_R - rho_L);
        rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*Splus*(rhou_R - rhou_L);
//能量方程增加热传导机制
        rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*Splus*(rhoe_R - rhoe_L) - Splus * rhobar*(T_R - T_L);
        
        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}


static inline void Lax_Flux_HeatConduction(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    double S_plus = 0;

    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        //计算声速
        double a_L = sqrt(gamma*p_L/rho_L);
        double a_R = sqrt(gamma*p_R/rho_R);
        //计算全局最大波速
        S_plus = max_of_three(fabs(u_L)+a_L, fabs(u_R)+a_R,S_plus);
    
    }   
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
     

        //计算需要使用的参数
        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);
        //计算温度变量
        double rhobar = sqrt(rho_L*rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;

        //Rusanov数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        rho_F = 0.5*(rho_FL + rho_FR) - 0.5*S_plus*(rho_R - rho_L);
        rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*S_plus*(rhou_R - rhou_L);
        //能量方程增加热传导机制
        rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*S_plus*(rhoe_R - rhoe_L) - S_plus * rhobar*(T_R - T_L);
        
        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLL_Flux_HeatConduction(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma ,double u_point,double dx) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        double a_L = sqrt(gamma * p_L / rho_L);
        double a_R = sqrt(gamma * p_R / rho_R);

        //计算需要使用的参数
        //计算声速

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算温度变量
        double rhobar = sqrt(rho_L*rho_R);
//        double rhobar = 0.5*(rho_L+rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double Splus = max_of_two(fabs(u_L)  + a_L, fabs(u_R) + a_R);

        //确定HLL数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        if (sleft > 0){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft <= 0 && sright >=0){
            rho_F = (sright*rho_FL - sleft*rho_FR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
            rhou_F = (sright*rhou_FL - sleft*rhou_FR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
            rhoe_F = (sright*rhoe_FL - sleft*rhoe_FR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
        }
        else if (sright < 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

         //针对新的相对波速添加相对的热扩散即可
        double S_L = min_of_two(u_L - a_L,u_R - a_R);
        double S_R = max_of_two(u_L + a_L,u_R + a_R);
        double s_L_P = sleft - u_point;
        double s_R_P = sright - u_point;
        double S_LW_L = u_L - a_L;
        double S_LW_R = u_R - a_R;
        double S_RW_L = u_L + a_L;
        double S_RW_R = u_R + a_R;

        if (s_L_P * s_R_P <= 0.1){
            //rho_F = rho_F + sleft * sright /(sright -sleft) * (rho_R - rho_L);
            rhoe_F = rhoe_F -(gamma/(gamma-1)) * fabs(Splus) * (T_R - T_L) * rhobar;
        }
            
        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLLC_Flux_HeatConduction(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma, double u_point) {
    
    for (int j = 1; j <= cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        double a_L = sqrt(gamma * p_L / rho_L);
        double a_R = sqrt(gamma * p_R / rho_R);

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

         //计算温度变量
        double rhobar = sqrt(rho_L*rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
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
            rhoe_F = rhoe_FL + sleft * (rho_star_L * fe_star_L -  rhoe_L);
        }
        else if(s_star < 0 && sright > 0){
            rho_F = rho_FR + sright * (rho_star_R - rho_R);
            rhou_F = rhou_FR + sright *(rho_star_R * s_star  - rhou_R);
            rhoe_F = rhoe_FR + sright * (rho_star_R * fe_star_R -  rhoe_R);
        }
        else if (sright  <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }


        double p_test = p_L * (1.0 + rho_L * u_L * (a_R+u_L)/p_L + gamma * pow(ma,3));
        double p_diff  = p_star-p_test;

        //针对新的相对波速添加相对的热扩散即可
        if (s_L_P * s_R_P < 0)
            rhoe_F = rhoe_F - ((gamma-1)/(gamma-1)) * rhobar *  (T_R - T_L ) * (fabs(s_L_P * s_R_P/(s_R_P-s_L_P)));
            //rho_F = 0.5 * (rho_FL + rho_FR) - 0.5 * Splus * (rho_R - rho_L);

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 

    }   
}



static inline void Roe_Flux_HeatConduction(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma,double u_point) {
    double epsilon = 1e-6;
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        double a_L = sqrt(gamma * p_L / rho_L);
        double a_R = sqrt(gamma * p_R / rho_R);

        //计算需要使用的参数

        //计算热力学变量
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);
        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;

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
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));

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
            rhoe_F = rhoe_F - (1.0/(gamma-1)) * rhobar * (T_R - T_L) * (fabs(s_L_P * s_R_P/(s_R_P-s_L_P)) + fabs(u_point));

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}


static inline void ER_Flux_PlusHeat(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        double rho_F = 0,u_F=0,p_F=0;
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];
        double c_L = sqrt(gamma*p_L/rho_L);
        double c_R = sqrt(gamma*p_R/rho_R);

        double T_L = p_L / rho_L;
        double T_R = p_R / rho_R;
        double Splus = max_of_two(fabs(u_L)  + c_L, fabs(u_R) + c_R);
        double rhoabr = sqrt(rho_L * rho_R);

        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;

        double p_star = ExactRieamnna_pstar(rho_L, rho_R, u_L, u_R, p_L, p_R, gamma, gamma ,1e-4, 1e3);
		double u_star = ExactRieamnna_ustar(p_star, rho_L, rho_R, u_L,u_R, p_L,p_R, gamma,gamma);
        
         
        //分别考虑左激波，左稀疏波的情况
        if (p_star >=  p_L){
            double part1 = (gamma + 1)* p_star + (gamma - 1) * p_L ;
            double part2 = (gamma - 1) * p_star + (gamma + 1) * p_L;
            double rho_star_L = rho_L * part1 / part2;
            double A_L = 2.0 / ((gamma + 1) * rho_L);
            //计算激波的速度
            double B_L = (gamma - 1) * p_L  / (gamma + 1);
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
            double rho_star_L = rho_L * pow(part1/part2, 1/gamma);
            //计算左稀疏波波头和波尾的速度
            double aL_star = c_L * pow((p_star / p_L), ((gamma - 1) / (2 *gamma)));
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
                rho_F = rho_L * pow(2/(gamma + 1) + (gamma - 1) * (u_L - si) / ((gamma + 1) * c_L), (2 / (gamma - 1))); 
                u_F = 2 * (c_L + (gamma - 1) * u_L/2 + si) / (gamma + 1);
                p_F = p_L * pow(2 / (gamma + 1) + (gamma - 1) * (u_L - si) / ((gamma + 1) * c_L),  2 * gamma / (gamma - 1));
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
            double part1 =  (gamma + 1) * p_star + (gamma - 1) * p_R;
            double part2 =  (gamma - 1) * p_star + (gamma + 1) * p_R;
            double rho_star_R = rho_R * part1 / part2;
            double A_R = 2 / ((gamma + 1) * rho_R);
            //计算激波的速度
            double B_R = (gamma - 1) * p_R / (gamma + 1);
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
            double rho_star_R = rho_R * pow(part1 / part2, 1/gamma);
            //计算左稀疏波波头和波尾的速度
            double aR_star = c_R * pow((p_star / p_R), ((gamma - 1) / (2 *gamma)));
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
                rho_F = rho_R * pow(2/(gamma + 1) - (gamma - 1) * (u_R - si) / ((gamma + 1) * c_R), (2 / (gamma - 1))); 
                u_F = 2 * (-c_R + (gamma - 1) * u_R / 2 + si) / (gamma + 1);
                p_F = p_R * pow(2 / (gamma + 1) - (gamma - 1) * (u_R - si) / ((gamma + 1) * c_R),  2 * gamma / (gamma - 1));
            }
            else if (si > S_HR){
                rho_F = rho_R;
                u_F = u_R;
                p_F = p_R;
            }
        }

        z[0][j] = rho_F * u_F; 
        z[1][j] = rho_F * u_F *u_F + p_F ; 
        z[2][j] = ((0.5*rho_F * u_F *u_F + p_F/(gamma-1)) + p_F) * u_F -  Splus * (rhoe_R - rhoe_L) * rhoabr; 
    }   
}


static inline void ER_Flux_Heat(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {

    const int h = 40;
    double Cell_Face_U[rows][h],Cell_Face_FU[rows][h];
    double CFU_L[rows][h],CFU_R[rows][h];
    double CFp_L[rows][h],CFp_R[rows][h];
    double dxi,dtau;
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //使用数值方法计算界面数值通量：
        for (int k = 0; k < h; k++)
        {
            if (k < h/2){
                Cell_Face_U[0][k] = rho_L;
                Cell_Face_U[1][k] = rho_L * u_L;
                Cell_Face_U[2][k] = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma - 1);
            }
            else{
                Cell_Face_U[0][k] = rho_R;
                Cell_Face_U[1][k] = rho_R * u_R;
                Cell_Face_U[2][k] = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma - 1);
            }
            
        }
        dxi = 1.0;
        dtau = 0.3;
        for (int Ite = 0; Ite < h; Ite++)
        {
            //计算数值通量：
            for (int m = 0; m < rows; m++){
                for (int n = 0; n < h; n++){
                    CFU_L[m][n] = Cell_Face_U[m][n];
                    CFU_R[m][n] = Cell_Face_U[m][n+1];
                }
            }

            Con_to_Pri_1D(3,h,CFp_L,CFU_L,gamma);
            Con_to_Pri_1D(3,h,CFp_R,CFU_R,gamma);
            HLL_Flux(3, h, CFp_L,CFp_R,Cell_Face_FU,gamma);

            for (int m = 0; m < rows; m++){
                for (int n = 1; n < h-1; n++){
                    Cell_Face_U[m][n] = Cell_Face_U[m][n] + dtau / dxi * (Cell_Face_FU[m][n] - Cell_Face_FU[m][n-1]);
                }
            }

            //执行边界条件：
            for (int m = 0; m < rows; m++){
                Cell_Face_U[m][0] = Cell_Face_U[m][1];
                Cell_Face_U[m][h-1] = Cell_Face_U[m][h-2];
            }
        }

        //提取数值通量
        double rho_cf =  Cell_Face_U[0][h/2-1] ;  
        double rhou_cf = Cell_Face_U[1][h/2-1] ;  
        double rhoe_cf = Cell_Face_U[2][h/2-1] ; 
        
        double u_cf = rhou_cf/rho_cf;
        double p_cf = (gamma-1) *(rhoe_cf - 0.5 * rhou_cf * u_cf); 
        double rho_F = rhou_cf;
        double rhou_F = rho_cf * pow(u_cf,2) + p_cf;
        double rhoe_F = (0.5 * rho_cf * pow(u_cf,2) + p_cf * gamma /(gamma-1))*u_cf;
        
        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 

    }   
}





/*                               ************************************                               */
/*                               ************************************                               */
/*                              Riemann Solver：Exact and Approximate                               */
/*                               ************************************                               */
/*                               ************************************                               */

/*                                      ******************                                          */
/*                                   在相对参考系下的黎曼解法器                                        */
/*                                      ******************                                          */
//Flux计算方法
static inline void HLL_Flux_XRela(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma, double u_r) {
    
    int h;
    for (int j = 1; j < cols-1; j++) {


        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j] - u_r;
        double u_R = y[1][j] - u_r;
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数
        //计算声速
        double a_L = sqrt(gamma * p_L/rho_L);
        double a_R = sqrt(gamma * p_R/rho_R);

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double S_LL = u_L - a_L;
        double S_LR = u_R - a_R;
        double S_RL = u_L + a_L;
        double S_RR = u_R + a_R;


        //变换到原本参考系
        double rho_FL = rho_L * (u_L + u_r);
        double rho_FR = rho_R * (u_R + u_r);
        double rhou_L = rho_L * (u_L + u_r);
        double rhou_R = rho_R * (u_R + u_r);
        double rhou_FL = rho_L * pow(u_L + u_r, 2) + p_L;
        double rhou_FR = rho_R * pow(u_R + u_r, 2) + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L+u_r, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R+u_r, 2) + p_R/(gamma-1);
        double rhoe_FL = (0.5 * rho_L * pow(u_L+u_r, 2) + p_L/(gamma-1) + p_L) * (u_L+u_r);
        double rhoe_FR = (0.5 * rho_R * pow(u_R+u_r, 2) + p_R/(gamma-1) + p_R) * (u_R+u_r);

        //计算一个最大波速
        double Splus = max_of_two(fabs(u_L) +fabs(u_r) + a_L, fabs(u_R) + fabs(u_r) + a_R);
        //计算最大流体速度
        double u_plus =  max_of_two(fabs(u_L) +fabs(u_r), fabs(u_R) + fabs(u_r));
 

        //确定HLL数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        if ( min_of_two(S_LL,S_LR) >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if (max_of_two(S_RL,S_RR) <= 0)
        {
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }
        
        else if(S_LL * S_LR < 0  && S_RL*S_RR > 0){
            if (S_LL < S_LR)
            {
                rho_F = rho_FL;
                rhou_F = rhou_FL;
                rhoe_F = rhoe_FL;
            }
            else {
                rho_F = 0.5*(rho_FL + rho_FR) - 0.5*Splus*(rho_R - rho_L);
                rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*Splus*(rhou_R - rhou_L);
                rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*Splus*(rhoe_R - rhoe_L);
            }

        }

        else if(S_LL * S_LR > 0  && S_RL*S_RR < 0){
            if (S_RL < S_RR)
            {
                rho_F = rho_FR;
                rhou_F = rhou_FR;
                rhoe_F = rhoe_FR;
            }
            else {
                rho_F = 0.5*(rho_FL + rho_FR) - 0.5*Splus*(rho_R - rho_L);
                rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*Splus*(rhou_R - rhou_L);
                rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*Splus*(rhoe_R - rhoe_L);
            }
        }
        else{
            rho_F = 0.5*(rho_FL + rho_FR) - 0.5*Splus*(rho_R - rho_L);
            rhou_F = 0.5*(rhou_FL + rhou_FR) - 0.5*Splus*(rhou_R - rhou_L);
            rhoe_F = 0.5*(rhoe_FL + rhoe_FR) - 0.5*Splus*(rhoe_R - rhoe_L);
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLLC_Flux_Rela(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double s_star = (p_R - p_L + rho_L*u_L*(sleft - u_L) - rho_R*u_R*(sright - u_R))\
                                 /(rho_L*(sleft - u_L) - rho_R*(sright - u_R));

        double u_stat_L = rho_L * (sleft-u_L)/(sleft-s_star);
        double u_stat_R = rho_R * (sright-u_R)/(sright-s_star);
        double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L)));
        double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R)));


        //HLLC数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        if (sleft >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && s_star >=0 ){
            rho_F = rho_FL + sleft * (u_stat_L - rho_L);
            rhou_F = rhou_FL + sleft *(u_stat_L * s_star  - rhou_L);
            rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L -  rhoe_L);

        }
        else if(s_star < 0 && sright > 0){
            rho_F = rho_FR + sright * (u_stat_R - rho_R);
            rhou_F = rhou_FR + sright *(u_stat_R * s_star  - rhou_R);
            rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R -  rhoe_R);
        }
        else if (sright <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



#endif  