#ifndef INITIALIZATION_H
#define INITIALIZATION_H


#include <stdio.h>
#include "Golbal.h"
#include "math.h"
#include "function.h"


// 定义圆周率PI常量
#ifndef PI
    #define PI 3.14159265358979323846
#endif


static inline void initialize(double *x,int n){
    int i;
    for(i = 0;i < n;i++){
        x[i]=0.0;
    }
}


static inline void initEuler1D(int rows, int cols,double (*x)[cols]) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            x[i][j] = 0;
        }
    }
}

static inline void initEuler2D(int rows, int cols, int depth, double (*x)[cols][depth]) {
    int i, j, k;
    for (i = 0; i < rows; i++) 
        for (j = 0; j < cols; j++) 
           for (k = 0; k < depth; k++) 
                x[i][j][k] = 0.0;

}

static inline void initEulerpri1D_Sod(int rows, int cols,double (*x)[cols]) {
    int i, j;
    for (i = 0; i < rows; i++) {
        if(i == 0){
            for (j = 0; j < cols; j++) {
                if( j <cols/2){
                    x[i][j] = 1;
                }
                else{
                    x[i][j] = 0.125;
                }
            }
        }
        if(i == 1){
            for (j = 0; j < cols; j++) {
                if( j <cols/2){
                    x[i][j] = 0;
                }
                else{
                    x[i][j] = 0;
                }
            }
        }
        if(i == 2){
            for (j = 0; j < cols; j++) {
                if( j <cols/2){
                    x[i][j] = 1.0;
                }
                else{
                    x[i][j] = 0.1;
                }
            }
        }
    }
}


static inline void initEulerpri1D_Smooth(int rows, int cols, int GC, double (*x)[cols], double Deltax) {

    for (int j = 0; j < cols; j++) {
        x[0][j] = 2.0 - 1.0/( 4.0 * PI )*(cos(4.0 * PI * (j+1-GC) * Deltax)-cos(4.0 * PI * (j-GC) * Deltax))/Deltax;
        x[1][j] = 1.0;
        x[2][j] = 1.0;
    }

}



static inline void initEulerpri1D_Osher(int rows, int cols, int GC, double (*pri)[cols], double Deltax) {
    int i, j;

    int NGC= 2*GC;

    for (i = 0; i < rows; i++) {
        if(i == 0){
            pri[i][0] = 3.857;
            pri[i][1] = 3.857;
            for (j = 0; j < cols; j++) {
                if( j < (cols - NGC)/10){
                    pri[i][j] = 3.857;
                }
                else{
                    pri[i][j] = 1. + 0.2 * sin(5*((j-GC)*Deltax-5));
                }
            }

        }
        if(i == 1){
            for (j = 0; j < cols; j++) {
                if( j < (cols-NGC)/10){
                    pri[i][j] = 2.629;
                }
                else{
                    pri[i][j] = 0;
                }
            }
        }
        if(i == 2){
            for (j = 0; j < cols; j++) {
                if( j < (cols - NGC)/10){
                    pri[i][j] = 10.333;
                }
                else{
                    pri[i][j] = 1.0;
                }
            }
        }
    }
}

static inline void initEulerpri1D_Shocktube(int rows, int cols,double (*x)[cols], double pri1[3],double pri2[3]) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if( j <cols/2)
                x[i][j] = pri1[i];
            else
                x[i][j] = pri2[i];
        }
    }
}

static inline void initEulerpri2D_Shocktube(int rows, int cols, int depth, double (*x)[cols][depth], double pri1[rows],double pri2[rows]) {
    int i, j, k;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            for ( k = 0; k < depth; k++){
                if( j < cols/2)
                    x[i][j][k] = pri1[i];
                else
                    x[i][j][k] = pri2[i];
            }
        }
    }
}


static inline void initEulerConser_Smooth(int rows, int cols, int GC, double (*y)[cols], double Deltax) {

    double u[10];
    for (int h = 0; h < 5; h++)
        u[h] = 0.0;

    for (int j = 0; j < cols; j++) {
        y[0][j] = 2.0 - 1.0/( 4.0 * PI )*(cos(4.0 * PI * (j+1-GC) * Deltax)-cos(4.0 * PI * (j-GC) * Deltax))/Deltax;
        y[1][j] = y[0][j] * 1.0;
        y[2][j] = 0.5*y[0][j] * 1.0 + 1.0/(M_gamma-1);
    }

}

static inline void initEulerconser1D(int rows, int cols,double (*x)[cols],double (*y)[cols]) {
    int j;
    for (j = 0; j < cols; j++) {
        y[0][j] = x[0][j];
        y[1][j] = x[1][j] * x[0][j];
        y[2][j] = x[2][j]/(M_gamma-1) + 0.5 * x[0][j] * x[1][j] * x[1][j];
    }
}

static inline void initEulerconser2D(int rows, int cols, int depth, double (*x)[cols][depth], double (*y)[cols][depth]) {
    int j,k;
    for (j = 0; j < cols; j++) {
        for ( k = 0; k < depth; k++){
            double rho, u, v, p;
            rho = x[0][j][k];
            u = x[1][j][k];
            v = x[2][j][k];
            p = x[3][j][k];
            y[0][j][k] = rho;
            y[1][j][k] = u * rho;
            y[2][j][k] = v * rho;
            y[3][j][k] = p /(M_gamma-1) + 0.5 * rho * (pow(u,2) + pow(v,2));
        }
    }
}



static inline void initEulerflux1D(int rows, int cols,double (*y)[cols],double (*z)[cols]) {
    int j;
    for (j = 0; j < cols; j++) {
        z[0][j] = y[1][j];
        z[1][j] = y[1][j]*y[1][j]/y[0][j] + (M_gamma - 1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j]);
        z[2][j] = (y[2][j] + (M_gamma - 1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j])) * y[1][j]/y[0][j];
    }
}

static inline void initEulerflux2D(int rows, int cols, int depth, double (*y)[cols][depth], double (*f)[cols][depth],double (*g)[cols][depth]) {
    int j,k;
    for (j = 0; j < cols; j++) {
        for ( k = 0; k < depth; k++){
            double rho, u, v, p;
            rho = y[0][j][k];
            u = y[1][j][k]/y[0][j][k];
            v = y[2][j][k]/y[0][j][k];
            p = (y[3][j][k] - 0.5 * rho * (pow(u,2) + pow(v,2)))*(M_gamma-1);
            //x dir initilize
            f[0][j][k] = rho * u;
            f[1][j][k] = rho * u * u + p;
            f[2][j][k] = rho * u * v;
            f[3][j][k] = (y[3][j][k] + p) * u;
            //y dir initilize
            g[0][j][k] = rho * v;
            g[1][j][k] = rho * u * v;
            g[2][j][k] = rho * v * v + p;
            g[3][j][k] = (y[3][j][k] + p) * v;

        }
    }
}



static inline void Con_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j <= cols-1; j++) {
        x[0][j] = y[0][j];
        x[1][j] = y[1][j]/y[0][j];
        x[2][j] = (M_gamma-1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j]);
    }
}

static inline void Con_to_Pri_2D(int rows, int cols, int depth, double (*x)[cols][depth], double (*y)[cols][depth]){
   int j,k;
    for (j = 0; j < cols; j++) {
        for ( k = 0; k < depth; k++){
            double rho, u, v, p;
            rho = y[0][j][k];
            u = y[1][j][k]/y[0][j][k];
            v = y[2][j][k]/y[0][j][k];
            p = (y[3][j][k] - 0.5 * rho * (pow(u,2) + pow(v,2)))*(M_gamma-1);
            //x dir initilize
            x[0][j][k] = rho;
            x[1][j][k] = u;
            x[2][j][k] = v;
            x[3][j][k] = p;
        }
    }
}


static inline void Pri_to_Con_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j <= cols-1; j++) {
        y[0][j] = x[0][j];
        y[1][j] = x[0][j] * x[1][j];
        y[2][j] = 0.5 * x[0][j] * pow(x[1][j],2) + x[2][j]/(M_gamma-1);
    }
}





//初始条件
static inline void Init_Sod(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 1.0;
    double rho2 = 0.125;
    double u1 = 0.0;
    double u2 = 0.0;
    double p1 = 1.0;
    double p2 = 0.1;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


static inline void Init_Sod_2D(double pri_Ver1[4], double pri_Ver2[4]) {
    double rho1 = 1.0;
    double rho2 = 0.125;
    double u1 = 0.0;
    double u2 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double p1 = 1.0;
    double p2 = 0.1;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = v1;
    pri_Ver1[3] = p1;


    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = v2;
    pri_Ver2[3] = p2;
}

static inline void Init_Sod_Rare(double pri_Ver1[3], double pri_Ver2[3]) {


    double rho1 = 1.0;
    double rho2 =  0.43757818061324344;
    //double rho2 =  0.5;
    double u1 = 0.0;
    double u2 =  0.9013775087441291 ;
    double p1 = 1.0;
    double p2 =  0.31439665844271514;
    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

static inline void Init_Sod_Shock(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 0.2655737117053071;
    double rho2 =  0.125;
    double u1 = 0.9274526200489499;
    double u2 =  0.0;
    double p1 =  0.30313017805064685;
    double p2 = 0.1;
    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


static inline void Init_DRare(double pri_Ver1[3], double pri_Ver2[3],double u_r) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = -2.0 + u_r;
    double u2 = 2.0 + u_r;
    double p1 = 1.0;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

static inline void Init_Shock_Impact(double pri_Ver1[3], double pri_Ver2[3],double u_r) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = 4.0 + u_r;
    double u2 = -4.0 + u_r;
    double p1 = 1.0;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}



static inline void Init_Shock1(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 24.0/11;
    double rho2 = 1.0;
    double u1 = 13.0/12;
    double u2 = 0;
    double p1 = 19.0/6;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


static inline void Init_Shock2(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 600.0/107.0;
    double rho2 = 1.0;
    double u1 = 493.0/60.0;
    double u2 = 0.0;
    double p1 = 499.0/6.0;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

static inline void Init_Shock3(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = 4.0;
    double u2 = -4.0;
    double p1 = 1.0;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


static inline void Init_Euler(int rows,int cols, int GC, double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double pri_Ver1[3], double pri_Ver2[3], double Deltax){
    //初始化
    initEuler1D(3, cols, x);
    initEuler1D(3, cols, y);
    initEuler1D(3, cols, z);
//    initEulerConser_Smooth(3, cols,GC, y, Deltax);
    //进一步初始化
    initEulerpri1D_Shocktube(3, cols, x,pri_Ver1,pri_Ver2);
//    initEulerpri1D_Osher(3, cols,6,x, 10.0/(cols-6));
    initEulerconser1D(3, cols,x, y);
    initEulerflux1D(3, cols, y,z);
}

static inline void Init_Euler_2D(int rows,int cols, int depth, int GC, double (*x)[cols][depth],double (*y)[cols][depth],\
                                    double (*f)[cols][depth], double (*g)[cols][depth], double Ver_a[rows], double Ver_b[rows], double Deltax){
    //初始化
    initEuler2D(rows, cols, depth,x);
    initEuler2D(rows, cols, depth,y);
    initEuler2D(rows, cols, depth,f);
    initEuler2D(rows, cols, depth,g);
   
//    initEulerConser_Smooth(3, cols,GC, y, Deltax);
    //进一步初始化
    initEulerpri2D_Shocktube(rows, cols,depth, x, Ver_a, Ver_b);
//    initEulerpri1D_Osher(3, cols,6,x, 10.0/(cols-6));
    initEulerconser2D(rows, cols, depth, x, y);
    initEulerflux2D(rows, cols, depth, y,f,g);
}



#endif