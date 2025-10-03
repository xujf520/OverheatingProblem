#ifndef FUNCTION_H
#define FUNCTION_H

#include <math.h>

#ifndef PI
    #define PI 3.14159265358979323846
#endif

static inline int sgn(double num) {
    if (num > 0) {
        return 1;
    } else if (num < 0) {
        return -1;
    } else {
        return 0;
    }
}


static inline double max_of_two(double a, double b) {
    if (a > b)
        return a;
    else 
        return b;
}



static inline double max_of_three(double a, double b, double c) {
    a = max_of_two(a,b);
    return max_of_two(a,c);
}


static inline double max_of_four(double a, double b, double c, double d) {
    a= max_of_three(a,b,c);
    return max_of_two(a,d);
}

static inline double min_of_two(double a, double b) {
    if (a > b)
        return b;
    else 
        return a;
}


static inline void Con_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j <= cols-1; j++) {
        x[0][j] = y[0][j];
        x[1][j] = y[1][j]/y[0][j];
        x[2][j] = (M_gamma-1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j]);
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



static inline double Get_Delta_T(int rows, int cols,double (*x)[cols], double dx, double CFL) {
    double S_plus = 0;
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho = x[0][j];
        double rhou = x[1][j];
        double rhoe = x[2][j];
        double u = rhou / rho;
        double p = (rhoe - 0.5 * rho * pow(u, 2))*(M_gamma-1);
        //计算声速
        double a= sqrt(M_gamma * p / rho);
        //计算全局最大波速
        S_plus = max_of_two (fabs(u) + a,S_plus);
    
    }   

//    return CFL * dx / S_plus;
//    return CFL * dx / S_plus;
    return  pow(dx,5.0/3.0);
}

//计算总守恒量
static inline void Total_Conser(int rows, int cols, double Ghost_Cell , double (*x)[cols], double Conser[rows],double delta_x) {
    
    for (int i = 0; i < rows; i++){
        for (int j = Ghost_Cell ; j < cols - Ghost_Cell; j++){
            Conser[i] += x[i][j] * delta_x;  
        } 
        
    }
}



#endif