#ifndef FUNCTION_H
#define FUNCTION_H

#include <math.h>

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



static inline double Get_Delta_T(int rows, int cols,double (*x)[cols], double dx, double CFL, double gamma) {
    double S_plus = 0;
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho = x[0][j];
        double rhou = x[1][j];
        double rhoe = x[2][j];
        double u = rhou / rho;
        double p = (rhoe - 0.5 * rho * pow(u, 2))*(gamma-1);
        //计算声速
        double a= sqrt(gamma * p / rho);
        //计算全局最大波速
        S_plus = max_of_two(fabs(u)+a,S_plus);
    
    }   

    return CFL * dx / S_plus;
}


#endif