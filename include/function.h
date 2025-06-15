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


static inline double Get_lamdat(int cols,double (*vec)[cols], double dx,double cfl,double gamma) {

    double dt;
    double lammax=0;
    double rho[cols], u[cols], p[cols], c[cols];

    for (int i = 0; i < cols; i++)
    {
        rho[i] = vec[0][i];
        u[i] = vec[1][i]/vec[0][i];
        p[i] = (vec[2][i]-0.5*rho[i]*pow(u[i],2))*(gamma-1);
        c[i] = sqrt (gamma*(p[i]/rho[i]));
    }

    for (int i = 0; i < cols; i++)
    {
        lammax = fmax(lammax,fabs(u[i])+fabs(c[i]));
    }

    dt = cfl*dx/lammax;
    
    return dt;
}

static inline double max_of_four(double a, double b, double c, double d) {
    double max_value = a;
    
    if (b > max_value) {
        max_value = b;
    }
    
    if (c > max_value) {
        max_value = c;
    }
    
    if (d > max_value) {
        max_value = d;
    }
    
    return max_value;
}

static inline double max_of_two(double a, double b) {
    if (a > b)
        return a;
    else 
        return b;
}

static inline double min_of_two(double a, double b) {
    if (a > b)
        return b;
    else 
        return a;
}


#endif