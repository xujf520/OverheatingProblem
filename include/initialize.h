#ifndef INITIALIZATION_H
#define INITIALIZATION_H


#include <stdio.h>
#include "math.h"
#include "function.h"

static inline void initialize(double *x,int n){
    int i;
    for(i = 0;i < n;i++){
        x[i]=0.0;
    }
}


static inline void temp_init(double *temp_init1,double *temp_init2,double *temp_init3,double n){
    int i;
    for(i = 0;i < n;i++){
        temp_init1[i]=0;
    }
    for(i = 0;i < n;i++){
        temp_init2[i]=0;
    }
    for(i = 0;i < n;i++){
        temp_init3[i]=0;
    }

}



static inline void initialize_1(double*u1,double*u2,double*u3,\
                        int n,double ro1,double e1, double ro2,double e2,double u){
    int i;
    for(i = 0;i < n;i++){
        if (i < n/2){
            u1[i] =ro1;
        }else{
            u1[i] =ro2;
        }    
    }
    for(i = 0;i < n;i++){
        if (i < n/2){
            u2[i] =ro1*u;
        }else{
            u2[i] =ro2*u;
        }    
    }
    for(i = 0;i < n;i++){
        if (i < n/2){
            u3[i] =e1;
        }else{
            u3[i] =e2;
        }    
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



static inline void initEulerpri1D_Osher(int rows, int cols,double (*pri)[cols], double (*x)) {
    int i, j;
    for (i = 0; i < rows; i++) {
        if(i == 0){
            pri[i][0] = 3.857;
            pri[i][1] = 3.857;
            for (j = 2; j < cols-2; j++) {
                if( j < (cols-4)/10){
                    pri[i][j] = 3.857;
                }
                else{
                    pri[i][j] = 1. + 0.2 * sin(5*(x[j-2]-5));
                }
            }
            pri[i][cols-2] = pri[i][cols-3];
            pri[i][cols-1] = pri[i][cols-3];
        }
        if(i == 1){
            for (j = 0; j < cols; j++) {
                if( j < (cols-4)/10){
                    pri[i][j] = 2.629;
                }
                else{
                    pri[i][j] = 0;
                }
            }
        }
        if(i == 2){
            for (j = 0; j < cols; j++) {
                if( j < (cols-4)/10){
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
        if(i == 0){
            for (j = 0; j < cols; j++) {
                if( j <cols/2){
                    x[i][j] = pri1[i];
                }
                else{
                    x[i][j] = pri2[i];
                }
            }
        }
        if(i == 1){
            for (j = 0; j < cols; j++) {
                if( j <cols/2){
                    x[i][j] = pri1[i];
                }
                else{
                    x[i][j] = pri2[i];
                }
            }
        }
        if(i == 2){
            for (j = 0; j < cols; j++) {
                if( j <cols/2){
                    x[i][j] = pri1[i];
                }
                else{
                    x[i][j] = pri2[i];
                }
            }
        }
    }
}


static inline void initEulerconser1D(int rows, int cols,double (*x)[cols],double (*y)[cols], double gamma) {
    int j;
    for (j = 0; j < cols; j++) {
        y[0][j] = x[0][j];
        y[1][j] = x[1][j] * x[0][j];
        y[2][j] = x[2][j]/(gamma-1) + 0.5 * x[0][j] * x[1][j] * x[1][j];
    }
}

static inline void init_ConS(int rows, int cols,double (*x)[cols],double (*y)[cols]) {
    int j;
    for (j = 0; j < cols; j++) {
        y[0][j] = x[0][j];
        y[1][j] = x[0][j]*x[1][j];
        y[2][j] = (x[2][j]/pow(x[0][j], 1.4))*x[0][j];
    }
}

static inline void initEulerflux1D(int rows, int cols,double (*y)[cols],double (*z)[cols],double gamma) {
    int j;
    for (j = 0; j < cols; j++) {
        z[0][j] = y[1][j];
        z[1][j] = y[1][j]*y[1][j]/y[0][j] + (gamma - 1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j]);
        z[2][j] = (y[2][j] + (gamma - 1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j])) * y[1][j]/y[0][j];
    }
}


static inline void init_FluxS(int rows, int cols,double (*x)[cols],double (*z)[cols]) {
    int j;
    for (j = 0; j < cols; j++) {
        z[0][j] = x[0][j] * x[1][j];
        z[1][j] = x[0][j] * pow(x[1][j],2) +x[2][j];
        z[2][j] = x[0][j] * x[1][j] * (x[2][j]/pow(x[0][j], 1.4));
    }
}





#endif