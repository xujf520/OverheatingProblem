#ifndef ERROR_H
#define ERROR_H

#include <stdio.h>
#include <math.h>
#include "Reconstruction.h"
#include "function.h"
#include "Golbal.h"



//计算有限体积格式误差L1误差
static inline void Scheme_Error_L1(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    

    for (int i = 0; i < rows; i++)
    {
        Error[i] = 0.0;
    }
    

    for (int i = 0; i < rows; i++){
        double E1 = 0;
        for (int j = GC ; j <= cols-GC-2; j++){
            double u[10];
            for (int h = 0; h < 10; h++)
                u[h] = x[i][j+h-2];

            double Inter = delta_x * numerical_integration(Smooth_function,(j-GC)*delta_x, (j-GC+1.0)*delta_x, 1000, 2, &u[2]);
            E1 = E1 + Inter;
        } 
        Error[i] = E1;
    }

    printf("Mass_Error = %f \n", Error[0]);
}


//计算有限体积格式误差L2误差
static inline void Scheme_Error_L2(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    
    int l =  cols- 2*GC;

    for (int i = 0; i < rows; i++)
    {
        Error[i] = 0.0;
    }
    

//    WENO3_Reconstruction(rows,cols,GC,x,Error_L,Error_R);

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

//计算有限体积格式误差L∞误差
static inline void Scheme_Error_LinFin(int rows, int cols,  int GC, double (*x)[cols], double delta_x) {

    double Error[rows];
    double Error_L[rows][cols],Error_R[rows][cols];
    
    int l =  cols- 2*GC;

    for (int i = 0; i < rows; i++)
        Error[i] = 0.0;


    for (int i = 0; i < rows; i++){
        double E1 = 0;
        double Error_max = 0.0;
        for (int j = GC ; j <= cols-GC-2; j++){
            Error_max = max_of_two(Error_max, fabs(x[i][j] + 1.0/( 4.0 * PI )*(cos(4.0 * PI * (j+1-GC) * delta_x)-cos(4.0 * PI * (j-GC) * delta_x))/delta_x - 2.0));
        } 
        Error[i] = Error_max;
    }

    printf("Mass_Error = %f \n", Error[0]);
}

#endif
