#ifndef TIME_ADVANCE_H
#define TIME_ADVANCE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "initialize.h"
#include "scheme.h"

void OutputFluxData_file();

                            /*1阶Euler时间推进方法*/


 //TVD重构配合不同Riemann Solver
static inline void RK1_TVD(int Recon_Accur, int AR_scheme, int rows, int cols , double (*x_d), double (*y)[cols]\
                            ,double (*z)[cols],double dt,double dx,double cfl ,double gamma, int k,double u_r) {
    int i,j;
    //计算时间步长
    //*dt =Get_lamdat(cols,y,dx,cfl,gamma);
    

    //AR_scheme is Approximate Riemann Solver
     //数值通量
    switch (AR_scheme) {
        case 0:
            RS_HLL(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma,u_r);
            break;
        case 1:
            RS_HLLC(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma);
            break;
        case 2:
            RS_Roe(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma);
            break;
        case 3:
            RS_ER(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma);
            break;
        case 4:
            RS_Marquina(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma);
            break;
        default:
            // 这里可以处理 scheme 不是 0、1、2 的情况
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }

    // 第一步计算
/*                                            *********                                          */
/*                                            输出数值通量                                        */
/*                                            *********                                          */
//    OutputFluxData_file(k, cols, dx, x_d, z);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            y[i][j] = y[i][j] - dt*(z[i][j]-z[i][j-1])/dx;
        }
    }

    //边界条件
    for ( i = 0; i < rows; i++)
    {
        if (i == 1)//条件意思是速度反射
        {
            y[i][0] = y[i][1];
            y[i][cols-1] = y[i][cols-2];
        }else{
            y[i][0] = y[i][1];
            y[i][cols-1] = y[i][cols-2];
        }
        
    }
}


static inline void RK1_TVD_FluxRela(int Recon_Accur, int Rela_scheme, int rows, int cols , double (*x_d), double (*y)[cols]\
                            ,double (*z)[cols],double dt,double dx,double cfl ,double gamma, int k, double u_r) {
    int i,j;
    //计算时间步长
    //*dt =Get_lamdat(cols,y,dx,cfl,gamma);
    

    //Rela_scheme ：是否开启相对运动状态的Riemann solver
    //数值通量
    switch (Rela_scheme) {
        case 0:
            RS_HLL_XRela(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma,u_r);
            break;
        case 1:
            printf("代码未完成\n");
            exit(1);
            break;
        default:
            RS_Marquina_XRela(Recon_Accur,3,cols,y,z,dt,dx,cfl,gamma,u_r);
            break;
    }

    // 第一步计算
/*                                            *********                                          */
/*                                            输出数值通量                                        */
/*                                            *********                                          */
//    OutputFluxData_file(k, cols, dx, x_d, z);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            y[i][j] = y[i][j] - dt*(z[i][j]-z[i][j-1])/dx;
        }
    }

    //边界条件
    for ( i = 0; i < rows; i++)
    {
        if (i == 1)//条件意思是速度反射
        {
            y[i][0] = y[i][1];
            y[i][cols-1] = y[i][cols-2];
        }else{
            y[i][0] = y[i][1];
            y[i][cols-1] = y[i][cols-2];
        }
        
    }
}


void OutputFluxData_file(int k, int cols, double dx, double* x,double (*z)[cols]) {
    char filename[100];
    sprintf(filename, "D:/Desktop/Single_Med_data/flux_output_%d.dat", k);
    // 打开文件
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("无法打开文件 %s\n", filename);
        return; // 返回错误代码
    }
    // 写入文件头部信息
    fprintf(file, "variables =  'x' 'flux_rho' 'flux_rhou' 'flux_E' \n");
    // 写入计算结果
    for (int i = 2; i < cols-2 ; i++) {
        fprintf(file,"%.8f\t%.8f\t%.8f\t%.8f\n",x[i-2]+0.5*dx, z[0][i],z[1][i], z[2][i]);
    }
    fclose(file);
    printf("numerical flux ：The %d th calculation ended\n", k);
}

#endif 
