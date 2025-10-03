#ifndef TIME_ADVANCE_H
#define TIME_ADVANCE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Golbal.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "initialize.h"
#include "scheme.h"

int Time_Step = 0;
//static double mass_error = 0.0;
//static double momentum_error = 0.0;
//static double energy_error = 0.0;

void OutputFluxData_file();
void OutputConservationErrors_file(); 

                            /*1阶Euler时间推进方法*/


//TVD重构配合不同Riemann Solver
static inline void RK1_TimeAd(int AR_scheme, int rows, int cols, int GC, double (*x_d), double (*y)[cols], double dt,double dx, double u_r) {
    int i,j;
    double Conser_1[3]={0.0}, Conser_2[3]={0.0};

    double Test[rows][cols];
    double Flux[rows][cols];
    
    //初始化参数
    for ( i = 0; i < rows; i++){
        for ( j = 0; j < cols; j++){
            Flux[i][j] = y[i][j];
        }
        
    }


    //施加边界条件，
//    BC_OutFlow(rows,cols,y,0,GC);                       //边界条件说明见具体子程序
//    BC_OutFlow(rows,cols,y,1,GC);
    BC_Periodicity(rows,cols,y,0,GC);
    BC_Periodicity(rows,cols,y,1,GC);
    for ( i = 0; i < rows; i++)
    {
        for ( j = 0; j < cols; j++)
        {
            Test[i][j] = y[i][j];
        }
        
    }
    

    Flux_Reconstruction_RP(AR_scheme,rows,cols,GC,y,Flux,dt,dx);

    // 第一步计算
/*                                            *********                                          */
/*                                            输出数值通量                                        */
/*                                            *********                                          */
//    OutputFluxData_file(k, cols, dx, x_d, z);

    Total_Conser(rows,cols,2,y,Conser_1,dx);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC; j++){
            y[i][j] = y[i][j] - dt*(Flux[i][j]-Flux[i][j-1])/dx;
        }
    }
    Total_Conser(rows,cols,2,y,Conser_2,dx);
//    OutputConservationErrors_file(Time_Step,Conser_1,Conser_2);
    Time_Step++;
}


//TVD重构配合不同Riemann Solver
static inline void RK1_TimeAd_RPHeat(int AR_scheme, int rows, int cols, int GC, double (*x_d), double (*y)[cols],double dt,double dx, double u_r) {
    int i,j;

    double Flux[rows][cols];
    
    //初始化参数
    for ( i = 0; i < rows; i++){
        for ( j = 0; j < cols; j++){
            Flux[i][j] = y[i][j];
        }
        
    }


    //施加边界条件，
//    BC_OutFlow(rows,cols,y,0,GC);                    
//    BC_OutFlow(rows,cols,y,1,GC);
    BC_Periodicity(rows,cols,y,0,GC);
    BC_Periodicity(rows,cols,y,1,GC);

    //AR_scheme is Approximate Riemann Solver
    //数值通量
    Flux_Reconstruction_RP_Heat(AR_scheme,rows,cols,GC,y,Flux,dt,dx,u_r);

    // 第一步计算
/*                                            *********                                          */
/*                                            输出数值通量                                        */
/*                                            *********                                          */
//    OutputFluxData_file(k, cols, dx, x_d, z);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            y[i][j] = y[i][j] - dt*(Flux[i][j] - Flux[i][j-1])/dx;
        }
    }
    
}



static inline void RK1_TVD_FluxRela(int Rela_scheme, int rows, int cols , double (*x_d), double (*y)[cols]\
                                            ,double (*z)[cols],double dt,double dx, int k, double u_r) {
    int i,j;
    //计算时间步长
    //*dt =Get_lamdat(cols,y,dx,CFL,gamma);
    

    //Rela_scheme ：是否开启相对运动状态的Riemann solver
    //数值通量
    switch (Rela_scheme) {
        case 0:
//            RS_HLL_XRela(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r,GC);
            break;
        case 1:
            printf("代码未完成\n");
            exit(1);
            break;
        default:
//            RS_Marquina_XRela(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r,GC);
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


static inline void RK3_TimeAd_RPHeat(int AR_scheme, int rows, int cols, int GC, double (*x_d), double (*y)[cols], double dt,double dx, double u_r) {
    int i,j;
    double Conser_U1[rows][cols],Conser_U2[rows][cols];
    double Flux[rows][cols];
    
    for ( i = 0; i < rows; i++){
        for ( j = 0; j < cols; j++){
            Conser_U1[i][j] = Conser_U2[i][j] = 0.0;
            Flux[i][j] = 0.0;
        }
       
    }

 // 第一步计算
    BC_Periodicity(rows,cols,y,0,GC);
    BC_Periodicity(rows,cols,y,1,GC);
//    BC_OutFlow(3,cols,y,0,GC);                       
//    BC_OutFlow(3,cols,y,1,GC);
    Flux_Reconstruction_RP_Heat(AR_scheme,rows,cols,GC,y,Flux,dt,dx,u_r);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC-1; j++){
            Conser_U1[i][j] = y[i][j] - dt*(Flux[i][j]-Flux[i][j-1])/dx;
        }
    }

    //第二步计算
    BC_Periodicity(rows,cols,Conser_U1,0,GC);
    BC_Periodicity(rows,cols,Conser_U1,1,GC);
//    BC_OutFlow(3,cols,Conser_U1,0,GC);                       
//    BC_OutFlow(3,cols,Conser_U1,1,GC);
    Flux_Reconstruction_RP_Heat(AR_scheme,rows,cols,GC,Conser_U1,Flux,dt,dx,u_r);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC-1; j++){
            Conser_U2[i][j] = (3.0/4.0) * y[i][j]  + (1.0/4.0)*(Conser_U1[i][j]- dt*(Flux[i][j]-Flux[i][j-1])/dx);
        }
    }

    //第三步计算
    BC_Periodicity(3,cols,Conser_U2,0,GC);                       
    BC_Periodicity(3,cols,Conser_U2,1,GC);
//    BC_OutFlow(3,cols,Conser_U2,0,GC);                       
//    BC_OutFlow(3,cols,Conser_U2,1,GC);
    Flux_Reconstruction_RP_Heat(AR_scheme,rows,cols,GC,Conser_U2,Flux,dt,dx,u_r);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC-1; j++){
            y[i][j] = (1.0/3.0) * y[i][j]  + (2.0/3.0)*(Conser_U2[i][j]- dt*(Flux[i][j]-Flux[i][j-1])/dx);
        }
    }

    BC_Periodicity(rows,cols,y,0,GC);
    BC_Periodicity(rows,cols,y,1,GC);
}


static inline void RK3_TimeAd(int AR_scheme, int rows, int cols, int GC, double (*x_d), double (*y)[cols], double dt,double dx, double u_r) {
    int i,j;
    double Conser_U1[rows][cols],Conser_U2[rows][cols];
    double Flux[rows][cols];
    
    for ( i = 0; i < rows; i++){
        for ( j = 0; j < cols; j++){
            Conser_U1[i][j] = Conser_U2[i][j] = 0.0;
            Flux[i][j] = 0.0;
        }
       
    }

 // 第一步计算
    BC_Periodicity(rows,cols,y,0,GC);
    BC_Periodicity(rows,cols,y,1,GC);
//    BC_OutFlow(3,cols,y,0,GC);                       
//    BC_OutFlow(3,cols,y,1,GC);
    Flux_Reconstruction_RP(AR_scheme,rows,cols,GC,y,Flux,dt,dx);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC-1; j++){
            Conser_U1[i][j] = y[i][j] - dt*(Flux[i][j]-Flux[i][j-1])/dx;
        }
    }

    //第二步计算
    BC_Periodicity(rows,cols,Conser_U1,0,GC);
    BC_Periodicity(rows,cols,Conser_U1,1,GC);
//    BC_OutFlow(3,cols,Conser_U1,0,GC);                       
//    BC_OutFlow(3,cols,Conser_U1,1,GC);
    Flux_Reconstruction_RP(AR_scheme,rows,cols,GC,Conser_U1,Flux,dt,dx);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC-1; j++){
            Conser_U2[i][j] = (3.0/4.0) * y[i][j]  + (1.0/4.0)*(Conser_U1[i][j]- dt*(Flux[i][j]-Flux[i][j-1])/dx);
        }
    }

    //第三步计算
    BC_Periodicity(3,cols,Conser_U2,0,GC);                       
    BC_Periodicity(3,cols,Conser_U2,1,GC);
//    BC_OutFlow(3,cols,Conser_U2,0,GC);                       
//    BC_OutFlow(3,cols,Conser_U2,1,GC);
    Flux_Reconstruction_RP(AR_scheme,rows,cols,GC,Conser_U2,Flux,dt,dx);
    for ( i = 0; i < rows; i++){
        for ( j = GC; j <= cols-GC-1; j++){
            y[i][j] = (1.0/3.0) * y[i][j]  + (2.0/3.0)*(Conser_U2[i][j]- dt*(Flux[i][j]-Flux[i][j-1])/dx);
        }
    }

    BC_Periodicity(rows,cols,y,0,GC);
    BC_Periodicity(rows,cols,y,1,GC);
}


void OutputConservationErrors_file(int time_step, double Conser_1[3], double Conser_2[3]) {
    // 创建结果目录（如果不存在）
    system("mkdir -p /mnt/d/Desktop/RP_FVM/data");
    
    // 打开文件（追加模式）
    FILE* file = fopen("/mnt/d/Desktop/RP_FVM/data/Conser_Error.dat", "a");  // 修正文件名拼写
    if (file == NULL) {
        printf("无法打开文件 /mnt/d/Desktop/RP_FVM/data/Conser_Error.dat\n");
        return;
    }
    
    // 计算误差（使用更高精度的中间变量）
    volatile double mass_error = Conser_2[0] - Conser_1[0];      // volatile防止编译器优化
    volatile double momentum_error = Conser_2[1] - Conser_1[1];
    volatile double energy_error = Conser_2[2] - Conser_1[2];
    
    // 如果是第一个时间步，写入文件头部
    if (time_step == 0) {
        fprintf(file, "Variables = \t TimeStep \t  Mass_Error  \t Momentum_Error  \t Energy_Error \n");
        fprintf(file, "Zone T=\"Conservation Errors\"\n");
    }
    
    // 写入数据（强制16位小数，使用科学计数法）
    fprintf(file, "%d\t%.16e\t%.16e\t%.16e\n", 
            time_step, 
            (double)mass_error, 
            (double)momentum_error, 
            (double)energy_error);
    
    // 立即刷新缓冲区确保数据写入
    fflush(file);
    
    // 关闭文件
    fclose(file);
    
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
