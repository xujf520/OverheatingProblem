#ifndef TIME_ADVANCE_H
#define TIME_ADVANCE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
static inline void RK1_TVD(int Recon_Accur, int AR_scheme, int rows, int cols , double (*x_d), double (*y)[cols]\
                                            ,double (*z)[cols],double dt,double dx, double gamma, int k, double u_r) {
    int i,j;
    double Conser_1[3]={0.0}, Conser_2[3]={0.0};


    //施加边界条件，
    BC_OutFlow(3,cols,y,0,2);                       //边界条件说明见具体子程序
//    BC_Reflect(3,cols,y,1,2);
    BC_OutFlow(3,cols,y,1,2);
//    BC_In_rho(3,cols,y,1,2);
//  BC_OutFlow(3,cols,y,1,1);

    //AR_scheme is Approximate Riemann Solver
     //数值通量
    switch (AR_scheme) {
        case 0:
            RS_Lax(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 1:
            RS_Rusanov(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 2:
            RS_HLL(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 3:
            RS_HLLC(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        case 4:
            RS_Roe(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        case 5:
            RS_Marquina(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        case 6:
            RS_StegerWarming(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        case 7:
            RS_VanLeer(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        case 8:
            RS_LiouSteffen(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        default:
            RS_ER(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }

    // 第一步计算
/*                                            *********                                          */
/*                                            输出数值通量                                        */
/*                                            *********                                          */
//    OutputFluxData_file(k, cols, dx, x_d, z);

    Total_Conser(rows,cols,2,y,Conser_1,dx);
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            y[i][j] = y[i][j] - dt*(z[i][j]-z[i][j-1])/dx;
        }
    }
    Total_Conser(rows,cols,2,y,Conser_2,dx);
//    OutputConservationErrors_file(Time_Step,Conser_1,Conser_2);
    Time_Step++;
}


//TVD重构配合不同Riemann Solver
static inline void RK1_TVD_RP_HeatConduction(int Recon_Accur, int AR_scheme, int rows, int cols , double (*x_d), double (*y)[cols]\
                                            ,double (*z)[cols],double dt,double dx, double gamma, int k,double u_r) {
    int i,j;
     //施加边界条件，
    BC_OutFlow(3,cols,y,0,2);                       //边界条件说明见具体子程序
    BC_Reflect(3,cols,y,1,2);

    //AR_scheme is Approximate Riemann Solver
     //数值通量
    switch (AR_scheme) {
        case 0:
            RS_Lax_HeatConduction(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 1:
            RS_Rusanov_HeatConduction(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 2:
            RS_HLL_HeatConduction(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 3:
            RS_HLLC_HeatConduction(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 4:
            RS_Roe_HeatConduction(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 5:
            RS_Marquina(Recon_Accur,3,cols,y,z,dt,dx,gamma);
            break;
        default:
            RS_ER_Heat(Recon_Accur,3,cols,y,z,dt,dx,gamma);
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
    
}



static inline void RK1_TVD_FluxRela(int Recon_Accur, int Rela_scheme, int rows, int cols , double (*x_d), double (*y)[cols]\
                                            ,double (*z)[cols],double dt,double dx, double gamma, int k, double u_r) {
    int i,j;
    //计算时间步长
    //*dt =Get_lamdat(cols,y,dx,CFL,gamma);
    

    //Rela_scheme ：是否开启相对运动状态的Riemann solver
    //数值通量
    switch (Rela_scheme) {
        case 0:
            RS_HLL_XRela(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
            break;
        case 1:
            printf("代码未完成\n");
            exit(1);
            break;
        default:
            RS_Marquina_XRela(Recon_Accur,3,cols,y,z,dt,dx,gamma,u_r);
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
