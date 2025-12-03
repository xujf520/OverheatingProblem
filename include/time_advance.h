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
static inline void RK1_TimeAd(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt, double dx, double dy) {
    int i,j,k;
    double Conser_1[3]={0.0}, Conser_2[3]={0.0};
    double (*Space_Item)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_F)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Source_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));

    // 检查内存分配是否成功
    if (Flux_F == NULL || Flux_G == NULL) {
        fprintf(stderr, "Memory allocation failed in RK1_TimeAd\n");
        // 释放已分配的内存
        free(Space_Item);
        free(Flux_F);
        free(Flux_G);
        free(Source_G);
        return;
    }
  
    //初始化参数
    for ( k = 0; k < var; k++){
        for ( i = 0; i < rows; i++){
            for ( j = 0; j < cols; j++){
                Space_Item[k][i][j] = 0.0;
                Flux_F[k][i][j] = 0.0;
                Flux_G[k][i][j] = 0.0;
                Source_G[k][i][j] = 0.0;
            }
        }
    }

    if (Source)
        Source_Gravity(var,rows,cols,GC,y,Source_G); 

   // 第一步计算
    // 施加边界条件
//    Boundary_Conditions(var, rows, cols, y, GC); 
    BC_BackwardStep_2D(var, rows, cols, y, GC,Time); 
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,y,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    
    for ( i = GC; i <= rows-GC-1; i++)
        for ( j = GC; j <= cols-GC-1; j++)
            for ( k = 0; k < var; k++)
                y[k][i][j] = y[k][i][j] + dt * Space_Item[k][i][j];
            
        
    free(Space_Item);
    free(Flux_F);
    free(Flux_G);
    free(Source_G);
    

}



/*static inline void RK3_TimeAd(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt,double dx, double dy) {
    int i,j,k;

    double (*Conser_U1)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U2)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_F)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Source_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    // 检查内存分配是否成功
    if (Flux_F == NULL || Flux_G == NULL) {
        fprintf(stderr, "Memory allocation failed in RK3_TimeAd\n");
        // 释放已分配的内存
        free(Conser_U1);
        free(Conser_U2);
        free(Flux_F);
        free(Flux_G);
        free(Source_G);
        return;
    }    
    //初始化参数
    for ( k = 0; k < var; k++){
        for ( i = 0; i < rows; i++){
            for ( j = 0; j < cols; j++){
                Conser_U1[k][i][j] = Conser_U2[k][i][j] = 0.0;
                Flux_F[k][i][j] = Flux_G[k][i][j] = 0.0;
                Source_G[k][i][j] = 0.0;
            }
        }
    }

    
    // 第一步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,y,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, y, GC,Time); 
    Boundary_Conditions(var, rows, cols, y, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,y,Flux_F,Flux_G,dt,dx,dy);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                Conser_U1[k][i][j] = y[k][i][j] - dt * (Flux_F[i][j][k]-Flux_F[i][j-1][k])/dx \
                                            - dt * (Flux_G[i][j][k]-Flux_G[i][j][k-1])/dy + dt * Source_G[k][i][j]; 


    //第二步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U1,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, Conser_U1, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U1, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U1,Flux_F,Flux_G,dt,dx,dy);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                Conser_U2[k][i][j] = (3.0/4.0) * y[k][i][j] +  (1.0/4.0) * (Conser_U1[i][j][k] \
                                            - dt * (Flux_F[i][j][k]-Flux_F[i][j-1][k])/dx \
                                            - dt * (Flux_G[i][j][k]-Flux_G[i][j][k-1])/dy + dt * Source_G[k][i][j]);
    
    
    //第三步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U2,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, Conser_U2, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U2, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U2,Flux_F,Flux_G,dt,dx,dy);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                y[k][i][j] = (1.0/2.0) * y[k][i][j] +  (2.0/3.0) * (Conser_U2[i][j][k] \
                                                    - dt * (Flux_F[i][j][k]-Flux_F[i][j-1][k])/dx \
                                                    - dt * (Flux_G[i][j][k]-Flux_G[i][j][k-1])/dy + dt * Source_G[k][i][j]);



    free(Conser_U1);
    free(Conser_U2);
    free(Flux_F);
    free(Flux_G);
    free(Source_G); 

}*/
static inline void RK3_TimeAd(int AR_scheme, int var, int rows,int cols, int GC, double (*y)[rows][cols], double dt,double dx, double dy) {
    int i,j,k;

    double (*Space_Item)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U1)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conser_U2)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_F)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Flux_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Source_G)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    // 检查内存分配是否成功
    if (Flux_F == NULL || Flux_G == NULL) {
        fprintf(stderr, "Memory allocation failed in RK3_TimeAd\n");
        // 释放已分配的内存
        free(Space_Item);
        free(Conser_U1);
        free(Conser_U2);
        free(Flux_F);
        free(Flux_G);
        free(Source_G);
        return;
    }    
    //初始化参数
    for ( k = 0; k < var; k++){
        for ( i = 0; i < rows; i++){
            for ( j = 0; j < cols; j++){
                Conser_U1[k][i][j] = Conser_U2[k][i][j] = 0.0;
                Flux_F[k][i][j] = Flux_G[k][i][j] = 0.0;
                Space_Item[k][i][j] = 0.0;
                Source_G[k][i][j] = 0.0;
            }
        }
    }

    // 第一步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,y,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, y, GC,Time); 
    Boundary_Conditions(var, rows, cols, y, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,y,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                Conser_U1[k][i][j] = y[k][i][j] + dt * Space_Item[k][i][j];


    //第二步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U1,Source_G); 
//    BC_DoubleMach_2D(var, rows, cols, Conser_U1, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U1, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U1,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                Conser_U2[k][i][j] = (3.0/4.0) * y[k][i][j] +  (1.0/4.0) * (Conser_U1[k][i][j] + dt * Space_Item[k][i][j]);
    
    
    //第三步计算
    if (Source)
        Source_Gravity(var,rows,cols,GC,Conser_U2,Source_G); 
    // 施加边界条件
//    BC_DoubleMach_2D(var, rows, cols, Conser_U2, GC,Time); 
    Boundary_Conditions(var, rows, cols, Conser_U2, GC);
    Flux_Reconstruction_RP(AR_scheme,var,rows,cols,GC,Conser_U2,Flux_F,Flux_G,dt,dx,dy);
    Space_Discrete_Item(var,rows,cols,GC,dx,dy,Flux_F,Flux_G,Source_G,Space_Item);

    for ( k = 0; k < var; k++)
        for ( i = GC; i <= rows-GC-1; i++)
            for ( j = GC; j <= cols-GC-1; j++)
                y[k][i][j] = (1.0/3.0) * y[k][i][j] +  (2.0/3.0) * (Conser_U2[k][i][j] + dt * Space_Item[k][i][j]);



    free(Space_Item);
    free(Conser_U1);
    free(Conser_U2);
    free(Flux_F);
    free(Flux_G);
    free(Source_G); 

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





void OutputFluxData_file(int k, int rows, double dx, double* x,double (*z)[rows]) {
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
    for (int i = 2; i < rows-2 ; i++) {
        fprintf(file,"%.8f\t%.8f\t%.8f\t%.8f\n",x[i-2]+0.5*dx, z[0][i],z[1][i], z[2][i]);
    }
    fclose(file);
    printf("numerical flux ：The %d th calculation ended\n", k);
}

#endif 
