#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include "initialize.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "Boundary_Condition.h"
#include "time_advance.h"
#include "scheme.h"
#include "Golbal.h"
#include "Error.h"

/*                                            *********                                          */
/*                                            声明子程序                                          */
/*                                            *********                                          */
// 控制计算步数子程序
void StepLoop();
void OutputData_file();
void OutputData_file_2D(bool Con_out, int rows, int cols, int GC, double *mx, double *my,
                            double (*U)[rows][cols], double (*FU)[rows][cols], double (*pri)[rows][cols],double now_time);

// 读取黎曼问题
void Get_RP_Parameters();

//简单的网格代码
void Mesh(int n, double deltax, double *x);
void Mesh_2D();
/*                                            *********                                          */
/*                                              主程序                                            */
/*                                            *********                                          */

//网格参数
const int L_nx = 100;                                                      //网格数量
const int L_ny = 400;                                                      //网格数量      
const int var = 4;                                        

int Time_ADM;                                                //Riemann Solver的具体方法 
//Riemann Solver
int scheme;
//边界条件
BoundaryConfig bc_config = {BC_OUTFLOW, BC_OUTFLOW, BC_OUTFLOW, BC_OUTFLOW};

int main(int argc, char *argv[]) {
    int LNX_ngc = L_nx + 2 * GhostCell; 
    int LNY_ngc = L_ny + 2 * GhostCell; 
    
    double Lx,Ly;
    double Tmax;       
    double Delta_x, Delta_y;
    double Delta_T;

    //
    double mesh_x[L_nx], mesh_y[L_ny];
    double pri_Ver1[4], pri_Ver2[4];
    double (*pri)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    double (*U)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    double (*FU)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    double (*GU)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    // 检查内存分配是否成功
    if (pri == NULL || U == NULL || FU == NULL|| GU == NULL) {
        fprintf(stderr, "Memory allocation failed in Main\n");
        // 释放已分配的内存
        free(pri);
        free(U);
        free(FU);
        free(GU);
    } 

    Init_Euler_2D(var,LNX_ngc,LNY_ngc,GhostCell,pri,U,FU,GU,&Lx,&Ly,&Tmax);            //初始化欧拉方程
    printf("Read initial conditions successfully!\n");
    printf("Euler equation initialization successful!\n");

    Delta_x = Lx / L_nx;
    Delta_y = Ly / L_ny; 
    double u_r = 0.0;
    //mesh
    Mesh_2D(L_nx, L_ny, mesh_x, mesh_y, Delta_x, Delta_y);


    double output_time = Tmax / Control_output;
    double next_output_time = output_time;
    printf("Mesh successfully!\n");
    Get_RP_Parameters(argc, argv, &Time_ADM, &scheme);
    switch (Control_Compution){
        case 0:
            for (Time = 0.0; Time < Tmax; Time = Time+Delta_T) {
                Delta_T = Get_Delta_T_2D(var,LNX_ngc,LNY_ngc,U,Delta_x,Delta_y);
                // 确保不会超过下一个输出时间点或Tmax
                if (Time + Delta_T > next_output_time) {
                    Delta_T = next_output_time - Time;
                }
                if (Time + Delta_T >= Tmax) {
                    Delta_T = Tmax - Time;
                }

                switch (Time_ADM) {
                    case 1:
                        RK1_TimeAd(scheme,var,LNX_ngc,LNY_ngc,GhostCell,U,Delta_T,Delta_x,Delta_y);
                        break;
                    case 3:
                        RK3_TimeAd(scheme,var,LNX_ngc,LNY_ngc,GhostCell,U,Delta_T,Delta_x,Delta_y);
                        break;
                    default:
                        //格式     
                        break;
                }
    
                double currentTime = Time + Delta_T;
    
                // 检查是否到达输出时间点（使用容差比较浮点数）
                if (fabs(currentTime - next_output_time) < 1e-10 * output_time) {
                    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, currentTime, Ite);
                    Con_to_Pri_2D(var,LNX_ngc,LNY_ngc,pri,U);
                    OutputData_file_2D(Control_Out,LNX_ngc,LNY_ngc,GhostCell,mesh_x,mesh_y,U,FU,pri,currentTime);
                    // 更新下一个输出时间点
                    next_output_time += output_time;
                }

                Ite++;
                
            }
            break;
        case 1:
            //迭代步数控制
            Time = 0.0;
            for (int m = 0; m <= 399; m++) {
                Delta_T = 0.1*Delta_x;
                switch (Time_ADM) {
                    case 1:
                        RK1_TimeAd(scheme,var,LNX_ngc,LNY_ngc,GhostCell,U,Delta_T,Delta_x,Delta_y);
                        break;
                    case 3:
                        RK3_TimeAd(scheme,var,LNX_ngc,LNY_ngc,GhostCell,U,Delta_T,Delta_x,Delta_y);
                        break;
                    default:
                        //格式     
                        break;
                }
                Ite++;
                Time = Time+Delta_T;
                if (Ite % Control_output == 0 ){
                    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, Time, Ite);
                }
            }
            break;
        default:
            break;
    }

    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, Time, Ite);
    printf("end of calculation!\n");

    //实现输出最后的结果
    Con_to_Pri_2D(var,LNX_ngc,LNY_ngc,pri,U);
    //打开文件并输出结果
    Control_Out = true;
    OutputData_file_2D(Control_Out,LNX_ngc,LNY_ngc,GhostCell,mesh_x,mesh_y,U,FU,pri,Tmax);

    printf("The program has completed its execution.\n");
    free(pri);
    free(U);
    free(FU);
    free(GU);
    printf("Memory Deallocation succesed in Main\n");
    return 0;

}

// 函数定义：同时获取RP_Method和scheme两个参数
void Get_RP_Parameters(int argc, char *argv[], int *Parameter1, int *Parameter2) {
    // 先显示所有可选择的内容
    printf("╔═══════════════════════════════════════════════════╗\n");
    printf("║               Available Options                   ║\n");
    printf("╠═══════════════════════════════════════════════════╣\n");
    printf("║ Time Advance Method:                              ║\n");
    printf("║   1: RK1                                          ║\n");
    printf("║   3: RK3                                          ║\n");
    printf("║                                                   ║\n");
    printf("║ Riemann Solvers:                                  ║\n");
    printf("║   1: HLL                                          ║\n");
    printf("║   2: HLLC                                         ║\n");
    printf("║   3: Roe                                          ║\n");
    printf("║   11: HLL with Heat Conduction    --HLLHC         ║\n");
    printf("║   22: HLLC with Heat Conduction   --HLLCHC        ║\n");
    printf("║   33: Roe with Heat Conduction    --RoeHC         ║\n");
    printf("║   Other: Exact Riemann                            ║\n");
    printf("╚═══════════════════════════════════════════════════╝\n\n");

    if (argc > 2) {
        // 有两个命令行参数
        *Parameter1 = atoi(argv[1]);
        *Parameter2 = atoi(argv[2]);
        printf("✓ Using Parameter1 %d and Parameter2 %d from command line arguments\n", *Parameter1, *Parameter2);
    } else if (argc > 1) {
        // 只有一个命令行参数，提示用户输入另一个
        *Parameter1 = atoi(argv[1]);
        printf("✓ Using Parameter1 %d from command line argument\n", *Parameter1);
        
        printf("Please enter the Parameter2 value (1-ll): ");
        scanf("%d", Parameter2);
        printf("✓ Using Parameter2 %d from user input\n", *Parameter2);
    } else {
        // 没有命令行参数，交互式输入两个参数（合并为一行）
        printf("Please enter Parameter1 and Parameter2 values (1-3, 1-ll): ");
        scanf("%d %d", Parameter1, Parameter2);
        
        printf("✓ Using Parameter1 %d and Parameter2 %d from user input\n", *Parameter1, *Parameter2);
    }

    // 验证输入的有效性
    if (*Parameter1 < 0 || *Parameter1 > 10) {
        printf("⚠ Warning: Parameter1 value %d is outside recommended range (1-3)\n", *Parameter1);
    }
    
    if (*Parameter2 < 0 || *Parameter2 > 100) {
        printf("⚠ Warning: Parameter2 value %d is outside recommended range (1-ll)\n", *Parameter2);
        printf("Using Exact Riemann solver as default\n");
    }
    
    printf("\n");
}



void Mesh_2D(int rows, int cols, double *mesh_x, double *mesh_y, double deltax, double deltay) {

    for (int j = 0; j < rows; j++) 
            mesh_x[j] = deltax * j + 0.5 * deltax;
    
    
    for (int k = 0; k < cols; k++) 
        mesh_y[k] = deltay * k + 0.5 * deltay;
        
}



// 控制计算步数子程序
void StepLoop(double* deltat, double deltax, double CFL, double t, int l, int rows, int vis, int scheme, 
                                            double* x, double (*U)[rows], double (*FU)[rows],double (*pri)[rows],double gamma) {
    for (int k = 0; k < l; k++) {
        //计算格式
        //RK1_TVD(vis, scheme, 3, rows, x, U, FU, 1e-10, deltax, CFL, gamma, k+1);
        //RK_2Roe(vis, scheme, 3, N_ngc, U, FU, deltat, deltax, CFL);

        t += *deltat;
//        Con_to_Pri_2D(var,LNX_ngc,LNY_ngc,pri,U);
//        OutputData_file(k+1,rows,x,U,FU,pri,t);
    }
}


void OutputData_file_2D(bool Con_out, int rows, int cols, int GC, double *mx, double *my, 
                            double (*U)[rows][cols], double (*FU)[rows][cols], double (*pri)[rows][cols],double now_time) {
    printf("Output result(rho, u, v, p, T, U_M, U_E)\n");
    
    
    // 修复：正确格式化文件名
    char filename[100];
    sprintf(filename, "/mnt/d/Desktop/RP_FVM/data/output_data_%.6f.plt", now_time);
    
    FILE* file = fopen(filename, "w");  
    if (file == NULL) {
        printf("File opening failed\n");
        printf("---------------Error----------------\n");
        return;
    }
    
    // Tecplot格式头信息
    fprintf(file, "TITLE = \"2D Fluid Dynamics Data\"\n");
    fprintf(file, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"rhou\", \"rhov\", \"rhoE\"\n");
    
    int output_rows = rows - 2*GC;
    int output_cols = cols - 2*GC;
    
    // 指定ZONE信息
    fprintf(file, "ZONE T=\"Time=%.6f\"\n", now_time);
    fprintf(file, "I=%d, J=%d\n", output_cols, output_rows);
    fprintf(file, "DATAPACKING=POINT\n");
    
    // 输出数据
    for (int j = 0; j < output_rows; j++) {
        for (int kk = 0; kk < output_cols; kk++) {
            int actual_j = j + GC;
            int actual_kk = kk + GC;
            
            double temperature = pri[3][actual_j][actual_kk] / pri[0][actual_j][actual_kk];
            
            fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
                    mx[j], my[kk], 
                    pri[0][actual_j][actual_kk], 
                    pri[1][actual_j][actual_kk], 
                    pri[2][actual_j][actual_kk], 
                    pri[3][actual_j][actual_kk],
                    temperature,
                    U[1][actual_j][actual_kk], 
                    U[2][actual_j][actual_kk], 
                    U[3][actual_j][actual_kk]);
        }
    }
    fclose(file);
    if (Con_out)
        printf("Output calculation result successful\n");
    else
        printf("The %d th calculation ended at %f\n", Ite, now_time);
    
        
}