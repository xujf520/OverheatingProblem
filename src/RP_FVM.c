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
void OutputData_file_2D();

// 读取黎曼问题
int Get_RP_Method(int argc, char *argv[]);
int Get_RP_input(int argc, char *argv[]);
void Get_RP_Parameters(int argc, char *argv[], int *RP_Method, int *scheme);

//简单的网格代码
void Mesh(int n, double deltax, double *x);
void Mesh_2D();
/*                                            *********                                          */
/*                                              主程序                                            */
/*                                            *********                                          */

//网格参数
const int L_nx = 400;                                                      //网格数量
const int L_ny = 100;                                                      //网格数量      
const int var = 4;                                        

int RP_Method;                                                //Riemann Solver的具体方法 
//Riemann Solver
int scheme;                                                  
int main(int argc, char *argv[]) {

    int LNX_ngc = L_nx + 2 * GhostCell; 
    int LNY_ngc = L_ny + 2 * GhostCell; 

    Get_RP_Parameters(argc, argv, &RP_Method, &scheme);

    double t = 0;
    double L = 1.,       Tmax = 0.15;       //计算域参数 
    double Delta_x, Delta_y;
    double Delta_T;

    // 添加声明，2个虚拟网格
    
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


    Delta_x = L / L_nx;
    Delta_y = L / L_ny; 

    int Ite = 0;
//初始条件
    double u_r = 0.0;

    Init_Sod_2D(pri_Ver1,pri_Ver2);                       
//    Init_Shock_Impact(pri_Ver1,pri_Ver2,u_r);                  //激波对撞 
//    Init_DRare(pri_Ver1,pri_Ver2,u_r);  
//   Init_Shock3(pri_Ver1,pri_Ver2);  
    printf("Read initial conditions successfully!\n");
  
    //mesh
//    Mesh(n,Delta_x,x);
    Mesh_2D(L_nx, L_ny, mesh_x, mesh_y, Delta_x, Delta_y);

    printf("Mesh successfully!\n");

    Init_Euler_2D(var,LNX_ngc,LNY_ngc,GhostCell,pri,U,FU,GU,pri_Ver1,pri_Ver2,Delta_x);            //初始化欧拉方程
    printf("Euler equation initialization successful!\n");


    switch (Control_Compution){
        case 0:
            //时间推进：时间一阶和时间二阶格式
            for (t = 0.0; t < Tmax; t = t+Delta_T) {
//               Delta_T = Get_Delta_T(3,N_ngc,U,Delta_x);
                Delta_T = Get_Delta_T_2D(var,LNX_ngc,LNY_ngc,U,Delta_x);
                if (t + Delta_T >= Tmax)
                    Delta_T = Tmax - t ;

                switch (RP_Method) {
                    case 0:
                        //RK1_TimeAd(scheme,3,N_ngc,GhostCell,x,U,FU,Delta_T,Delta_x, u_r);
//                        RK3_TimeAd(scheme,3,N_ngc,GhostCell,x,U,Delta_T,Delta_x, u_r);
                        RK1_TimeAd(scheme,var,LNX_ngc,LNY_ngc,GhostCell,U,Delta_T,Delta_x,Delta_y);
                        break;
                    case 1:
//                        RK1_TimeAd_RPHeat(scheme,3,N_ngc,GhostCell,x,U,Delta_T,Delta_x,u_r);
//                        RK1_TimeAd_RPHeat(scheme,3,N_ngc,GhostCell,x,U,Delta_T,Delta_x,u_r);
//                        RK3_TimeAd_RPHeat(scheme,3,N_ngc,GhostCell,x,U,Delta_T,Delta_x, u_r);
                        break;
                    default:
//                        RK1_TVD_FluxRela(scheme,3,N_ngc,x,U,FU,Delta_T,Delta_x,1,u_r);     
                        break;
                }
                Ite++;
                if (Ite % Control_output == 0 )
                    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t+Delta_T, Ite);
            }
            break;
        case 1:
            //迭代步数控制
            t = 0.0;
            for (int m = 0; m <= 399; m++) {
                Delta_T = 0.1*Delta_x;
                switch (RP_Method) {
                    case 0:
//                        RK1_TimeAd(scheme,3,N_ngc,GhostCell,x,U,Delta_T,Delta_x,u_r);
                        break;
                    case 1:
 //                       RK1_TimeAd_RPHeat(scheme,3,N_ngc,GhostCell,x,U,Delta_T,Delta_x,u_r);
                        break;
                    default:
 //                       RK1_TVD_FluxRela(scheme,3,N_ngc,x,U,FU,Delta_T,Delta_x,1,u_r);     
                        break;
                }
                Ite++;
                t = t+Delta_T;
                if (Ite % Control_output == 0 ){
                    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t, Ite);
                }
            }
            break;
        default:
            break;
    }

    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t, Ite);
    printf("end of calculation!\n");

    //实现输出最后的结果
    Con_to_Pri_2D(var,LNX_ngc,LNY_ngc,pri,U);
    //打开文件并输出结果
    OutputData_file_2D(0,LNX_ngc,LNY_ngc,GhostCell,mesh_x,mesh_y,U,FU,pri,t);

//    Scheme_Error_L1(3,N_ngc,GhostCell,U,Delta_x);
    
    printf("The program has completed its execution.\n");


    free(pri);
    free(U);
    free(FU);
    free(GU);
    printf("Memory Deallocation succesed in Main\n");

    return 0;

}


// 函数定义：同时获取RP_Method和scheme两个参数
// 函数定义：同时获取RP_Method和scheme两个参数
void Get_RP_Parameters(int argc, char *argv[], int *RP_Method, int *scheme) {
    // 先显示所有可选择的内容
    printf("╔═══════════════════════════════════════════════════╗\n");
    printf("║               Available Options                   ║\n");
    printf("╠═══════════════════════════════════════════════════╣\n");
    printf("║ RP Methods:                                       ║\n");
    printf("║   0: Origin                                       ║\n");
    printf("║   1: Heat Conduction                              ║\n");
    printf("║                                                   ║\n");
    printf("║ Riemann Solvers:                                  ║\n");
    printf("║   0: Lax                                          ║\n");
    printf("║   1: Rusanov                                      ║\n");
    printf("║   2: HLL                                          ║\n");
    printf("║   3: HLLC                                         ║\n");
    printf("║   4: Roe                                          ║\n");
    printf("║   5: Marquina                                     ║\n");
    printf("║   6: StegerWarming                                ║\n");
    printf("║   7: VanLeer                                      ║\n");
    printf("║   8: LiouSteffen                                  ║\n");
    printf("║   Other: Exact Riemann                            ║\n");
    printf("╚═══════════════════════════════════════════════════╝\n\n");

    if (argc > 2) {
        // 有两个命令行参数
        *RP_Method = atoi(argv[1]);
        *scheme = atoi(argv[2]);
        printf("✓ Using RP_Method %d and scheme %d from command line arguments\n", *RP_Method, *scheme);
    } else if (argc > 1) {
        // 只有一个命令行参数，提示用户输入另一个
        *RP_Method = atoi(argv[1]);
        printf("✓ Using RP_Method %d from command line argument\n", *RP_Method);
        
        printf("Please enter the scheme value (0-8): ");
        scanf("%d", scheme);
        printf("✓ Using scheme %d from user input\n", *scheme);
    } else {
        // 没有命令行参数，交互式输入两个参数（合并为一行）
        printf("Please enter RP_Method and scheme values (0-1, 0-8): ");
        scanf("%d %d", RP_Method, scheme);
        
        printf("✓ Using RP_Method %d and scheme %d from user input\n", *RP_Method, *scheme);
    }

    // 验证输入的有效性
    if (*RP_Method < 0 || *RP_Method > 1) {
        printf("⚠ Warning: RP_Method value %d is outside recommended range (0-1)\n", *RP_Method);
    }
    
    if (*scheme < 0 || *scheme > 8) {
        printf("⚠ Warning: Scheme value %d is outside recommended range (0-8)\n", *scheme);
        printf("Using Exact Riemann solver as default\n");
    }
    
    printf("\n");
}




//划分网格代码
void Mesh(int n, double deltax, double *x) {
    // 初始化网格节点坐标
    for (int i = 0; i < n; i++) {
        x[i] = deltax * i + 0.5 * deltax;
    }
}

void Mesh_2D(int cols, int depth, double *mesh_x, double *mesh_y, double deltax, double deltay) {

    for (int j = 0; j < cols; j++) 
            mesh_x[j] = deltax * j + 0.5 * deltax;
    
    
    for (int k = 0; k < depth; k++) 
        mesh_y[k] = deltay * k + 0.5 * deltay;
        
}



// 控制计算步数子程序
void StepLoop(double* deltat, double deltax, double CFL, double t, int l, int cols, int vis, int scheme, 
                                            double* x, double (*U)[cols], double (*FU)[cols],double (*pri)[cols],double gamma) {
    for (int k = 0; k < l; k++) {
        //计算格式
        //RK1_TVD(vis, scheme, 3, cols, x, U, FU, 1e-10, deltax, CFL, gamma, k+1);
        //RK_2Roe(vis, scheme, 3, N_ngc, U, FU, deltat, deltax, CFL);

        t += *deltat;
        Con_to_Pri_1D(3, cols, pri, U);
//        OutputData_file(k+1,cols,x,U,FU,pri,t);
    }
}


void OutputData_file(int k, int cols, int GC, double * x, double (*U)[cols], double (*FU)[cols],double (*pri)[cols], double t) {
    // 文件输出
    printf("Output result(rho, u, p, T, U_M, U_E)\n");
    if (k == 0)
    {
        FILE* file = fopen("/mnt/d/Desktop/RP_FVM/data/output_data.dat","w");
        if (file == NULL) {
            
            printf("File opening failed\n");
            printf("---------------Error----------------\n");
            return; // 返回错误代码
        }
        fprintf(file, "variables=x \t rho \t u\t p\t T\t rhou\t rhoE\t rhoE_K\t rhoE_I\t F_rho\t F_rhou\t F_rhoE\n");
        for (int i = GC; i <= cols-GC-1 ; i++) {
            fprintf(file,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",\
                        x[i-GC], pri[0][i], pri[1][i], pri[2][i], pri[2][i]/pri[0][i],\
                        U[1][i], U[2][i], 0.5*pri[1][i]*U[1][i], pri[2][i]/(1.4-1) , \
                        U[1][i], U[1][i] * pri[1][i] + pri[2][i] ,pri[1][i] * (U[2][i] + pri[2][i]));
        }
        fclose(file);
        printf("Output calculation result successful\n");
       // printf("%f\n",t);
    }
    else
    {
        char filename[100];
        sprintf(filename, "/mnt/d/Desktop/RP_FVM/data/output_%d.dat", k);
        FILE* file = fopen(filename, "w");
        if (file == NULL) {
            printf("无法打开文件 %s\n", filename);
            return;
        }
        fprintf(file, "variables =  'x' 'rho' 'u' 'p' 'rhou' \n");
        for (int i = 2; i < cols-2; i++) {
            fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", 
            x[i-2], pri[0][i], pri[1][i], pri[2][i], U[1][i]);
        }
        fclose(file);
        printf("The %d th calculation ended at %f\n", k, t);
    }
}


void OutputData_file_2D(int k, int cols, int depth, int GC, 
                       double (*mx), double (*my), double (*U)[cols][depth], double (*FU)[cols],
                       double (*pri)[cols][depth], double t) {
    
    printf("Output result(rho, u, v, p, T, U_M, U_E)\n");
    
    if (k == 0) {
        FILE* file = fopen("/mnt/d/Desktop/RP_FVM/data/output_data.dat", "w");
        if (file == NULL) {
            printf("File opening failed\n");
            printf("---------------Error----------------\n");
            return;
        }
        
        // Tecplot格式头信息
        fprintf(file, "TITLE = \"2D Fluid Dynamics Data\"\n");
        fprintf(file, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"rhou\", \"rhov\", \"rhoE\"\n");
        
        int output_cols = cols - 2*GC;
        int output_depth = depth - 2*GC;
        
        // 指定ZONE信息 - 关键修复
        fprintf(file, "ZONE T=\"Time=%.6f\"\n", t);
        fprintf(file, "I=%d, J=%d\n", output_depth, output_cols);  // 注意：Tecplot中J是行数，I是列数
        fprintf(file, "DATAPACKING=POINT\n");
        
        // 输出数据 - 注意循环顺序
        for (int j = 0; j < output_cols; j++) {
            for (int kk = 0; kk < output_depth; kk++) {  // 避免变量名冲突
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
        printf("Output calculation result successful\n");
    }
    else {
        char filename[100];
        sprintf(filename, "/mnt/d/Desktop/RP_FVM/data/output_%d.dat", k);
        FILE* file = fopen(filename, "w");
        if (file == NULL) {
            printf("无法打开文件 %s\n", filename);
            return;
        }
        
        // 为后续时间步也添加完整的2D格式
        fprintf(file, "TITLE = \"2D Fluid Dynamics Data - Step %d\"\n", k);
        fprintf(file, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"rhou\", \"rhov\", \"rhoE\"\n");
        
        int output_cols = cols - 2*GC;
        int output_depth = depth - 2*GC;
        
        fprintf(file, "ZONE T=\"Time=%.6f\"\n", t);
        fprintf(file, "I=%d, J=%d\n", output_depth, output_cols);
        fprintf(file, "DATAPACKING=POINT\n");
        
        for (int j = 0; j < output_cols; j++) {
            for (int kk = 0; kk < output_depth; kk++) {
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
        printf("The %d th calculation ended at %f\n", k, t);
    }
}
