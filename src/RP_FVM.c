#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>  // 新增，用于目录操作
#include <time.h>
#include "initialize.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "Boundary_Condition.h"
#include "time_advance.h"
#include "scheme.h"
#include "Golbal.h"
#include "Error.h"

#include <omp.h>


// 全局变量存储当前算例
TestCase2D current_test_case = TEST_1D_SHOCKTUBE; // 默认值

// 算例相关函数声明
const char* getTestCaseName(TestCase2D test_case);
const char* getTestCaseShortName(TestCase2D test_case);
void printAvailableTestCases();
void set_current_test_case(TestCase2D test_case);

/*                                            *********                                          */
/*                                            声明子程序                                          */
/*                                            *********                                          */
// 控制计算步数子程序
void StepLoop();
void OutputData_file();
void OutputData_file_2D(bool Con_out, int rows, int cols, int GC, double *mx, double *my, 
                        double (*U)[rows][cols], double (*FU)[rows][cols], 
                        double (*pri)[rows][cols], double now_time, int scheme_type);

// 读取黎曼问题
void Get_RP_Parameters();

//简单的网格代码
void Mesh(int n, double deltax, double *x);
void Mesh_2D();
/*                                            *********                                          */
/*                                              主程序                                            */
/*                                            *********                                          */

//网格参数
const int L_nx = 200;                                                      //网格数量
const int L_ny = 200;                                                      //网格数量      
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

    // 增加计时变量
    clock_t start_time, end_time;
    double total_cpu_time = 0.0;
    double average_time_per_step = 0.0;

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

    // 显示可用算例
    printAvailableTestCases();

    // 选择测试算例
    TestCase2D selected_test = TEST_TAYLOR_GREEN_VORTEX;

    // 设置当前算例（全局变量）
    set_current_test_case(selected_test);

    // 初始化
    Init_Euler_2D(selected_test,var,LNX_ngc,LNY_ngc,GhostCell,pri,U,FU,GU,&Lx,&Ly,&Tmax);

    printf("Read initial conditions successfully!\n");
    printf("Euler equation initialization successful!\n");

    Delta_x = Lx / L_nx;
    Delta_y = Ly / L_ny; 
    double u_r = 0.0;
    
    // mesh
    Mesh_2D(L_nx, L_ny, mesh_x, mesh_y, Delta_x, Delta_y);

    double output_time = Tmax / Control_output;
    double next_output_time = output_time;
    printf("Mesh successfully!\n");
    
    Get_RP_Parameters(argc, argv, &Time_ADM, &scheme);
    
    // 记录程序开始时间
    start_time = clock();
    
    switch (Control_Compution){
        case 0:
            for (Time = 0.0; Time < Tmax; Time = Time+Delta_T) {
                // 记录单步开始时间
                clock_t step_start_time = clock();
                
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
                
                // 记录单步结束时间并累加
                clock_t step_end_time = clock();
                total_cpu_time += ((double)(step_end_time - step_start_time)) / CLOCKS_PER_SEC;
    
                // 检查是否到达输出时间点（使用容差比较浮点数）
                if (fabs(currentTime - next_output_time) < 1e-10 * output_time) {
                    // 计算平均每步耗时
                    if (Ite > 0) {
                        average_time_per_step = total_cpu_time / Ite;
                        printf("Step = %d     Time = %f     Avg time per step = %.6f seconds\n", 
                               Ite, currentTime, average_time_per_step);
                    } else {
                        printf("Step = %d     Time = %f  \n", Ite, currentTime);
                    }
                    
                    printf("Calculation of step %d is completed \n", Ite);
                    
                    Con_to_Pri_2D(var,LNX_ngc,LNY_ngc,pri,U);
                    OutputData_file_2D(Control_Out,LNX_ngc,LNY_ngc,GhostCell,mesh_x,mesh_y,U,FU,pri,currentTime,scheme);
                    // 更新下一个输出时间点
                    next_output_time += output_time;
                }

                Ite++;
                
            }
            break;
        case 1:
            // 迭代步数控制
            Time = 0.0;
            for (int m = 0; m <= 399; m++) {
                // 记录单步开始时间
                clock_t step_start_time = clock();
                
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
                
                // 记录单步结束时间并累加
                clock_t step_end_time = clock();
                total_cpu_time += ((double)(step_end_time - step_start_time)) / CLOCKS_PER_SEC;
                
                Ite++;
                Time = Time + Delta_T;
                
                if (Ite % Control_output == 0 ){
                    // 计算平均每步耗时
                    if (Ite > 0) {
                        average_time_per_step = total_cpu_time / Ite;
                        printf("Step = %d     Time = %f     Avg time per step = %.6f seconds\n", 
                               Ite, Time, average_time_per_step);
                    } else {
                        printf("Step = %d     Time = %f  \n", Ite, Time);
                    }
                    printf("Calculation of step %d is completed \n", Ite);
                }
            }
            break;
        default:
            break;
    }
    
    // 记录程序结束时间
    end_time = clock();
    double total_program_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    // 计算最终的平均每步耗时
    if (Ite > 0) {
        average_time_per_step = total_cpu_time / Ite;
    }

       printf("\n");
    printf("════════════════════════════════════════════════════════════════════\n");
    printf("                   CALCULATION SUMMARY                              \n");
    printf("════════════════════════════════════════════════════════════════════\n");
    printf("\n");
    printf("  ▸ Final Statistics:\n");
    printf("    • Steps Completed:   %d\n", Ite);
    printf("    • Final Time:        %.6f\n", Time);
    printf("\n");
    printf("  ▸ Performance Metrics:\n");
    printf("    • Total Time:        %.3f seconds\n", total_program_time);
    printf("    • Avg Time/Step:     %.6f seconds\n", average_time_per_step);
    if (Ite > 0) {
        printf("    • Steps per Second:  %.2f steps/sec\n", Ite / total_program_time);
    }
    printf("\n");
    printf("════════════════════════════════════════════════════════════════════\n");
    printf("✓ Calculation completed successfully!\n");
    printf("\n");

    // 实现输出最后的结果
    Con_to_Pri_2D(var,LNX_ngc,LNY_ngc,pri,U);
    // 打开文件并输出结果
    Control_Out = true;
    OutputData_file_2D(Control_Out,LNX_ngc,LNY_ngc,GhostCell,mesh_x,mesh_y,U,FU,pri,Tmax,scheme);

    printf("The program has completed its execution.\n");
    free(pri);
    free(U);
    free(FU);
    free(GU);
    printf("Memory Deallocation succesed in Main\n");
    
    return 0;
}





/*                                            *********                                          */
/*                                            算例名称函数实现                                    */
/*                                            *********                                          */

// 获取完整算例名称
const char* getTestCaseName(TestCase2D test_case) {
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: return "1D_Sod_Shocktube";
        case TEST_2D_SHOCKTUBE_CASE1: return "2D_Riemann_Case1";
        case TEST_2D_SHOCKTUBE_CASE2: return "2D_Riemann_Case2";
        case TEST_2D_SHOCKTUBE_CASE3: return "2D_Riemann_Case3";
        case TEST_2D_SHOCKTUBE_CASE4: return "2D_Riemann_Case4";
        case TEST_2D_SHOCKTUBE_CASE5: return "2D_Riemann_Case5";
        case TEST_TAYLOR_GREEN_VORTEX: return "Taylor_Green_Vortex";
        case TEST_GAUSSIAN_PULSE: return "Gaussian_Pulse";
        case TEST_KELVIN_HELMHOLTZ: return "Kelvin_Helmholtz_Smooth";
        case TEST_KELVIN_HELMHOLTZ_SHARP: return "Kelvin_Helmholtz_Sharp";
        case TEST_KELVIN_HELMHOLTZ_VP: return "Kelvin_Helmholtz_VP";
        case TEST_RAYLEIGH_TAYLOR: return "Rayleigh_Taylor";
        case TEST_DOUBLE_MACH_REFLECTION: return "Double_Mach_Reflection";
        case TEST_BACKWARD_STEP: return "Backward_Step";
        case TEST_BLAST_WAVE: return "Blast_Wave";
        case TEST_NOH_PROBLEM: return "Noh_Problem";
        default: return "Unknown_Case";
    }
}

// 获取简短算例名称（用于文件名）
const char* getTestCaseShortName(TestCase2D test_case) {
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: return "Sod";
        case TEST_2D_SHOCKTUBE_CASE1: return "Riemann1";
        case TEST_2D_SHOCKTUBE_CASE2: return "Riemann2";
        case TEST_2D_SHOCKTUBE_CASE3: return "Riemann3";
        case TEST_2D_SHOCKTUBE_CASE4: return "Riemann4";
        case TEST_2D_SHOCKTUBE_CASE5: return "Riemann5";
        case TEST_TAYLOR_GREEN_VORTEX: return "TGV";
        case TEST_GAUSSIAN_PULSE: return "Gaussian";
        case TEST_KELVIN_HELMHOLTZ: return "KH_Smooth";
        case TEST_KELVIN_HELMHOLTZ_SHARP: return "KH_Sharp";
        case TEST_KELVIN_HELMHOLTZ_VP: return "KH_VP";
        case TEST_RAYLEIGH_TAYLOR: return "RT";
        case TEST_DOUBLE_MACH_REFLECTION: return "DMR";
        case TEST_BACKWARD_STEP: return "BackStep";
        case TEST_BLAST_WAVE: return "Blast";
        case TEST_NOH_PROBLEM: return "Noh";
        default: return "Unknown";
    }
}
// 设置当前算例
void set_current_test_case(TestCase2D test_case) {
    current_test_case = test_case;
    printf("Current test case set to: %s\n", getTestCaseName(test_case));
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

    // 默认值设为1
    *Parameter1 = 1;
    *Parameter2 = 1;

    if (argc > 2) {
        // 有两个命令行参数
        *Parameter1 = atoi(argv[1]);
        *Parameter2 = atoi(argv[2]);
        printf("✓ Using Parameter1 %d and Parameter2 %d from command line arguments\n", *Parameter1, *Parameter2);
    } else if (argc > 1) {
        // 只有一个命令行参数
        *Parameter1 = atoi(argv[1]);
        printf("✓ Using Parameter1 %d from command line argument\n", *Parameter1);
        
        // 提示用户输入第二个参数，默认值为1
        printf("Please enter the Parameter2 value (1-ll, default: 1): ");
        char input[100];
        fgets(input, sizeof(input), stdin);
        
        if (strlen(input) == 1) {  // 只输入了回车
            *Parameter2 = 1;
            printf("✓ Using default Parameter2 = 1 (HLL)\n");
        } else {
            *Parameter2 = atoi(input);
            printf("✓ Using Parameter2 %d from user input\n", *Parameter2);
        }
    } else {
        // 没有命令行参数，交互式输入两个参数
        char input1[100], input2[100];
        
        // 输入第一个参数
        printf("Please enter Parameter1 value (1 or 3, default: 1): ");
        fgets(input1, sizeof(input1), stdin);
        
        if (strlen(input1) == 1) {  // 只输入了回车
            *Parameter1 = 1;
        } else {
            *Parameter1 = atoi(input1);
        }
        
        // 输入第二个参数
        printf("Please enter Parameter2 value (1-ll, default: 1): ");
        fgets(input2, sizeof(input2), stdin);
        
        if (strlen(input2) == 1) {  // 只输入了回车
            *Parameter2 = 1;
        } else {
            *Parameter2 = atoi(input2);
        }
        
        printf("✓ Using Parameter1 %d and Parameter2 %d from user input\n", *Parameter1, *Parameter2);
    }

    // 验证输入的有效性，如果无效则使用默认值1
    if (*Parameter1 != 1 && *Parameter1 != 3) {
        printf("⚠ Warning: Parameter1 value %d is invalid, using default 1 (RK1)\n", *Parameter1);
        *Parameter1 = 1;
    }
    
    // 检查Parameter2是否在常见范围内
    if (*Parameter2 < 0 || *Parameter2 > 100) {
        printf("⚠ Warning: Parameter2 value %d is outside recommended range\n", *Parameter2);
        printf("Using default HLL solver (1)\n");
        *Parameter2 = 1;
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

// 更简单的文件名生成，避免缓冲区溢出
// 修改后的OutputData_file_2D函数，包含scheme信息
void OutputData_file_2D(bool Con_out, int rows, int cols, int GC, double *mx, double *my, 
                        double (*U)[rows][cols], double (*FU)[rows][cols], 
                        double (*pri)[rows][cols], double now_time, int scheme_type) {
    printf("Output result(rho, u, v, p, T, U_M, U_E)\n");
    
    // 获取算例名称
    const char* case_name;
    switch(current_test_case) {
        case TEST_NOH_PROBLEM: 
            case_name = "NOH"; 
            break;
        case TEST_BACKWARD_STEP: 
            case_name = "BACKSTEP"; 
            break;
        case TEST_DOUBLE_MACH_REFLECTION: 
            case_name = "DMR"; 
            break;
        case TEST_BLAST_WAVE: 
            case_name = "BLAST"; 
            break;
        case TEST_1D_SHOCKTUBE: 
            case_name = "SOD"; 
            break;
        case TEST_2D_SHOCKTUBE_CASE1: 
            case_name = "RIEMANN1"; 
            break;
        case TEST_2D_SHOCKTUBE_CASE2: 
            case_name = "RIEMANN2"; 
            break;
        case TEST_2D_SHOCKTUBE_CASE3: 
            case_name = "RIEMANN3"; 
            break;
        case TEST_2D_SHOCKTUBE_CASE4: 
            case_name = "RIEMANN4"; 
            break;
        case TEST_2D_SHOCKTUBE_CASE5: 
            case_name = "RIEMANN5"; 
            break;
        case TEST_TAYLOR_GREEN_VORTEX: 
            case_name = "TGV"; 
            break;
        case TEST_GAUSSIAN_PULSE: 
            case_name = "GAUSSIAN"; 
            break;
        case TEST_KELVIN_HELMHOLTZ: 
            case_name = "KH"; 
            break;
        case TEST_KELVIN_HELMHOLTZ_SHARP: 
            case_name = "KH_SHARP"; 
            break;
        case TEST_KELVIN_HELMHOLTZ_VP: 
            case_name = "KH_VP"; 
            break;
        case TEST_RAYLEIGH_TAYLOR: 
            case_name = "RT"; 
            break;
        default: 
            case_name = "UNKNOWN";
    }
    
    // 获取Riemann Solver名称
    const char* scheme_name;
    switch(scheme_type) {
        case 1: scheme_name = "HLL"; break;
        case 2: scheme_name = "HLLC"; break;
        case 3: scheme_name = "Roe"; break;
        case 11: scheme_name = "HLLHC"; break;
        case 22: scheme_name = "HLLCHC"; break;
        case 33: scheme_name = "RoeHC"; break;
        default: scheme_name = "ExactRiemann"; break;
    }
    
    // 生成文件名：算例名称 + scheme类型 + 原来的文件名格式
    char filename[200];
    sprintf(filename, "/mnt/d/Desktop/RP_FVM/data/%s_%s_output_data_%.6f.plt", 
            case_name, scheme_name, now_time);
    
    printf("Writing to file: %s\n", filename);
    
    FILE* file = fopen(filename, "w");  
    if (file == NULL) {
        printf("File opening failed\n");
        printf("---------------Error----------------\n");
        
        // 尝试创建目录
        system("mkdir -p /mnt/d/Desktop/RP_FVM/data");
        
        // 再次尝试打开文件
        file = fopen(filename, "w");
        if (file == NULL) {
            printf("Still failed to open file. Check permissions.\n");
            return;
        }
    }
    
    // Tecplot格式头信息
    fprintf(file, "TITLE = \"2D Fluid Dynamics Data - %s with %s\"\n", case_name, scheme_name);
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
    
    printf("Test Case: %s\n", case_name);
    printf("Riemann Solver: %s\n", scheme_name);
    printf("Time: %.6f\n", now_time);
    printf("File: %s\n", filename);
    
    if (Con_out) {
        printf("Output calculation result successful\n");
    } else {
        printf("The %d th calculation ended at %f\n", Ite, now_time);
    }
    printf("----------------------------------------\n");
}