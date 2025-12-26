#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>  // For Directory Operations
#include <omp.h>
#include <time.h>
#include "mesh.h"
#include "initialize.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "Boundary_Condition.h"
#include "time_advance.h"
#include "scheme.h"
#include "Golbal.h"
#include "Error.h"



// Unified High Precision Timer
static inline double get_wall_time(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
#endif
}

// Global Variable Stores The Current Case
TestCase2D current_test_case = TEST_1D_SHOCKTUBE; // 默认值

// Function Declarations Related To Examples
const char* getTestCaseName(TestCase2D test_case);
const char* getTestCaseShortName(TestCase2D test_case);
void ReadControlFile(const char* filename);
void set_material_parameters(double gamma, bool source, double gravity, double cfl);
void set_current_test_case(TestCase2D test_case);

/*                                            *********                                          */
/*                                         Declare Subroutine                                    */
/*                                            *********                                          */
// Subroutine For Controlling Calculation Steps
void StepLoop();
void OutputData_file();
void OutputData_file_2D(bool Con_out, int rows, int cols, int GC, double *mx, double *my, 
                        double (*U)[rows][cols], double (*FU)[rows][cols], 
                        double (*pri)[rows][cols], double now_time, int scheme_type);
void OutputData_file_2D_BS(bool Con_out, int rows, int cols, int GC, double *mx, double *my, 
                        double (*U)[rows][cols], double (*FU)[rows][cols], 
                        double (*pri)[rows][cols], double now_time, int scheme_type);

// Reading The Riemann Problem
void Get_RP_Parameters();

//Simple Grid Code
void Mesh(int n, double deltax, double *x);
void Mesh_2D();
/*                                            *********                                          */
/*                                        Partial Parameters                                      */
/*                                            *********                                          */
// Control Parameter Structure
typedef struct {
    TestCase2D test_case;
    int L_nx;
    int L_ny;
    int Time_ADM;
    int scheme;
    double M_gamma;
    bool Source;
    double Gravity;
    double CFL;
    double Tmax;
    int Control_Compution;
    int Control_output;
    int Recon_Accur;      
    bool Characteriz;     
    char output_dir[256];
} ControlParams;
// Global Control Parameter Variable
ControlParams ctrl_params;

BoundaryConfig bc_config = {BC_OUTFLOW, BC_OUTFLOW, BC_OUTFLOW, BC_OUTFLOW};

/*                                            *********                                          */
/*                                         Main Program                                            */
/*                                            *********                                          */

int main(int argc, char *argv[]) {
   //read Control File
    if (argc > 1) {
        ReadControlFile(argv[1]);
    } else {
        printf("Using default parameters (no control file specified)\n");
        // Set Default Control Parameters
        ctrl_params.test_case = TEST_BLAST_WAVE;
        ctrl_params.L_nx = 100;
        ctrl_params.L_ny = 100;
        ctrl_params.Time_ADM = 3;
        ctrl_params.scheme = 1;
        ctrl_params.M_gamma = 1.4;
        ctrl_params.Source = false;
        ctrl_params.Gravity = 1.0;
        ctrl_params.CFL = 0.4;
        ctrl_params.Tmax = 0.1;
        ctrl_params.Control_Compution = 0;
        ctrl_params.Control_output = 2;
        strcpy(ctrl_params.output_dir, "/mnt/d/Desktop/RP_FVM/data");
    }
    // Set Global Parameters According To The Control File
    set_current_test_case(ctrl_params.test_case);
    
    int LNX_ngc = ctrl_params.L_nx + 2 * GhostCell; 
    int LNY_ngc = ctrl_params.L_ny + 2 * GhostCell; 
    
    double Lx, Ly;
    double Tmax = ctrl_params.Tmax;       
    double Delta_x, Delta_y;
    double Delta_T;

    // Timing Variable
    double program_start_time, program_end_time;
    double total_wall_time = 0.0;
    double average_time_per_step = 0.0;

    double mesh_x[ctrl_params.L_nx], mesh_y[ctrl_params.L_ny];
    double pri_Ver1[4], pri_Ver2[4];
    double (*pri)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    double (*U)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    double (*FU)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    double (*GU)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
    // Check If The Memory Allocation Was Successful
    if (pri == NULL || U == NULL || FU == NULL || GU == NULL) {
        fprintf(stderr, "Memory allocation failed in Main\n");
        free(pri); free(U); free(FU); free(GU);
        return 1;
    } 
    // Set Material Parameters
    set_material_parameters(ctrl_params.M_gamma, ctrl_params.Source, ctrl_params.Gravity, ctrl_params.CFL);
    // Initialize
    Init_Euler_2D(ctrl_params.test_case, var, LNX_ngc, LNY_ngc, GhostCell, pri, U, FU, GU, &Lx, &Ly, &Tmax);

    printf("Read initial conditions successfully!\n");
    printf("Euler equation initialization successful!\n");

    Delta_x = Lx / ctrl_params.L_nx;
    Delta_y = Ly / ctrl_params.L_ny;
    
    // mesh
    Mesh_2D(ctrl_params.L_nx, ctrl_params.L_ny, mesh_x, mesh_y, Delta_x, Delta_y);

    double output_time = Tmax / ctrl_params.Control_output;
    double next_output_time = output_time;
    printf("Mesh successfully!\n");

    // Record Program Start Time
    program_start_time = omp_get_wtime();
    
    switch (ctrl_params.Control_Compution) {
        case 0:
            for (Time = 0.0; Time < Tmax; Time = Time + Delta_T) {
                // Record Single Step Start Time
                double step_start_time = omp_get_wtime();
                Delta_T = Get_Delta_T_2D(var, LNX_ngc, LNY_ngc, U, Delta_x, Delta_y);
                
                // Ensure It Does Not Exceed The Next Output Time Point Or Tmax
                if (Time + Delta_T > next_output_time) {
                    Delta_T = next_output_time - Time;
                }
                if (Time + Delta_T >= Tmax) {
                    Delta_T = Tmax - Time;
                }

                switch (ctrl_params.Time_ADM) {
                    case 1:
                        RK1_TimeAd_Unified(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y, ctrl_params.test_case);
                        break;
                    case 3:
                        RK3_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y);
                        //RK3_TimeAd_Unified(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y, ctrl_params.test_case);
                        break;
                    default:
                        fprintf(stderr, "Error: Invalid time advancement method\n");
                        break;
                }

                double currentTime = Time + Delta_T;
                
                // 记录单步结束时间并累加
                double step_end_time = omp_get_wtime();
                double step_wall_time = step_end_time - step_start_time;
                total_wall_time += step_wall_time;

                // 检查是否到达输出时间点
                if (fabs(currentTime - next_output_time) < 1e-10 * output_time) {
                    // 计算平均每步耗时
                    if (Ite > 0) {
                        average_time_per_step = total_wall_time / Ite;
                        printf("Step = %d     Time = %.6f     Wall time per step = %.6f seconds", 
                               Ite, currentTime, average_time_per_step);
                        printf(" (Threads: %d)\n", omp_get_max_threads());
                    } else {
                        printf("Step = %d     Time = %.6f  \n", Ite, currentTime);
                    }
                    
                    printf("Calculation of step %d is completed \n", Ite);
                    
                    Con_to_Pri_2D(var, LNX_ngc, LNY_ngc, pri, U);
                    OutputData_file_2D(Control_Out, LNX_ngc, LNY_ngc, GhostCell, 
                                       mesh_x, mesh_y, U, FU, pri, currentTime, ctrl_params.scheme);
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
                double step_start_time = omp_get_wtime();
                
                Delta_T = 0.1 * Delta_x;
                switch (ctrl_params.Time_ADM) {
                    case 1:
                        RK1_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, 
                                   U, Delta_T, Delta_x, Delta_y);
                        break;
                    case 3:
                        RK3_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, 
                                   U, Delta_T, Delta_x, Delta_y);
                        break;
                    default:
                        fprintf(stderr, "Error: Invalid time advancement method\n");
                        break;
                }
                
                // 记录单步结束时间并累加
                double step_end_time = omp_get_wtime();
                double step_wall_time = step_end_time - step_start_time;
                total_wall_time += step_wall_time;
                
                Ite++;
                Time = Time + Delta_T;
                
                if (Ite % ctrl_params.Control_output == 0) {
                    // 计算平均每步耗时
                    if (Ite > 0) {
                        average_time_per_step = total_wall_time / Ite;
                        printf("Step = %d     Time = %.6f     Wall time per step = %.6f seconds", 
                               Ite, Time, average_time_per_step);
                        printf(" (Threads: %d)\n", omp_get_max_threads());
                    } else {
                        printf("Step = %d     Time = %.6f  \n", Ite, Time);
                    }
                    printf("Calculation of step %d is completed \n", Ite);
                }
            }
            break;
            
        default:
            fprintf(stderr, "Error: Invalid control computation mode\n");
            break;
    }
    
    // 记录程序结束时间
    program_end_time = omp_get_wtime();
    double total_program_time = program_end_time - program_start_time;
    
    // 计算最终的平均每步耗时
    if (Ite > 0) {
        average_time_per_step = total_wall_time / Ite;
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
    printf("    • Total Wall Time:   %.3f seconds\n", total_program_time);
    printf("    • Avg Time/Step:     %.6f seconds\n", average_time_per_step);
    if (Ite > 0) {
        printf("    • Steps per Second:  %.2f steps/sec\n", Ite / total_program_time);
    }
    printf("    • Threads Used:      %d\n", omp_get_max_threads());
    printf("\n");
    printf("════════════════════════════════════════════════════════════════════\n");
    printf("✓ Calculation completed successfully!\n");
    printf("\n");

    // 实现输出最后的结果
    Con_to_Pri_2D(var, LNX_ngc, LNY_ngc, pri, U);
    // 打开文件并输出结果
    Control_Out = true;
    OutputData_file_2D(Control_Out, LNX_ngc, LNY_ngc, GhostCell, 
                       mesh_x, mesh_y, U, FU, pri, Tmax, ctrl_params.scheme);

    printf("The program has completed its execution.\n");
    free(pri);
    free(U);
    free(FU);
    free(GU);
    printf("Memory Deallocation succesed in Main\n");
    
    return 0;
}



void ReadControlFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open control file %s\n", filename);
        exit(1);
    }
    
    char line[256];
    char key[100];
    char value[150];
    char test_case_name[100];
    
    // 设置默认值
    ctrl_params.test_case = TEST_RAYLEIGH_TAYLOR;
    ctrl_params.L_nx = 100;
    ctrl_params.L_ny = 400;
    ctrl_params.Time_ADM = 3;
    ctrl_params.scheme = 1;
    ctrl_params.M_gamma = 1.667;
    ctrl_params.Source = true;
    ctrl_params.Gravity = 1.0;
    ctrl_params.CFL = 0.4;
    ctrl_params.Tmax = 1.0;
    ctrl_params.Control_Compution = 0;
    ctrl_params.Control_output = 2;
    ctrl_params.Recon_Accur = 5;        // 新增：默认5阶WENO重构
    ctrl_params.Characteriz = false;    // 新增：默认不开启特征重构
    strcpy(ctrl_params.output_dir, "/mnt/d/Desktop/RP_FVM/data");
    
    while (fgets(line, sizeof(line), file)) {
        // 跳过注释行和空行
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }
        
        // 移除行尾的换行符
        line[strcspn(line, "\n")] = 0;
        line[strcspn(line, "\r")] = 0;
        
        // 解析键值对
        if (sscanf(line, "%99[^:]: %149[^\n]", key, value) == 2) {
            // 去除键值两端的空格
            char *key_trim = key;
            char *value_trim = value;
            while (*key_trim == ' ') key_trim++;
            while (*key_trim && key_trim[strlen(key_trim)-1] == ' ') 
                key_trim[strlen(key_trim)-1] = 0;
            while (*value_trim == ' ') value_trim++;
            while (*value_trim && value_trim[strlen(value_trim)-1] == ' ') 
                value_trim[strlen(value_trim)-1] = 0;
            
            if (strcmp(key_trim, "TestCase") == 0) {
                if (strcmp(value_trim, "Sod_Shocktube") == 0) ctrl_params.test_case = TEST_1D_SHOCKTUBE;
                else if (strcmp(value_trim, "Riemann_Case1") == 0) ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE1;
                else if (strcmp(value_trim, "Riemann_Case2") == 0) ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE2;
                else if (strcmp(value_trim, "Riemann_Case3") == 0) ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE3;
                else if (strcmp(value_trim, "Riemann_Case4") == 0) ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE4;
                else if (strcmp(value_trim, "Riemann_Case5") == 0) ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE5;
                else if (strcmp(value_trim, "Taylor_Green_Vortex") == 0) ctrl_params.test_case = TEST_TAYLOR_GREEN_VORTEX;
                else if (strcmp(value_trim, "Gaussian_Pulse") == 0) ctrl_params.test_case = TEST_GAUSSIAN_PULSE;
                else if (strcmp(value_trim, "Kelvin_Helmholtz") == 0) ctrl_params.test_case = TEST_KELVIN_HELMHOLTZ;
                else if (strcmp(value_trim, "Rayleigh_Taylor") == 0) ctrl_params.test_case = TEST_RAYLEIGH_TAYLOR;
                else if (strcmp(value_trim, "Double_Mach_Reflection") == 0) ctrl_params.test_case = TEST_DOUBLE_MACH_REFLECTION;
                else if (strcmp(value_trim, "Backward_Step") == 0) ctrl_params.test_case = TEST_BACKWARD_STEP;
                else if (strcmp(value_trim, "Blast_Wave") == 0) ctrl_params.test_case = TEST_BLAST_WAVE;
                else if (strcmp(value_trim, "Noh_Problem") == 0) ctrl_params.test_case = TEST_NOH_PROBLEM;
                else {
                    fprintf(stderr, "Warning: Unknown test case '%s', using default Rayleigh_Taylor\n", value_trim);
                }
            }
            else if (strcmp(key_trim, "L_nx") == 0) ctrl_params.L_nx = atoi(value_trim);
            else if (strcmp(key_trim, "L_ny") == 0) ctrl_params.L_ny = atoi(value_trim);
            else if (strcmp(key_trim, "Time_ADM") == 0) ctrl_params.Time_ADM = atoi(value_trim);
            else if (strcmp(key_trim, "scheme") == 0) ctrl_params.scheme = atoi(value_trim);
            else if (strcmp(key_trim, "M_gamma") == 0) ctrl_params.M_gamma = atof(value_trim);
            else if (strcmp(key_trim, "Source") == 0) {
                if (strcmp(value_trim, "true") == 0 || strcmp(value_trim, "1") == 0) 
                    ctrl_params.Source = true;
                else 
                    ctrl_params.Source = false;
            }
            else if (strcmp(key_trim, "Gravity") == 0) ctrl_params.Gravity = atof(value_trim);
            else if (strcmp(key_trim, "CFL") == 0) ctrl_params.CFL = atof(value_trim);
            else if (strcmp(key_trim, "Tmax") == 0) ctrl_params.Tmax = atof(value_trim);
            else if (strcmp(key_trim, "Control_Compution") == 0) ctrl_params.Control_Compution = atoi(value_trim);
            else if (strcmp(key_trim, "Control_output") == 0) ctrl_params.Control_output = atoi(value_trim);
            else if (strcmp(key_trim, "output_dir") == 0) strcpy(ctrl_params.output_dir, value_trim);
            // 新增：重构精度参数
            else if (strcmp(key_trim, "Recon_Accur") == 0) {
                int recon_val = atoi(value_trim);
                // 验证重构精度值的有效性
                if (recon_val == 1 || recon_val == 2 || recon_val == 3 || recon_val == 5) {
                    ctrl_params.Recon_Accur = recon_val;
                } else {
                    fprintf(stderr, "Warning: Invalid Recon_Accur value %d. Valid values are 1, 2, 3, 5. Using default 5.\n", recon_val);
                    ctrl_params.Recon_Accur = 5;
                }
            }
            // 新增：特征重构参数
            else if (strcmp(key_trim, "Characteriz") == 0) {
                if (strcmp(value_trim, "true") == 0 || strcmp(value_trim, "1") == 0) 
                    ctrl_params.Characteriz = true;
                else if (strcmp(value_trim, "false") == 0 || strcmp(value_trim, "0") == 0)
                    ctrl_params.Characteriz = false;
                else {
                    fprintf(stderr, "Warning: Invalid Characteriz value '%s'. Using default false.\n", value_trim);
                    ctrl_params.Characteriz = false;
                }
            }
        }
    }
    
    fclose(file);
    
    // 打印读取的参数
    printf("╔═══════════════════════════════════════════════════╗\n");
    printf("║            Control Parameters Loaded              ║\n");
    printf("╠═══════════════════════════════════════════════════╣\n");
    printf("║ Test Case:         %-30s ║\n", getTestCaseName(ctrl_params.test_case));
    printf("║ Grid Size (nx×ny): %-5d × %-5d               ║\n", ctrl_params.L_nx, ctrl_params.L_ny);
    printf("║ Time Integration:  %-30s ║\n", ctrl_params.Time_ADM == 1 ? "RK1" : "RK3");
    printf("║ Riemann Solver:    %-30s ║\n", 
           ctrl_params.scheme == 1 ? "HLL" : 
           ctrl_params.scheme == 2 ? "HLLC" : 
           ctrl_params.scheme == 3 ? "Roe" : 
           ctrl_params.scheme == 11 ? "HLLHC" : 
           ctrl_params.scheme == 22 ? "HLLCHC" : 
           ctrl_params.scheme == 33 ? "RoeHC" : "ExactRiemann");
    printf("║ Gamma (γ):         %-30.3f ║\n", ctrl_params.M_gamma);
    printf("║ Gravity Source:    %-30s ║\n", ctrl_params.Source ? "ON" : "OFF");
    printf("║ Gravity Constant:  %-30.3f ║\n", ctrl_params.Gravity);
    printf("║ CFL Number:        %-30.3f ║\n", ctrl_params.CFL);
    printf("║ Final Time (Tmax): %-30.3f ║\n", ctrl_params.Tmax);
    printf("║ Output Steps:      %-30d ║\n", ctrl_params.Control_output);
    printf("║ Output Directory:  %-30s ║\n", ctrl_params.output_dir);
    // 新增：重构参数显示
    printf("║ Reconstruction:    %-30s ║\n", 
           ctrl_params.Recon_Accur == 1 ? "0th Order (Constant)" :
           ctrl_params.Recon_Accur == 2 ? "2nd Order TVD" :
           ctrl_params.Recon_Accur == 3 ? "3rd Order WENO" :
           ctrl_params.Recon_Accur == 5 ? "5th Order WENO" : "Unknown");
    printf("║ Characteristic:    %-30s ║\n", ctrl_params.Characteriz ? "Characteristic" : "Primitive");
    printf("╚═══════════════════════════════════════════════════╝\n\n");
}

// 读取控制文件的函数
void set_material_parameters(double gamma, bool source, double gravity, double cfl) {

    M_gamma = gamma;
    Source = source;
    Gravity = gravity;

    CFL = cfl;
    
    printf("Material parameters set:\n");
    printf("  Gamma: %.3f\n", M_gamma);
    printf("  Gravity source: %s\n", Source ? "ON" : "OFF");
    if (Source) {
        printf("  Gravity constant: %.3f\n", Gravity);
    }
    printf("  CFL number: %.3f\n", CFL);
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
        case TEST_KELVIN_HELMHOLTZ: return "Kelvin_Helmholtz";
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






// 辅助函数：获取格式名称
const char* getSchemeName(int scheme_type) {
    switch(scheme_type) {
        case 1: return "HLL";
        case 2: return "HLLC";
        case 3: return "Roe";
        case 11: return "HLLHC";
        case 22: return "HLLCHC";
        case 33: return "RoeHC";
        default: return "ExactRiemann";
    }
}

// 辅助函数：判断点是否在空白区域（后台阶算例）
// 空白区域：x在[0.6, 3.0]且y在[0.0, 0.2]
int isBlankRegion(double x, double y) {
    return (x >= 0.6 && x <= 3.0 && y >= 0.0 && y <= 0.2) ? 0 : 1;
}

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
    
    // 生成文件名：使用 .dat 扩展名表示ASCII格式
    char filename[512];
    sprintf(filename, "%s/%s_%s_output_data_%.6f.plt", 
            ctrl_params.output_dir, 
            getTestCaseShortName(ctrl_params.test_case),
            scheme_name, 
            now_time);
    
    printf("Writing to file: %s\n", filename);
    
    // 创建目录（如果不存在）
    char mkdir_cmd[1024];
    sprintf(mkdir_cmd, "mkdir -p %s", ctrl_params.output_dir);
    system(mkdir_cmd);
    
    FILE* file = fopen(filename, "w");  
    if (file == NULL) {
        printf("Error: Failed to open file: %s\n", filename);
        printf("Check permissions and path.\n");
        return;
    }
    
    // Tecplot格式头信息 - 修正变量列表，添加IBLANK
    fprintf(file, "TITLE = \"2D Fluid Dynamics Data - %s with %s\"\n", case_name, scheme_name);
    fprintf(file, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"rhou\", \"rhov\", \"rhoE\", \"IBLANK\"\n");
    
    int output_rows = rows - 2*GC;
    int output_cols = cols - 2*GC;
    
    // 检查是否为后台阶算例
    if (current_test_case == TEST_BACKWARD_STEP) {
        printf("Processing BACKWARD_STEP case with IBLANK blanking...\n");
        
        // 对于后台阶算例，输出整个区域，使用IBLANK标记空白区域
        // 空白区域规则：x在[0.6, 3.0]且y在[0.0, 0.2]的是空区域
        
        // 修正ZONE语法：所有参数在一行，使用POINT格式
        fprintf(file, "ZONE T=\"Backward Step: Time=%.6f\", I=%d, J=%d, DATAPACKING=POINT\n", 
                now_time, output_cols, output_rows);
        
        // 输出整个区域的数据
        for (int i = 0; i < output_rows; i++) {
            for (int j = 0; j < output_cols; j++) {
                int actual_i = i + GC;
                int actual_j = j + GC;
                
                // 判断是否为空白区域
                int iblank = 1; // 默认1表示有效区域
                if (mx[i] >= 0.6 && mx[i] <= 3.0 && my[j] >= 0.0 && my[j] <= 0.2) {
                    iblank = 0; // 0表示空白区域
                }
                
                double temperature = pri[3][actual_i][actual_j] / pri[0][actual_i][actual_j];
                
                // 输出数据，包括IBLANK值
                fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\n",
                        mx[i], my[j], 
                        pri[0][actual_i][actual_j], 
                        pri[1][actual_i][actual_j], 
                        pri[2][actual_i][actual_j], 
                        pri[3][actual_i][actual_j],
                        temperature,
                        U[1][actual_i][actual_j], 
                        U[2][actual_i][actual_j], 
                        U[3][actual_i][actual_j],
                        iblank);
            }
        }
        
        printf("Backward Step region output with IBLANK: %d x %d points\n", output_cols, output_rows);
        printf("IBLANK=0 for blanked region (x=[0.6,3.0], y=[0.0,0.2]), IBLANK=1 for active region.\n");
        
    } else {
        // 其他算例保持原样输出，IBLANK全部设为1
        fprintf(file, "ZONE T=\"Time=%.6f\", I=%d, J=%d, DATAPACKING=POINT\n", 
                now_time, output_cols, output_rows);
        
        // 输出数据
        for (int i = 0; i < output_rows; i++) {
            for (int j = 0; j < output_cols; j++) {
                int actual_i = i + GC;
                int actual_j = j + GC;
                
                double temperature = pri[3][actual_i][actual_j] / pri[0][actual_i][actual_j];
                
                // IBLANK=1 表示所有区域有效
                fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\n",
                        mx[i], my[j], 
                        pri[0][actual_i][actual_j], 
                        pri[1][actual_i][actual_j], 
                        pri[2][actual_i][actual_j], 
                        pri[3][actual_i][actual_j],
                        temperature,
                        U[1][actual_i][actual_j], 
                        U[2][actual_i][actual_j], 
                        U[3][actual_i][actual_j],
                        1); // IBLANK=1
            }
        }
    }
    
    fclose(file);
    
    printf("Test Case: %s\n", case_name);
    printf("Riemann Solver: %s\n", scheme_name);
    printf("Time: %.6f\n", now_time);
    printf("File: %s\n", filename);
    printf("Grid size: %d x %d (excluding ghost cells)\n", output_cols, output_rows);
    
    if (Con_out) {
        printf("Output calculation result successful\n");
    } else {
        printf("The %d th calculation ended at %f\n", Ite, now_time);
    }
    printf("----------------------------------------\n");
}