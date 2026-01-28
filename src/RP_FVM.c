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
void Precision_test();
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
        ctrl_params.test_case = TEST_BACKWARD_STEP;
        ctrl_params.L_nx = 300;
        ctrl_params.L_ny = 100;
        ctrl_params.Time_ADM = 1;
        ctrl_params.scheme = 10; 
        ctrl_params.M_gamma = 1.4;
        ctrl_params.Source = false;
        ctrl_params.Gravity = 1.0;
        ctrl_params.CFL = 0.4;
        ctrl_params.Tmax = 0.1;
        ctrl_params.Control_Compution = 0;
        ctrl_params.Control_output = 4;
        ctrl_params.Recon_Accur = 1;        // 添加默认值
        ctrl_params.Characteriz = true;    // 添加默认值
        strcpy(ctrl_params.output_dir, "/mnt/d/Desktop/RP_FVM/data");
    }
    // Set Global Parameters According To The Control File
    set_current_test_case(ctrl_params.test_case);

    //精度测试：
    if(ctrl_params.test_case == TEST_PRECISION){
        Precision_test();
        return 0;
    }
    
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
        {
            // Declare variables inside the case block with braces
            int last_percent = -1; // Record last percentage to avoid duplicate printing
            double progress_prev = 0.0;
            char spinner[4] = {'|', '/', '-', '\\'};
            int spinner_index = 0;
            
            for (Time = 0.0; Time < Tmax; Time = Time + Delta_T) {
                // Record start time of current step
                double step_start_time = omp_get_wtime();
                Delta_T = Get_Delta_T_2D(var, LNX_ngc, LNY_ngc, U, Delta_x, Delta_y);
                //Delta_T = Get_Delta_T_2D_X(var, LNX_ngc, LNY_ngc, U, Delta_x, Delta_y);
                
                // Ensure we don't exceed next output time point or Tmax
                if (Time + Delta_T > next_output_time) {
                    Delta_T = next_output_time - Time;
                }
                if (Time + Delta_T >= Tmax) {
                    Delta_T = Tmax - Time;
                }

                switch (ctrl_params.Time_ADM) {
                    case 1:
                        RK1_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y);
                        //RK1_TimeAd_Unified(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y, ctrl_params.test_case);
                        break;
                    case 2:
                        RK2_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y);
                        //RK2_TimeAd_Unified(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y, ctrl_params.test_case);
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
                
                // Record end time of current step and accumulate
                double step_end_time = omp_get_wtime();
                double step_wall_time = step_end_time - step_start_time;
                total_wall_time += step_wall_time;
                
                // Calculate progress percentage
                double progress = currentTime / Tmax * 100.0;
                
                // Simple animated progress bar - update each step
                int bar_width = 30;
                int pos = (int)(bar_width * progress / 100.0);
                
                // Update spinner animation
                spinner_index = (spinner_index + 1) % 4;
                
                // Display progress bar (overwrites same line using carriage return)
                printf("\r[");
                for (int i = 0; i < bar_width; ++i) {
                    if (i < pos) printf("█");
                    else if (i == pos) printf("%c", spinner[spinner_index]);
                    else printf(" ");
                }
                printf("] %.2f%% (Step: %d, Time: %.4f)", progress, Ite, currentTime);
                fflush(stdout); // Force flush output buffer
                
                // Check if we've reached output time point
                if (fabs(currentTime - next_output_time) < 1e-10 * output_time) {
                    // Calculate average time per step
                    if (Ite > 0) {
                        average_time_per_step = total_wall_time / Ite;
                        
                        // Print detailed information on new line
                        printf("\n");
                        printf("Step %d completed | Time: %.6f | Avg time/step: %.6f s\n", 
                            Ite, currentTime, average_time_per_step);
                        
                        // Save data
                        Con_to_Pri_2D(var, LNX_ngc, LNY_ngc, pri, U);
                        OutputData_file_2D(Control_Out, LNX_ngc, LNY_ngc, GhostCell, 
                                        mesh_x, mesh_y, U, FU, pri, currentTime, ctrl_params.scheme);
                        printf("Data saved.\n");
                        
                        // Reset progress bar display
                        printf("\n"); // New line for next progress bar
                    } else {
                        printf("\nStep = %d     Time = %.6f  \n", Ite, currentTime);
                        printf("Calculation of step %d is completed \n", Ite);
                    }
                    
                    // Update next output time point
                    next_output_time += output_time;
                }

                Ite++;
            }
            
            // Display after simulation completion
            printf("\r[");
            for (int i = 0; i < 30; ++i) printf("█");
            printf("] 100.00%% (Simulation completed!)\n");
            break;
        }
        
        // 迭代步数控制
        case 1:
           
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


/*==============================================================================
 * Function: ReadControlFile
 * Description: Reads simulation control parameters from a configuration file
 *              and sets up the simulation accordingly.
 * Parameters:
 *   - filename: Path to the control file
 *============================================================================*/
void ReadControlFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open control file %s\n", filename);
        exit(1);
    }
    
    char line[256];
    char key[100];
    char value[150];
    
    // Set default values for all control parameters
    ctrl_params.test_case = TEST_RAYLEIGH_TAYLOR;
    ctrl_params.L_nx = 100;
    ctrl_params.L_ny = 400;
    ctrl_params.Time_ADM = 3;               // Default: RK3 time integration
    ctrl_params.scheme = 1;                 // Default: HLL Riemann solver
    ctrl_params.M_gamma = 1.4;              // Default: Air (γ = 1.4)
    ctrl_params.Source = false;             // Default: No source terms
    ctrl_params.Gravity = 1.0;              // Default gravity constant
    ctrl_params.CFL = 0.4;                  // Default CFL number
    ctrl_params.Tmax = 1.0;                 // Default maximum simulation time
    ctrl_params.Control_Compution = 0;      // Default: Time-controlled simulation
    ctrl_params.Control_output = 2;         // Default output interval
    ctrl_params.Recon_Accur = 2;            // Default: 2nd order TVD reconstruction
    ctrl_params.Characteriz = false;        // Default: Primitive variable reconstruction
    strcpy(ctrl_params.output_dir, "/mnt/d/Desktop/RP_FVM/data");  // Default output directory
    
    // Read and parse control file line by line
    while (fgets(line, sizeof(line), file)) {
        // Skip comment lines and empty lines
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }
        
        // Remove trailing newline characters
        line[strcspn(line, "\n")] = 0;
        line[strcspn(line, "\r")] = 0;
        
        // Parse key-value pairs (format: "key: value")
        if (sscanf(line, "%99[^:]: %149[^\n]", key, value) == 2) {
            // Trim whitespace from key and value
            char *key_trim = key;
            char *value_trim = value;
            while (*key_trim == ' ') key_trim++;
            while (*key_trim && key_trim[strlen(key_trim)-1] == ' ') 
                key_trim[strlen(key_trim)-1] = 0;
            while (*value_trim == ' ') value_trim++;
            while (*value_trim && value_trim[strlen(value_trim)-1] == ' ') 
                value_trim[strlen(value_trim)-1] = 0;
            
            // Process recognized parameters
            if (strcmp(key_trim, "TestCase") == 0) {
                // Map test case string to enumeration value
                if (strcmp(value_trim, "Precision_Test") == 0) 
                    ctrl_params.test_case = TEST_PRECISION;
                else if (strcmp(value_trim, "1D_Sod_Shocktube") == 0) 
                    ctrl_params.test_case = TEST_1D_SHOCKTUBE;
                else if (strcmp(value_trim, "1D_Contact_Wave") == 0) 
                    ctrl_params.test_case = TEST_1D_CONTACTWAVE;
                else if (strcmp(value_trim, "1D_Shock_Impact") == 0) 
                    ctrl_params.test_case = TEST_1D_SHOCKIMPACT;
                else if (strcmp(value_trim, "1D_Impact_Wall") == 0) 
                    ctrl_params.test_case = TEST_1D_IMPACTWALL;
                else if (strcmp(value_trim, "1D_Double_Rarefaction") == 0) 
                    ctrl_params.test_case = TEST_1D_DOUBLERARE;
                else if (strcmp(value_trim, "1D_Noh_Problem") == 0) 
                    ctrl_params.test_case = TEST_1D_NOHPROBLEM;
                else if (strcmp(value_trim, "2D_Riemann_Case1") == 0) 
                    ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE1;
                else if (strcmp(value_trim, "2D_Riemann_Case2") == 0) 
                    ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE2;
                else if (strcmp(value_trim, "2D_Riemann_Case3") == 0) 
                    ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE3;
                else if (strcmp(value_trim, "2D_Riemann_Case4") == 0) 
                    ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE4;
                else if (strcmp(value_trim, "2D_Riemann_Case5") == 0) 
                    ctrl_params.test_case = TEST_2D_SHOCKTUBE_CASE5;
                else if (strcmp(value_trim, "Taylor_Green_Vortex") == 0) 
                    ctrl_params.test_case = TEST_TAYLOR_GREEN_VORTEX;
                else if (strcmp(value_trim, "Gaussian_Pulse") == 0) 
                    ctrl_params.test_case = TEST_GAUSSIAN_PULSE;
                else if (strcmp(value_trim, "Kelvin_Helmholtz") == 0) 
                    ctrl_params.test_case = TEST_KELVIN_HELMHOLTZ;
                else if (strcmp(value_trim, "Rayleigh_Taylor") == 0) 
                    ctrl_params.test_case = TEST_RAYLEIGH_TAYLOR;
                else if (strcmp(value_trim, "Double_Mach_Reflection") == 0) 
                    ctrl_params.test_case = TEST_DOUBLE_MACH_REFLECTION;
                else if (strcmp(value_trim, "Odd_Even_Decoupling") == 0) 
                    ctrl_params.test_case = TEST_OddEven_Decoupling;
                else if (strcmp(value_trim, "Backward_Step") == 0) 
                    ctrl_params.test_case = TEST_BACKWARD_STEP;
                else if (strcmp(value_trim, "Blast_Wave") == 0) 
                    ctrl_params.test_case = TEST_BLAST_WAVE;
                else if (strcmp(value_trim, "2D_Noh_Problem") == 0) 
                    ctrl_params.test_case = TEST_NOH_PROBLEM;
                else {
                    fprintf(stderr, "Warning: Unknown test case '%s', using default Rayleigh_Taylor\n", value_trim);
                }
            }
            else if (strcmp(key_trim, "L_nx") == 0) 
                ctrl_params.L_nx = atoi(value_trim);
            else if (strcmp(key_trim, "L_ny") == 0) 
                ctrl_params.L_ny = atoi(value_trim);
            else if (strcmp(key_trim, "Time_ADM") == 0) 
                ctrl_params.Time_ADM = atoi(value_trim);
            else if (strcmp(key_trim, "scheme") == 0) 
                ctrl_params.scheme = atoi(value_trim);
            else if (strcmp(key_trim, "M_gamma") == 0) 
                ctrl_params.M_gamma = atof(value_trim);
            else if (strcmp(key_trim, "Source") == 0) {
                if (strcmp(value_trim, "true") == 0 || strcmp(value_trim, "1") == 0) 
                    ctrl_params.Source = true;
                else 
                    ctrl_params.Source = false;
            }
            else if (strcmp(key_trim, "Gravity") == 0) 
                ctrl_params.Gravity = atof(value_trim);
            else if (strcmp(key_trim, "CFL") == 0) 
                ctrl_params.CFL = atof(value_trim);
            else if (strcmp(key_trim, "Tmax") == 0) 
                ctrl_params.Tmax = atof(value_trim);
            else if (strcmp(key_trim, "Control_Compution") == 0) 
                ctrl_params.Control_Compution = atoi(value_trim);
            else if (strcmp(key_trim, "Control_output") == 0) 
                ctrl_params.Control_output = atoi(value_trim);
            else if (strcmp(key_trim, "output_dir") == 0) 
                strcpy(ctrl_params.output_dir, value_trim);
            else if (strcmp(key_trim, "Recon_Accur") == 0) {
                int recon_val = atoi(value_trim);
                // Validate reconstruction accuracy value
                if (recon_val == 1 || recon_val == 2 || recon_val == 3 || recon_val == 5) {
                    ctrl_params.Recon_Accur = recon_val;
                } else {
                    fprintf(stderr, "Warning: Invalid Recon_Accur value %d. Valid values are 1, 2, 3, 5. Using default 2.\n", recon_val);
                    ctrl_params.Recon_Accur = 2;
                }
            }
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
            else {
                fprintf(stderr, "Warning: Unknown parameter '%s' in control file\n", key_trim);
            }
        }
    }
    
    fclose(file);
    

    // Display loaded control parameters in formatted output
    printf("╔═══════════════════════════════════════════════════╗\n");
    printf("║            Control Parameters Loaded              ║\n");
    printf("╠═══════════════════════════════════════════════════╣\n");
    printf("║ Test Case:         %-30s ║\n", getTestCaseName(ctrl_params.test_case));
    printf("║ Grid Size (nx×ny): %-5d × %-5d               ║\n", ctrl_params.L_nx, ctrl_params.L_ny);

    // Display time integration method
    printf("║ Time Integration:  %-30s ║\n", 
        ctrl_params.Time_ADM == 1 ? "RK1 (Euler)" : 
        (ctrl_params.Time_ADM == 2 ? "RK2 (Midpoint)" : 
        (ctrl_params.Time_ADM == 3 ? "RK3 (TVD)" : "Unknown")));
        
    // Display Riemann solver type
    const char* scheme_name;
    switch (ctrl_params.scheme) {
        case 1: scheme_name = "HLL"; break;
        case 2: scheme_name = "HLLC"; break;
        case 3: scheme_name = "Roe"; break;
        case 11: scheme_name = "HLL with Heat Conduction"; break;
        case 22: scheme_name = "HLLC with Heat Conduction"; break;
        case 33: scheme_name = "Roe with Heat Conduction"; break;
        default: scheme_name = "Exact Riemann"; break;
    }
    printf("║ Riemann Solver:    %-30s ║\n", scheme_name);
    
    printf("║ Gamma (γ):         %-30.3f ║\n", ctrl_params.M_gamma);
    printf("║ Gravity Source:    %-30s ║\n", ctrl_params.Source ? "ON" : "OFF");
    printf("║ Gravity Constant:  %-30.3f ║\n", ctrl_params.Gravity);
    printf("║ CFL Number:        %-30.3f ║\n", ctrl_params.CFL);
    printf("║ Final Time (Tmax): %-30.3f ║\n", ctrl_params.Tmax);
    printf("║ Control Mode:      %-30s ║\n", 
           ctrl_params.Control_Compution == 0 ? "Time Control" : "Step Control");
    printf("║ Output Interval:   %-30d ║\n", ctrl_params.Control_output);
    printf("║ Output Directory:  %-30s ║\n", ctrl_params.output_dir);
    
    // Display reconstruction method
    const char* recon_name;
    switch (ctrl_params.Recon_Accur) {
        case 1: recon_name = "1st Order (Constant)"; break;
        case 2: recon_name = "2nd Order TVD"; break;
        case 3: recon_name = "3rd Order WENO"; break;
        case 5: recon_name = "5th Order WENO"; break;
        default: recon_name = "Unknown"; break;
    }
    printf("║ Reconstruction:    %-30s ║\n", recon_name);
    
    printf("║ Variable Type:     %-30s ║\n", 
           ctrl_params.Characteriz ? "Characteristic" : "Primitive");
    printf("╚═══════════════════════════════════════════════════╝\n\n");
}

/*==============================================================================
 * Function: set_material_parameters
 * Description: Sets global material parameters for the simulation
 * Parameters:
 *   - gamma: Specific heat ratio (γ)
 *   - source: Enable/disable source terms
 *   - gravity: Gravity constant value
 *   - cfl: CFL stability number
 *============================================================================*/
void set_material_parameters(double gamma, bool source, double gravity, double cfl) {
    M_gamma = gamma;
    Source = source;
    Gravity = gravity;
    CFL = cfl;
}


// Set current test case and log information
void set_current_test_case(TestCase2D test_case) {
    current_test_case = test_case;
    printf("================================================\n");
    printf("Test case configured:\n");
    printf("  Full Name:  %s\n", getTestCaseName(test_case));
    printf("  Short Name: %s\n", getTestCaseShortName(test_case));
    printf("  Description: %s\n", getTestCaseDescription(test_case));
    printf("================================================\n");
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
    printf("Output result(rho, u, v, p, T, Entropy, U_M, U_E)\n");
    
   // Get Case Name Systematic Abbreviation
    const char* case_name;
    switch(current_test_case) {
        case TEST_PRECISION:           case_name = "Error"; break;
        case TEST_1D_SHOCKTUBE:        case_name = "1D_SOD"; break;
        case TEST_1D_CONTACTWAVE:      case_name = "1D_CONTACT"; break;
        case TEST_1D_SHOCKIMPACT:      case_name = "1D_SHOCKIMPACT"; break;
        case TEST_1D_IMPACTWALL:       case_name = "1D_IMPACTWALL"; break;
        case TEST_1D_DOUBLERARE:       case_name = "1D_DOUBLERARE"; break;
        case TEST_1D_NOHPROBLEM:       case_name = "1D_NOH"; break;
        case TEST_2D_SHOCKTUBE_CASE1:  case_name = "2D_RIEMANN1"; break;
        case TEST_2D_SHOCKTUBE_CASE2:  case_name = "2D_RIEMANN2"; break;
        case TEST_2D_SHOCKTUBE_CASE3:  case_name = "2D_RIEMANN3"; break;
        case TEST_2D_SHOCKTUBE_CASE4:  case_name = "2D_RIEMANN4"; break;
        case TEST_2D_SHOCKTUBE_CASE5:  case_name = "2D_RIEMANN5"; break;
        case TEST_TAYLOR_GREEN_VORTEX: case_name = "TGV"; break;
        case TEST_GAUSSIAN_PULSE:      case_name = "GAUSSIAN"; break;
        case TEST_KELVIN_HELMHOLTZ:    case_name = "KH"; break;
        case TEST_RAYLEIGH_TAYLOR:     case_name = "RT"; break;
        case TEST_DOUBLE_MACH_REFLECTION: case_name = "DMR"; break;
        case TEST_BLAST_WAVE:          case_name = "BLASTWAVE"; break;
        case TEST_NOH_PROBLEM:         case_name = "2D_NOH"; break;
        case TEST_OddEven_Decoupling:  case_name = "OEDC"; break;
        case TEST_BACKWARD_STEP:       case_name = "BACKSTEP"; break;
        default:                       case_name = "UNKNOWN"; break;
    }
    
    // Get Riemann Solver Name
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
    
    int output_rows = rows - 2*GC;
    int output_cols = cols - 2*GC;
    
    // 根据维数选择不同的文件名和输出格式
    char filename[512];
    if (output_cols == 1) {
        // 一维剖面数据 - 使用.dat扩展名
        sprintf(filename, "%s/%s_%s_profile_%.6f.dat", 
                ctrl_params.output_dir, 
                getTestCaseShortName(ctrl_params.test_case),
                scheme_name, 
                now_time);
    } else {
        // 二维数据 - 使用.plt扩展名
        sprintf(filename, "%s/%s_%s_output_data_%.6f.plt", 
                ctrl_params.output_dir, 
                getTestCaseShortName(ctrl_params.test_case),
                scheme_name, 
                now_time);
    }
    
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
    
    if (output_cols == 1) {
        // 一维剖面输出 - 优化为Origin格式
        
        // Origin格式：使用注释行开头，然后是表头
        fprintf(file, "! 1D Profile Data - %s with %s\n", case_name, scheme_name);
        fprintf(file, "! Time = %.6f\n", now_time);
        fprintf(file, "! Number of points = %d\n", output_rows);
        // 表头行 - Origin可以识别简单的表头
        // 使用空格分隔的表头，Origin可以自动识别为列名
        fprintf(file, "X\tDensity\tU\tV\tPressure\tTemperature\tEntropy\tMomentumX\tMomentumY\tEnergy\n");
        
        // 一维剖面数据输出（每行对应一个x位置）
        for (int i = 0; i < output_rows; i++) {
            int actual_i = i + GC;
            int actual_j = GC;  // 对于一维情况，只取中间的y位置
            
            double temperature = pri[3][actual_i][actual_j] / pri[0][actual_i][actual_j];
            double entropy = pri[3][actual_i][actual_j] / pow(pri[0][actual_i][actual_j], M_gamma);
            
            // 输出数据，使用制表符分隔，Origin可以正确识别
            // 格式：x坐标 + 8个物理量
            fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
                    mx[i],  // x坐标
                    pri[0][actual_i][actual_j],  // 密度 rho
                    pri[1][actual_i][actual_j],  // x方向速度 u
                    pri[2][actual_i][actual_j],  // y方向速度 v
                    pri[3][actual_i][actual_j],  // 压力 p
                    temperature,                 // 温度 T
                    entropy,                     // 熵
                    U[1][actual_i][actual_j],    // x方向动量 rhou
                    U[2][actual_i][actual_j],    // y方向动量 rhov
                    U[3][actual_i][actual_j]);   // 总能量 rhoE
        }
        
        printf("1D profile data output successful (Origin compatible format)\n");
        
    } else {
        // 二维输出 - 保持原有Tecplot格式
        // Tecplot格式头信息 - 修正变量列表，添加IBLANK
        fprintf(file, "TITLE = \"2D Fluid Dynamics Data - %s with %s\"\n", case_name, scheme_name);
        fprintf(file, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"Entropy\", \"rhou\", \"rhov\", \"rhoE\", \"IBLANK\"\n");
        
        // 检查是否为后台阶算例
        if (current_test_case == TEST_BACKWARD_STEP) {
            printf("Processing BACKWARD_STEP case with IBLANK blanking...\n");
            
            // 对于后台阶算例，输出整个区域，使用IBLANK标记空白区域
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
                    double entropy = pri[3][actual_i][actual_j] / pow(pri[0][actual_i][actual_j], M_gamma);
                    
                    // 输出数据，包括Entropy和IBLANK值
                    fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\n",
                            mx[i], my[j], 
                            pri[0][actual_i][actual_j], 
                            pri[1][actual_i][actual_j], 
                            pri[2][actual_i][actual_j], 
                            pri[3][actual_i][actual_j],
                            temperature,
                            entropy,  // 新增的绝热指数（熵）
                            U[1][actual_i][actual_j], 
                            U[2][actual_i][actual_j], 
                            U[3][actual_i][actual_j],
                            iblank);
                }
            }
            
            printf("Backward Step region output with IBLANK and Entropy\n");
            
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
                    double entropy = pri[3][actual_i][actual_j] / pow(pri[0][actual_i][actual_j], M_gamma);
                    
                    // IBLANK=1 表示所有区域有效
                    fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\n",
                            mx[i], my[j], 
                            pri[0][actual_i][actual_j], 
                            pri[1][actual_i][actual_j], 
                            pri[2][actual_i][actual_j], 
                            pri[3][actual_i][actual_j],
                            temperature,
                            entropy,  // 新增的绝热指数（熵）
                            U[1][actual_i][actual_j], 
                            U[2][actual_i][actual_j], 
                            U[3][actual_i][actual_j],
                            1); // IBLANK=1
                }
            }
        }
    }
    
    fclose(file);
    
    printf("Test Case: %s\n", case_name);
    printf("Riemann Solver: %s\n", scheme_name);
    printf("Time: %.6f\n", now_time);
    printf("File: %s\n", filename);
    
    if (output_cols == 1) {
        printf("1D profile data with %d points (x-direction)\n", output_rows);
        printf("Format optimized for Origin software\n");
        printf("Columns: X, Density, U, V, Pressure, Temperature, Entropy, MomentumX, MomentumY, Energy\n");
    } else {
        printf("Grid size: %d x %d (excluding ghost cells)\n", output_cols, output_rows);
    }
    
    if (Con_out) {
        printf("Output calculation result successful\n");
    } else {
        printf("The %d th calculation ended at %f\n", Ite, now_time);
    }
    printf("----------------------------------------\n");
}

/**
 * Perform precision test for convergence analysis
 * Tests the numerical scheme with different grid resolutions
 * Calculates L1, L2, and L∞ errors for each grid
 * Prints convergence table showing errors and orders of accuracy
 */
void Precision_test() {
    
    // Define grid point sequence (consistent with Table 3.3 in the reference paper)
    int nx_values[] = {10, 20, 40, 80, 160, 320, 640};
    int num_nx = sizeof(nx_values) / sizeof(nx_values[0]);
    
    // Initialize error counter
    error_count = 0;
    
    // Loop over different grid resolutions
    for (int idx = 0; idx < num_nx; idx++) {
        int nx = nx_values[idx];
        int LNX_ngc = nx + 2 * GhostCell;      // Total rows including ghost cells
        int LNY_ngc = ctrl_params.L_ny + 2 * GhostCell;  // Total columns including ghost cells

        double Lx, Ly;                         // Domain dimensions
        double Tmax = ctrl_params.Tmax;        // Final simulation time
        double Delta_x, Delta_y;               // Grid spacing in x and y directions
        double Delta_T;                        // Time step

        double mesh_x[nx], mesh_y[ctrl_params.L_ny];     // Mesh coordinates
        double pri_Ver1[4], pri_Ver2[4];                  // Temporary primitive variables
        
        // Allocate memory for primitive variables, conservative variables, and fluxes
        double (*pri)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
        double (*U)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
        double (*FU)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
        double (*GU)[LNX_ngc][LNY_ngc] = malloc(4 * sizeof(double[LNX_ngc][LNY_ngc]));
        
        // Check if memory allocation was successful
        if (pri == NULL || U == NULL || FU == NULL || GU == NULL) {
            fprintf(stderr, "Memory allocation failed in Main\n");
            free(pri); free(U); free(FU); free(GU);
            continue;  // Skip to next grid resolution
        } 
        
        // Set material parameters (gamma, source terms, gravity, CFL number)
        set_material_parameters(ctrl_params.M_gamma, ctrl_params.Source, ctrl_params.Gravity, ctrl_params.CFL);
        
        // Initialize flow field with test case configuration
        Init_Precision(ctrl_params.test_case, var, LNX_ngc, LNY_ngc, GhostCell, pri, U, FU, GU, &Lx, &Ly, &Tmax);

        // Calculate grid spacing
        Delta_x = Lx / nx;
        Delta_y = Ly / ctrl_params.L_ny;

        // Generate computational mesh
        Mesh_2D(nx, ctrl_params.L_ny, mesh_x, mesh_y, Delta_x, Delta_y);

        // Time integration loop
        double Time = 0.0;
        for (Time = 0.0; Time < Tmax; Time = Time + Delta_T) {
            // Calculate time step based on CFL condition
            Delta_T = Precision_Get_Delta_T_2D(var, LNX_ngc, LNY_ngc, U, Delta_x, Delta_y);
            
            // Adjust final time step to exactly reach Tmax
            if (Time + Delta_T >= Tmax) {
                Delta_T = Tmax - Time;
            }
            
            // Perform time advancement using selected Runge-Kutta method
            switch (ctrl_params.Time_ADM) {
                case 1:
                    RK1_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y);
                    break;
                case 2:
                    RK2_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y);
                    break;
                case 3:
                    RK3_TimeAd(ctrl_params.scheme, var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_T, Delta_x, Delta_y);
                    break;
                default:
                    fprintf(stderr, "Error: Invalid time advancement method");
                    break;
            }
        }
        
        // Store current grid index for error data
        int current_index = error_count;
        error_count++;
        
        // Calculate and store all error norms
        TEST_Scheme_Error_L1(var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_x, Time, nx, current_index);
        TEST_Scheme_Error_L2(var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_x, Time, current_index);
        TEST_Scheme_Error_Linf(var, LNX_ngc, LNY_ngc, GhostCell, U, Delta_x, Time, current_index);
        
        printf("\n--- Completed calculation for Nx = %d, Δx = %.6f ---", nx, Delta_x);
        // Free allocated memory
        free(pri);
        free(U);
        free(FU);
        free(GU);
    }
    // Calculate convergence orders from error data
    calculate_orders();
    // Print formatted error table
    print_error_table();
    printf("The program has completed its execution.");
}