#ifndef INITIALIZATION_H
#define INITIALIZATION_H


#include <stdio.h>
#include "Golbal.h"
#include "math.h"
#include "function.h"
#include "Boundary_Condition.h"


// 在头文件中定义测试算例枚举类型
typedef enum {
    TEST_1D_SHOCKTUBE,          // 1D激波管（Sod问题）
    TEST_2D_SHOCKTUBE_CASE1,    // 二维黎曼问题 Case 1
    TEST_2D_SHOCKTUBE_CASE2,    // 二维黎曼问题 Case 2
    TEST_2D_SHOCKTUBE_CASE3,    // 二维黎曼问题 Case 3
    TEST_2D_SHOCKTUBE_CASE4,    // 二维黎曼问题 Case 4
    TEST_2D_SHOCKTUBE_CASE5,    // 二维黎曼问题 Case 5
    TEST_TAYLOR_GREEN_VORTEX,   // 泰勒格林涡
    TEST_GAUSSIAN_PULSE,        // 高斯脉冲
    TEST_KELVIN_HELMHOLTZ,      // KH不稳定性（平滑界面）
    TEST_KELVIN_HELMHOLTZ_SHARP,// KH不稳定性（锐利界面）
    TEST_KELVIN_HELMHOLTZ_VP,   // KH不稳定性（速度扰动）
    TEST_RAYLEIGH_TAYLOR,       // RT不稳定性
    TEST_DOUBLE_MACH_REFLECTION,// 双马赫反射
    TEST_BACKWARD_STEP,         // 后台阶流动
    TEST_BLAST_WAVE,           // 球爆炸问题
    TEST_NOH_PROBLEM           // Noh问题
} TestCase2D;



// 定义圆周率PI常量
#ifndef PI
    #define PI 3.14159265358979323846
#endif

static inline void initEuler2D(int var, int rows, int cols, double (*x)[rows][cols]) {
    int i, j, k;

    #pragma omp parallel for collapse(3)
    for (i = 0; i < var; i++) 
        for (j = 0; j < rows; j++) 
           for (k = 0; k < cols; k++) 
                x[i][j][k] = 0.0;

}

//初始条件
static inline void initEulerpri2D_1DShocktube(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 1D Shock Tube Test Case ===\n");
    
    // Set boundary conditions for 1D shock tube
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_REFLECTION; // Bottom: reflective wall
    bc_config.top = BC_REFLECTION;    // Top: reflective wall

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度
    
    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 ){
                x[0][i][j] = 1.0;
                x[1][i][j] = 0.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 1.0;
            }
            else {
                x[0][i][j] = 0.125;
                x[1][i][j] = 0.0;
                x[2][i][j] = 0.0;
                x[3][i][j] = 0.1;
            }
        }
    }

    *Time = 0.15;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Reflective Wall\n");
    printf("  Top:    Reflective Wall\n");
}

// Case 1
static inline void initEulerpri2D_Shocktube1(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 1 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_OUTFLOW;    // Bottom: outflow
    bc_config.top = BC_OUTFLOW;       // Top: outflow

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 2;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 3;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
        }
    }
    *Time = 0.23;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 2
static inline void initEulerpri2D_Shocktube2(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    
    printf("=== 2D Shock Tube Case 2 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_OUTFLOW;    // Bottom: outflow
    bc_config.top = BC_OUTFLOW;       // Top: outflow

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 2;       // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = -0.3;    // v
                x[3][i][j] = 1;       // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1;       // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = -0.4;    // v
                x[3][i][j] = 1;       // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1.0625;  // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = 0.2145;  // v
                x[3][i][j] = 0.4;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 0.5197;  // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = -1.125;  // v
                x[3][i][j] = 0.4;     // p
            }
        }
    }

    *Time = 0.30;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 3
static inline void initEulerpri2D_Shocktube3(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 3 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_OUTFLOW;    // Bottom: outflow
    bc_config.top = BC_OUTFLOW;       // Top: outflow

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 0.5065;  // ρ
                x[1][i][j] = 0.8939;  // u
                x[2][i][j] = 0;       // v
                x[3][i][j] = 0.35;    // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1.5;     // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = 0;       // v
                x[3][i][j] = 1.5;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1.1;     // ρ
                x[1][i][j] = 0.8939;  // u
                x[2][i][j] = 0.8939;  // v
                x[3][i][j] = 1.1;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 0.5065;  // ρ
                x[1][i][j] = 0;       // u
                x[2][i][j] = 0.8939;  // v
                x[3][i][j] = 0.35;    // p
            }
        }
    }

    *Time = 0.25;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 4
static inline void initEulerpri2D_Shocktube4(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 4 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_OUTFLOW;    // Bottom: outflow
    bc_config.top = BC_OUTFLOW;       // Top: outflow

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 2;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = 0.75;  // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 1;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = 0.5;   // v
                x[3][i][j] = 1;     // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 3;     // ρ
                x[1][i][j] = -0.75; // u
                x[2][i][j] = -0.5;  // v
                x[3][i][j] = 1;     // p
            }
        }
    }

    *Time = 0.30;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// Case 5
static inline void initEulerpri2D_Shocktube5(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {

    printf("=== 2D Shock Tube Case 5 ===\n");

    // Set boundary conditions for 2D Riemann problem
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_OUTFLOW;    // Bottom: outflow
    bc_config.top = BC_OUTFLOW;       // Top: outflow

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度、

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++){
            if( i < rows/2 && j < cols/2) {
                x[0][i][j] = 1;       // ρ
                x[1][i][j] = -0.6259; // u
                x[2][i][j] = 0.1;     // v
                x[3][i][j] = 1;       // p
            }
            else if (i < rows/2 && j >= cols/2) {
                x[0][i][j] = 0.5197;  // ρ
                x[1][i][j] = 0.1;     // u
                x[2][i][j] = 0.1;     // v
                x[3][i][j] = 0.4;     // p
            }
            else if (i >= rows/2 && j < cols/2) {
                x[0][i][j] = 0.8;     // ρ
                x[1][i][j] = 0.1;     // u
                x[2][i][j] = 0.1;     // v
                x[3][i][j] = 1;       // p
            }
            else if (i >= rows/2 && j >= cols/2) {
                x[0][i][j] = 1;       // ρ
                x[1][i][j] = 0.1;     // u
                x[2][i][j] = -0.6259; // v
                x[3][i][j] = 1;       // p
            }
        }
    }

    *Time = 0.25;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}

// 泰勒格林涡问题
// 泰勒格林涡问题
static inline void initEulerpri2D_TaylorGreenVortex(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Taylor-Green Vortex Test Case ===\n");

    // Set boundary conditions for periodic flow
    bc_config.left = BC_PERIODICITY;   // Left: periodic
    bc_config.right = BC_PERIODICITY;  // Right: periodic  
    bc_config.bottom = BC_PERIODICITY; // Bottom: periodic
    bc_config.top = BC_PERIODICITY;    // Top: periodic

    *Lx = 2.0 * M_PI;  // x方向周期
    *Ly = 2.0 * M_PI;  // y方向周期
    double rho0 = 1.0;       // 参考密度
    double p0 = 1.0;         // 参考压力
    double U0 = 1.0;         // 参考速度
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    
    // 预计算常数
    double inv_nx = 1.0 / nx;
    double inv_ny = 1.0 / ny;
    double prefactor = rho0 * U0 * U0 / 4.0;
    
    // 并行初始化整个网格
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // 计算物理坐标位置
            double x_pos = (double)i * inv_nx * (*Lx);
            double y_pos = (double)j * inv_ny * (*Ly);
            
            // 泰勒格林涡速度场
            double sin_x = sin(x_pos);
            double cos_x = cos(x_pos);
            double sin_y = sin(y_pos);
            double cos_y = cos(y_pos);
            
            double u = U0 * sin_x * cos_y;      // u速度
            double v = -U0 * cos_x * sin_y;     // v速度
            double pressure = p0 + prefactor * (cos(2.0 * x_pos) + cos(2.0 * y_pos)); // 压力
            
            // 存储原始变量
            x[0][i][j] = rho0;    // 密度
            x[1][i][j] = u;       // u速度
            x[2][i][j] = v;       // v速度
            x[3][i][j] = pressure;// 压力
        }
    }

    *Time = 5.0;

    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Boundary Conditions:\n");
    printf("  Left:   Periodic\n");
    printf("  Right:  Periodic\n");
    printf("  Bottom: Periodic\n");
    printf("  Top:    Periodic\n");
}

// 二维高斯密度脉冲问题 - 所有边界都为出口边界
static inline void initEulerpri2D_GaussianPulse(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== 2D Gaussian Density Pulse Test Case ===\n");

    // Set boundary conditions for open domain
    bc_config.left = BC_OUTFLOW;      // Left: outflow
    bc_config.right = BC_OUTFLOW;     // Right: outflow  
    bc_config.bottom = BC_OUTFLOW;    // Bottom: outflow
    bc_config.top = BC_OUTFLOW;       // Top: outflow

    int i, j;
    *Lx = 10.0;    // 计算域长度
    *Ly = 10.0;    // 计算域宽度
    double x0 = (*Lx)/2.0;  // 脉冲中心x坐标
    double y0 = (*Ly)/2.0;  // 脉冲中心y坐标
    double A = 1.0;      // 脉冲幅度
    double sigma = 1.0;  // 脉冲宽度
    double rho0 = 1.0;   // 背景密度
    double p0 = 1.0;     // 背景压力
    double u0 = 0.5;     // x方向背景速度
    double v0 = 0.3;     // y方向背景速度
    double x_pos, y_pos, r2;
    double nx = rows - 2* GhostCell;
    double ny = cols - 2* GhostCell;


    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            // 计算物理坐标位置
            x_pos = (double)i / nx * (*Lx);
            y_pos = (double)j / ny * (*Ly);
            
            // 计算到脉冲中心的距离平方
            r2 = (x_pos - x0)*(x_pos - x0) + (y_pos - y0)*(y_pos - y0);
            
            // 设置初始条件：高斯密度脉冲 + 均匀背景流动
            x[0][i][j] = rho0 + A * exp(-r2 / (sigma * sigma));  // 密度：背景 + 高斯脉冲
            x[1][i][j] = u0;      // x方向速度
            x[2][i][j] = v0;      // y方向速度
            x[3][i][j] = p0;      // 压力
        }
    }

    *Time = 4.0;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Gaussian Pulse: Center at (%f, %f), Amplitude = %f\n", x0, y0, A);
    printf("Boundary Conditions:\n");
    printf("  Left:   Outflow\n");
    printf("  Right:  Outflow\n");
    printf("  Bottom: Outflow\n");
    printf("  Top:    Outflow\n");
}


// KH不稳定性测试算例
static inline void initEulerpri2D_KelvinHelmholtz(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Kelvin-Helmholtz Instability Test Case ===\n");

    // Set boundary conditions for KH instability
    bc_config.left = BC_PERIODICITY;   // Left: periodic
    bc_config.right = BC_PERIODICITY;  // Right: periodic  
    bc_config.bottom = BC_REFLECTION;  // Bottom: reflective wall
    bc_config.top = BC_REFLECTION;     // Top: reflective wall

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度
    
    double rho1 = 2.0;    // 下层流体密度
    double rho2 = 1.0;    // 上层流体密度
    double u1 = -0.5;     // 下层流体速度
    double u2 = 0.5;      // 上层流体速度
    double p0 = 2.5;      // 背景压力
    double amplitude = 0.01;  // 扰动幅度
    double width = 0.05;      // 剪切层宽度
    
    double y_pos, y_center, perturbation;
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    

    #pragma omp parallel for collapse(2)
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            // 计算物理坐标位置
            y_pos = (double)j / ny * (*Ly);
            y_center = 0.5 * (*Ly);  // 计算域中心
            
            // 添加正弦扰动
            perturbation = amplitude * sin(4.0 * M_PI * (double)i / nx * (*Lx));
            
            // 计算平滑的密度和速度过渡
            double transition = 0.5 * (1.0 + tanh((y_pos - y_center + perturbation) / width));
            
            // 设置初始条件
            x[0][i][j] = rho1 + (rho2 - rho1) * transition;  // 密度
            x[1][i][j] = u1 + (u2 - u1) * transition;        // x方向速度
            x[2][i][j] = 0.0;                                // y方向速度（初始为零）
            x[3][i][j] = p0;                                 // 压力
        }
    }

    *Time = 1.0;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Fluid Properties:\n");
    printf("  Lower layer: density = %.1f, velocity = %.1f\n", rho1, u1);
    printf("  Upper layer: density = %.1f, velocity = %.1f\n", rho2, u2);
    printf("  Background pressure: %.1f\n", p0);
    printf("  Perturbation amplitude: %.3f\n", amplitude);
    printf("  Shear layer width: %.3f\n", width);
    printf("Boundary Conditions:\n");
    printf("  Left:   Periodic\n");
    printf("  Right:  Periodic\n");
    printf("  Bottom: Reflective Wall\n");
    printf("  Top:    Reflective Wall\n");
}

static inline void initEulerpri2D_KelvinHelmholtz_Sharp(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Kelvin-Helmholtz Instability Test Case (Sharp Interface) ===\n");

    // Set boundary conditions for KH instability
    bc_config.left = BC_PERIODICITY;   // Left: periodic
    bc_config.right = BC_PERIODICITY;  // Right: periodic  
    bc_config.bottom = BC_REFLECTION;  // Bottom: reflective wall
    bc_config.top = BC_REFLECTION;     // Top: reflective wall

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度
    
    double rho1 = 2.0;    // 下层流体密度
    double rho2 = 1.0;    // 上层流体密度
    double u1 = -0.5;     // 下层流体速度
    double u2 = 0.5;      // 上层流体速度
    double p0 = 2.5;      // 背景压力
    double amplitude = 0.01;  // 扰动幅度
    
    double y_pos, y_center, interface_position;
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            // 计算物理坐标位置
            y_pos = (double)j / ny * (*Ly);
            y_center = 0.5 * (*Ly);  // 计算域中心
            
            // 计算界面位置（带扰动）
            interface_position = y_center + amplitude * sin(4.0 * M_PI * (double)i / nx * (*Lx));
            
            // 锐利界面 - 直接判断位置
            if (y_pos < interface_position) {
                // 下层流体
                x[0][i][j] = rho1;  // 密度
                x[1][i][j] = u1;    // x方向速度
            } else {
                // 上层流体
                x[0][i][j] = rho2;  // 密度
                x[1][i][j] = u2;    // x方向速度
            }
            
            x[2][i][j] = 0.0;       // y方向速度（初始为零）
            x[3][i][j] = p0;        // 压力
        }
    }

    *Time = 1.0;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Fluid Properties:\n");
    printf("  Lower layer: density = %.1f, velocity = %.1f\n", rho1, u1);
    printf("  Upper layer: density = %.1f, velocity = %.1f\n", rho2, u2);
    printf("  Density ratio: %.1f\n", rho1/rho2);
    printf("  Velocity difference: %.1f\n", u2 - u1);
    printf("  Background pressure: %.1f\n", p0);
    printf("  Perturbation amplitude: %.3f\n", amplitude);
    printf("  Interface type: Sharp (discontinuous)\n");
    printf("Boundary Conditions:\n");
    printf("  Left:   Periodic\n");
    printf("  Right:  Periodic\n");
    printf("  Bottom: Reflective Wall\n");
    printf("  Top:    Reflective Wall\n");
}

static inline void initEulerpri2D_KelvinHelmholtz_VelocityPerturbation(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Kelvin-Helmholtz Instability Test Case (Velocity Perturbation) ===\n");

    // Set boundary conditions for KH instability
    bc_config.left = BC_PERIODICITY;   // Left: periodic
    bc_config.right = BC_PERIODICITY;  // Right: periodic  
    bc_config.bottom = BC_REFLECTION;  // Bottom: reflective wall
    bc_config.top = BC_REFLECTION;     // Top: reflective wall

    int i, j;
    *Lx = 1.0;    // 计算域长度
    *Ly = 1.0;    // 计算域宽度
    
    double rho1 = 2.0;    // 下层流体密度
    double rho2 = 1.0;    // 上层流体密度
    double u1 = -0.5;     // 下层流体速度
    double u2 = 0.5;      // 上层流体速度
    double p0 = 2.5;      // 背景压力
    double amplitude = 0.05;  // 扰动幅度
    
    double y_pos, y_center, interface_thickness = 0.02;
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            // 计算物理坐标位置
            y_pos = (double)j / ny * (*Ly);
            y_center = 0.5 * (*Ly);  // 计算域中心
            
            // 计算速度扰动 - 在y方向添加正弦扰动
            double vy_perturb = amplitude * sin(4.0 * M_PI * (double)i / nx * (*Lx));
            
            // 设置密度分布（仍然保持两层结构）
            if (y_pos < y_center) {
                // 下层流体
                x[0][i][j] = rho1;  // 密度
                x[1][i][j] = u1;    // x方向速度
            } else {
                // 上层流体
                x[0][i][j] = rho2;  // 密度
                x[1][i][j] = u2;    // x方向速度
            }
            
            // 设置速度扰动 - 主要修改在这里
            x[2][i][j] = vy_perturb;  // y方向速度（添加扰动）
            x[3][i][j] = p0;          // 压力
        }
    }

    *Time = 1.0;
    printf("Computational domain dimensions: Lx = %f, Ly = %f \n", *Lx, *Ly);
    printf("Final simulation time: t = %f \n", *Time);
    printf("Fluid Properties:\n");
    printf("  Lower layer: density = %.1f, velocity = %.1f\n", rho1, u1);
    printf("  Upper layer: density = %.1f, velocity = %.1f\n", rho2, u2);
    printf("  Density ratio: %.1f\n", rho1/rho2);
    printf("  Velocity difference: %.1f\n", u2 - u1);
    printf("  Background pressure: %.1f\n", p0);
    printf("  Velocity perturbation amplitude: %.3f\n", amplitude);
    printf("  Perturbation type: Vertical velocity perturbation\n");
    printf("  Perturbation wavenumber: 4\n");
    printf("Boundary Conditions:\n");
    printf("  Left:   Periodic\n");
    printf("  Right:  Periodic\n");
    printf("  Bottom: Reflective Wall\n");
    printf("  Top:    Reflective Wall\n");
}


double L_fixed_value_state[4] = {0};
double R_fixed_value_state[4] = {0};
double B_fixed_value_state[4] = {2.0,0,0,1.5};
double T_fixed_value_state[4] = {1.0,0,0,15.0/4};
static inline void initEulerpri2D_RayleighTaylor(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Rayleigh-Taylor Instability Test Case ===\n");

     // 设置边界条件
    bc_config.left = BC_REFLECTION;      // 左右边界：反射
    bc_config.right = BC_REFLECTION;
    bc_config.bottom = BC_FIXED_VALUE;  // 上下边界：固定值
    bc_config.top = BC_FIXED_VALUE;

    // 设置计算域
    *Lx = 0.25;
    *Ly = 1.0;
    *Time = 1.95;

    // 初始化流场
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // 计算物理坐标
            double x_pos = (double)(i + 0.5 - GhostCell) / nx * (*Lx);
            double y_pos = (double)(j + 0.5 - GhostCell) / ny * (*Ly);
            
            // 根据y坐标位置设置不同的初始状态
            if (y_pos < 0.5) {
                // 下层流体: 0 < y < 1/2
                x[0][i][j] = 2.0;                           // 密度 ρ = 2
                x[1][i][j] = 0.0;                           // x动量 ρu = 0
                x[2][i][j] = -0.025 *sqrt((2.0 * y_pos + 1.0)*M_gamma/2.0) * cos(8.0 * M_PI * x_pos); // y动量 ρv = ρ * (-0.025 cos(8πx))
                //x[2][i][j] = 0.0; // y动量 ρv = ρ * (-0.025 cos(8πx))
                x[3][i][j] = 2.0 * y_pos + 1.0;             // 压力 p = 2y + 1
            } else {
                // 上层流体: 1/2 < y ≤ 1
                x[0][i][j] = 1.0;                           // 密度 ρ = 1
                x[1][i][j] = 0.0;                           // x动量 ρu = 0
                x[2][i][j] = -0.025 *sqrt((y_pos + 1.5)*M_gamma/1.0) * cos(8.0 * M_PI * x_pos); // y动量 ρv = ρ * (-0.025 cos(8πx))
                //x[2][i][j] = 0.0; // y动量 ρv = ρ * (-0.025 cos(8πx))
                x[3][i][j] = y_pos + 1.5;                   // 压力 p = y + 3/2
            }
        }
    }

    printf("Computational domain: [0, %f] x [0, %f]\n", *Lx, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Density ratio: 2:1 (heavy fluid below, light fluid above)\n");
    printf("Initial perturbation: v = -0.025 cos(8πx)\n");
    printf("Gravity direction: -y (implied by pressure gradient)\n");

}

// 二维双马赫反射问题
// 全局入口状态定义（在某个源文件中）
double inflow_state[4] = {1.4, 0.0, 0.0, 1.0}; // 默认值
// 双马赫反射初始化函数
double post_shock_state[4] = {8.0, 7.145, -4.125, 116.83333};
double pre_shock_state[4] = {1.4, 0.0, 0.0, 1.0};
double shock_slope = 1.732;  // tan(60°)
double shock_start_x = 1.0/6.0;
double current_time = 0.0;   // 当前时间，需要在主循环中更新
// 双马赫反射问题初始化

static inline void initEulerpri2D_DoubleMachReflection(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Double Mach Reflection Test Case ===\n");

    // 设置边界条件
    bc_config.left = BC_INFLOW;                 // 左边界：入流（激波后状态）
    bc_config.right = BC_OUTFLOW;               // 右边界：出流
//    bc_config.bottom = BC_DOUBLE_MACH_BOTTOM;   // 下边界：特殊处理
//    bc_config.top = BC_DOUBLE_MACH_TOP;         // 上边界：特殊处理

    // 设置计算域
    *Lx = 4.0;
    *Ly = 1.0;
    *Time = 0.2;

    // 设置激波状态（根据论文公式22）
    post_shock_state[0] = 8.0;          // 密度
    post_shock_state[1] = 7.145;        // x方向速度
    post_shock_state[2] = -4.125;       // y方向速度
    post_shock_state[3] = 116.83333;    // 压力

    pre_shock_state[0] = 1.4;           // 密度
    pre_shock_state[1] = 0.0;           // x方向速度
    pre_shock_state[2] = 0.0;           // y方向速度
    pre_shock_state[3] = 1.0;           // 压力

    shock_slope = 1.732;                // tan(60°)
    shock_start_x = 1.0/6.0;

    // 初始化流场
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows-GhostCell; i++) {
        for (int j = GhostCell; j < cols-GhostCell; j++) {
            // 计算物理坐标
            double x_pos = (double)(i + 0.5 - GhostCell) / nx * (*Lx);
            double y_pos = (double)(j + 0.5 - GhostCell) / ny * (*Ly);
            
            // 根据激波位置初始化流场
            // 激波线方程: y = shock_slope * (x - shock_start_x)
            if (y_pos <= shock_slope * (x_pos - shock_start_x)) {
                // 激波后区域
                x[0][i][j] = pre_shock_state[0];   // 密度
                x[1][i][j] = pre_shock_state[1];   // x速度
                x[2][i][j] = pre_shock_state[2];   // y速度
                x[3][i][j] = pre_shock_state[3];   // 压力
                
            } else {
                x[0][i][j] = post_shock_state[0];  // 密度
                x[1][i][j] = post_shock_state[1];  // x速度
                x[2][i][j] = post_shock_state[2];  // y速度
                x[3][i][j] = post_shock_state[3];  // 压力
               
            }
        }
    }

    printf("Computational domain: [0, %f] x [0, %f]\n", *Lx, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Shock Mach number: 10\n");
    printf("Shock angle: 60 degrees\n");
    printf("Pre-shock state: (ρ, u, v, p) = (%f, %f, %f, %f)\n", 
           pre_shock_state[0], pre_shock_state[1], pre_shock_state[2], pre_shock_state[3]);
    printf("Post-shock state: (ρ, u, v, p) = (%f, %f, %f, %f)\n", 
           post_shock_state[0], post_shock_state[1], post_shock_state[2], post_shock_state[3]);
    printf("Reflecting wall from x = %f to x = %f\n", 1.0/6.0, 4.0);
}




static inline void initEulerpri2D_BackwardStep(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Woodward & Colella Backward Step Flow Test Case (1984) ===\n");

    // 设置边界条件（根据论文描述）
    bc_config.left = BC_INFLOW;           // 左边界：超声速入流
    bc_config.right = BC_OUTFLOW;         // 右边界：出流
    bc_config.bottom = BC_OUTFLOW;     // 下边界：反射壁面
    bc_config.top = BC_OUTFLOW;        // 上边界：反射壁面

    // 计算域尺寸（根据论文图12）
    *Lx = 3.0;     // 长度：从x=0到x=3
    *Ly = 1.0;     // 高度：从y=0到y=1
    *Time = 4.0;   // 计算时间：t=4.0（论文中使用这个时间）

    // 台阶几何参数（根据论文描述）
    double step_height = 0.2;           // 台阶高度：通道高度的1/5
    double step_position = 0.6;         // 台阶位置：x=0.6（距离入口）
    double channel_height = 1.0;        // 通道总高度
    double entrance_length = 0.6;       // 入口段长度（台阶前）

    // 来流条件（根据论文：马赫数3的超声速流动）
    double gamma = 1.4;                 // 比热比
    double inflow_mach = 3.0;           // 马赫数3
    double inflow_density = 1.4;        // 密度（标准化）
    double inflow_pressure = 1.0;       // 压力（标准化）
    
    // 计算声速和速度
    double sound_speed = sqrt(gamma * inflow_pressure / inflow_density);
    double inflow_velocity_x = inflow_mach * sound_speed;  // x方向速度
    double inflow_velocity_y = 0.0;                        // y方向速度
    
    // 计算总能量
    double inflow_energy = inflow_pressure / (gamma - 1.0) + 
                          0.5 * inflow_density * (inflow_velocity_x * inflow_velocity_x);

    printf("Woodward & Colella (1984) setup:\n");
    printf("Computational domain: [0, %f] x [0, %f]\n", *Lx, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Step configuration:\n");
    printf("  - Step height: %f (20%% of channel height)\n", step_height);
    printf("  - Step position: x = %f\n", step_position);
    printf("  - Entrance length before step: %f\n", entrance_length);
    printf("Inflow conditions (Mach 3):\n");
    printf("  - Mach number: %f\n", inflow_mach);
    printf("  - Density: %f\n", inflow_density);
    printf("  - Pressure: %f\n", inflow_pressure);
    printf("  - Velocity: (%f, %f)\n", inflow_velocity_x, inflow_velocity_y);
    printf("  - Sound speed: %f\n", sound_speed);
    printf("  - Gamma: %f\n", gamma);

    // 计算网格参数
    int nx = rows - 2 * GhostCell;      // x方向有效网格数
    int ny = cols - 2 * GhostCell;      // y方向有效网格数
    double dx = *Lx / nx;
    double dy = *Ly / ny;

    // 初始化整个流场为来流状态
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // 计算物理坐标（网格中心坐标）
            double x_pos = (i + 0.5 - GhostCell) * dx;
            double y_pos = (j + 0.5 - GhostCell) * dy;
            
            // 判断是否在台阶固体区域内
            int is_in_solid = 0;
            
            // 台阶固体区域定义（根据论文图12）：
            // 1. x < 0.6 且 y < 0.2 的区域（台阶本身）
            // 2. x = 0.6 且 y < 0.2 的垂直面（台阶前缘）
            if (x_pos < step_position && y_pos < step_height) {
                is_in_solid = 1;  // 台阶主体区域
            } else if (fabs(x_pos - step_position) < 0.5 * dx && y_pos < step_height) {
                is_in_solid = 1;  // 台阶垂直面
            }
            
            if (is_in_solid) {
                // 固体区域：设置为高密度静止流体（壁面）
                // 实际上这些单元格会被边界条件覆盖，但初始化为壁面状态
                x[0][i][j] = 10.0 * inflow_density;  // 高密度，表示固体
                x[1][i][j] = 0.0;                    // 零速度
                x[2][i][j] = 0.0;
                x[3][i][j] = inflow_pressure;        // 压力与来流相同
            } else {
                // 流体区域：设置为均匀来流条件
                x[0][i][j] = inflow_density;
                x[1][i][j] = inflow_density * inflow_velocity_x;  // 动量
                x[2][i][j] = 0.0;                                // y方向动量
                x[3][i][j] = inflow_energy;                      // 总能量
                
                // 在台阶拐角处设置小的扰动（可选，帮助流动发展）
                double dx_to_corner = x_pos - step_position;
                double dy_to_corner = y_pos - step_height;
                double dist_to_corner = sqrt(dx_to_corner*dx_to_corner + dy_to_corner*dy_to_corner);
                
                if (dist_to_corner < 0.1) {
                    // 在拐角附近添加微小扰动（模拟真实流动的初始分离）
                    double perturbation = 0.01 * exp(-dist_to_corner / 0.02);
                    x[2][i][j] = inflow_density * perturbation * inflow_velocity_x;
                }
            }
        }
    }

    // 验证初始条件的物理合理性
    printf("\nInitial condition verification:\n");
    printf("Total cells in solid region (step): approximately %d\n", 
           (int)(step_position / dx * step_height / dy));
    printf("Total cells in fluid region: approximately %d\n", 
           nx * ny - (int)(step_position / dx * step_height / dy));
    printf("Flow expansion ratio: %f\n", 
           channel_height / (channel_height - step_height));
    
    // 计算并显示一些关键参数
    double reynolds_number = 0.0;  // 无粘欧拉方程，雷诺数无穷大
    printf("Reynolds number: infinite (inviscid Euler equations)\n");
    printf("Flow features expected (from Woodward & Colella 1984):\n");
    printf("  1. Strong shock from step corner\n");
    printf("  2. Expansion fan at step corner\n");
    printf("  3. Recompression shock downstream\n");
    printf("  4. Separation region behind step\n");
    printf("  5. Shear layer development\n");
}


static inline void initEulerpri2D_BlastWave(int var, int rows, int cols, 
                                           double (*x)[rows][cols], 
                                           double *Lx, double *Ly, double *Time) {
    printf("=== 2D Blast Wave Test Case ===\n");
    printf("Strong spherical explosion with initial pressure discontinuity\n");
    
    // 设置边界条件
    bc_config.left = BC_REFLECTION;    // 反射边界
    bc_config.right = BC_REFLECTION;   // 反射边界
    bc_config.bottom = BC_REFLECTION;  // 反射边界
    bc_config.top = BC_REFLECTION;     // 反射边界
    
    // 计算域尺寸
    *Lx = 1.0;    // [-0.5, 0.5]
    *Ly = 1.0;    // [-0.5, 0.5]
    *Time = 0.038; // 文献中的典型时间
    
    // 物理参数
    double gamma = 1.4;  // 比热比
    
    // 高压区（中心爆炸）
    double high_pressure = 10.0;   // 高压
    double high_density = 1.0;     // 密度
    
    // 低压区（外部环境）
    double low_pressure = 0.1;     // 低压
    double low_density = 1.0;      // 密度相同，便于观察压力效应
    
    // 初始速度为零
    double init_velocity_x = 0.0;
    double init_velocity_y = 0.0;
    
    // 爆炸半径
    double blast_radius = 0.1;
    
    // 初始化流场
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    double dx = *Lx / nx;
    double dy = *Ly / ny;
    
    // 计算域中心坐标（用于球对称）
    double center_x = 0.0;
    double center_y = 0.0;
    
    int high_pressure_cells = 0;
    int low_pressure_cells = 0;
    

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // 计算物理坐标（从-0.5到0.5）
            double x_pos = (i + 0.5 - GhostCell) * dx - 0.5;
            double y_pos = (j + 0.5 - GhostCell) * dy - 0.5;
            
            // 计算到中心的距离
            double dx_center = x_pos - center_x;
            double dy_center = y_pos - center_y;
            double radius = sqrt(dx_center * dx_center + dy_center * dy_center);
            
            if (radius <= blast_radius) {
                // 高压爆炸区域
                x[0][i][j] = high_density;
                x[1][i][j] = high_density * init_velocity_x;
                x[2][i][j] = high_density * init_velocity_y;
                x[3][i][j] = high_pressure/(gamma-1.0) + 
                            0.5 * high_density * (init_velocity_x * init_velocity_x + 
                                                init_velocity_y * init_velocity_y);
                high_pressure_cells++;
            } else {
                // 低压环境区域
                x[0][i][j] = low_density;
                x[1][i][j] = low_density * init_velocity_x;
                x[2][i][j] = low_density * init_velocity_y;
                x[3][i][j] = low_pressure/(gamma-1.0) + 
                            0.5 * low_density * (init_velocity_x * init_velocity_x + 
                                               init_velocity_y * init_velocity_y);
                low_pressure_cells++;
            }
        }
    }
    
    printf("Computational domain: [%f, %f] x [%f, %f]\n", -0.5, 0.5, -0.5, 0.5);
    printf("Final time: t = %f\n", *Time);
    printf("High pressure region (r <= %f):\n", blast_radius);
    printf("  Density: ρ = %f, Pressure: p = %f\n", high_density, high_pressure);
    printf("Low pressure region (r > %f):\n", blast_radius);
    printf("  Density: ρ = %f, Pressure: p = %f\n", low_density, low_pressure);
    printf("Initial velocity: (u, v) = (%f, %f)\n", init_velocity_x, init_velocity_y);
    printf("Gamma: γ = %f\n", gamma);
    printf("Grid cells in high pressure region: %d\n", high_pressure_cells);
    printf("Grid cells in low pressure region: %d\n", low_pressure_cells);
    printf("Pressure ratio: %f\n", high_pressure / low_pressure);
    printf("Expected features:\n");
    printf("  1. Strong outward propagating shock wave\n");
    printf("  2. Contact discontinuity separating fluids\n");
    printf("  3. Rarefaction wave moving inward\n");
    printf("  4. Complex shock interactions at late times\n");
}


static inline void initEulerpri2D_NohProblem(int var, int rows, int cols,
                                            double (*x)[rows][cols],
                                            double *Lx, double *Ly, double *Time) {
    printf("=== 2D Noh Problem Test Case ===\n");
    printf("Strong converging shock wave with infinite strength\n");
    
    // 设置边界条件（球对称问题的近似）
    bc_config.left = BC_REFLECTION;    // 反射边界
    bc_config.right = BC_REFLECTION;   // 反射边界
    bc_config.bottom = BC_REFLECTION;  // 反射边界
    bc_config.top = BC_REFLECTION;     // 反射边界
    
    // 计算域尺寸（第一象限，利用对称性）
    *Lx = 1.0;    // [0, 1]
    *Ly = 1.0;    // [0, 1]
    *Time = 0.6;  // 典型计算时间
    
    // 物理参数
    double gamma = 5.0/3.0;  // Noh问题常用γ=5/3
    
    // 初始条件（均匀介质向内坍塌）
    double init_density = 1.0;
    double init_pressure = 1.0e-6;  // 非常小的压力，近似冷气体
    
    // 向内径向速度（大小为1，指向原点）
    double init_speed = 1.0;  // 向内速度大小
    
    // 初始化流场
    double nx = rows - 2 * GhostCell;
    double ny = cols - 2 * GhostCell;
    double dx = *Lx / nx;
    double dy = *Ly / ny;
    
    // 原点坐标（左下角）
    double origin_x = 0.0;
    double origin_y = 0.0;
    
    // 添加微小扰动避免完全对称（可选，有助于数值稳定性）
    double perturbation_amplitude = 0.001;
    

    #pragma omp parallel for collapse(2)
    for (int i = GhostCell; i < rows - GhostCell; i++) {
        for (int j = GhostCell; j < cols - GhostCell; j++) {
            // 计算物理坐标
            double x_pos = (i + 0.5 - GhostCell) * dx;
            double y_pos = (j + 0.5 - GhostCell) * dy;
            
            // 计算到原点的距离
            double distance = sqrt(x_pos * x_pos + y_pos * y_pos);
            
            // 计算径向单位向量（指向原点）
            double u_comp, v_comp;
            if (distance > 1.0e-10) {
                u_comp = -x_pos / distance;  // 指向原点的x分量
                v_comp = -y_pos / distance;  // 指向原点的y分量
            } else {
                u_comp = 0.0;
                v_comp = 0.0;
            }
            
            // 添加微小随机扰动避免完全对称
            double perturbation = perturbation_amplitude * 
                                 (2.0 * ((double)rand()/RAND_MAX) - 1.0);
            u_comp += perturbation;
            v_comp += perturbation;
            
            // 归一化保持速度大小
            double norm = sqrt(u_comp * u_comp + v_comp * v_comp);
            if (norm > 1.0e-10) {
                u_comp = u_comp * init_speed / norm;
                v_comp = v_comp * init_speed / norm;
            }
            
            // 设置守恒变量
            x[0][i][j] = init_density;
            x[1][i][j] = init_density * u_comp;
            x[2][i][j] = init_density * v_comp;
            x[3][i][j] = init_pressure/(gamma-1.0) + 
                        0.5 * init_density * (u_comp * u_comp + v_comp * v_comp);
            
            // 在中心附近区域添加微小密度扰动（可选）
            if (distance < 0.05) {
                x[0][i][j] = init_density * (1.0 + perturbation_amplitude * 
                                            ((double)rand()/RAND_MAX - 0.5));
            }
        }
    }
    
    // Noh问题的理论解参数
    double shock_speed = (gamma - 1.0) / 2.0;  // 激波速度
    double density_jump = (gamma + 1.0) / (gamma - 1.0);  // 密度跳跃比
    double post_shock_density = init_density * density_jump;
    double post_shock_pressure = 0.5 * init_density * init_speed * init_speed * 
                                (gamma + 1.0);
    
    printf("Computational domain: [%f, %f] x [%f, %f]\n", 0.0, *Lx, 0.0, *Ly);
    printf("Final time: t = %f\n", *Time);
    printf("Initial conditions:\n");
    printf("  Density: ρ₀ = %f\n", init_density);
    printf("  Pressure: p₀ = %e (approximately zero)\n", init_pressure);
    printf("  Radial velocity: v_r = -%f (inward)\n", init_speed);
    printf("  Gamma: γ = %f (5/3)\n", gamma);
    printf("\nTheoretical post-shock state (self-similar solution):\n");
    printf("  Shock speed: v_s = %f\n", shock_speed);
    printf("  Shock position at t=%f: r_s = %f\n", *Time, shock_speed * (*Time));
    printf("  Density jump ratio: %f\n", density_jump);
    printf("  Post-shock density: ρ₁ = %f\n", post_shock_density);
    printf("  Post-shock pressure: p₁ = %f\n", post_shock_pressure);
    printf("  Post-shock velocity: v₁ = 0 (stagnation)\n");
    printf("\nNumerical challenges:\n");
    printf("  1. Infinite strength shock (p₀ ≈ 0)\n");
    printf("  2. Strong convergence to center singularity\n");
    printf("  3. Preservation of spherical symmetry\n");
    printf("  4. Shock overheating problem verification\n");
    
    // 计算激波最终位置
    double shock_final_radius = shock_speed * (*Time);
    printf("\nExpected shock front at t=%f: radius = %f\n", *Time, shock_final_radius);
    
    if (shock_final_radius > 0.5) {
        printf("Warning: Shock may reach boundary before final time!\n");
        printf("Consider reducing final time or increasing domain size.\n");
    }
}


static inline void initEulerconser2D(int var, int rows, int cols, double (*x)[rows][cols], double (*y)[rows][cols]) {
    int j,k;

    #pragma omp parallel for collapse(2)
    for (j = 0; j < rows; j++) {
        for ( k = 0; k < cols; k++){
            double rho, u, v, p;
            rho = x[0][j][k];
            u = x[1][j][k];
            v = x[2][j][k];
            p = x[3][j][k];
            y[0][j][k] = rho;
            y[1][j][k] = u * rho;
            y[2][j][k] = v * rho;
            y[3][j][k] = p /(M_gamma-1) + 0.5 * rho * (pow(u,2) + pow(v,2));
        }
    }
}



static inline void initEulerflux2D(int var, int rows, int cols, double (*y)[rows][cols], double (*f)[rows][cols],double (*g)[rows][cols]) {
    int j,k;

    #pragma omp parallel for collapse(2)

    for (j = 0; j < rows; j++) {
        for ( k = 0; k < cols; k++){
            double rho, u, v, p;
            rho = y[0][j][k];
            u = y[1][j][k]/y[0][j][k];
            v = y[2][j][k]/y[0][j][k];
            p = (y[3][j][k] - 0.5 * rho * (pow(u,2) + pow(v,2)))*(M_gamma-1);
            //x dir initilize
            f[0][j][k] = rho * u;
            f[1][j][k] = rho * u * u + p;
            f[2][j][k] = rho * u * v;
            f[3][j][k] = (y[3][j][k] + p) * u;
            //y dir initilize
            g[0][j][k] = rho * v;
            g[1][j][k] = rho * u * v;
            g[2][j][k] = rho * v * v + p;
            g[3][j][k] = (y[3][j][k] + p) * v;

        }
    }
}

static inline void Con_to_Pri_2D(int var, int rows, int cols, double (*x)[rows][cols], double (*y)[rows][cols]){
   int j,k;

    #pragma omp parallel for collapse(2)
    for (j = 0; j < rows; j++) {
        for ( k = 0; k < cols; k++){
            double rho, u, v, p;
            rho = y[0][j][k];
            u = y[1][j][k]/y[0][j][k];
            v = y[2][j][k]/y[0][j][k];
            p = (y[3][j][k] - 0.5 * rho * (pow(u,2) + pow(v,2)))*(M_gamma-1);
            //x dir initilize
            x[0][j][k] = rho;
            x[1][j][k] = u;
            x[2][j][k] = v;
            x[3][j][k] = p;
        }
    }
}


// 统一的欧拉方程2D测试算例初始化函数
static inline void initEulerTestCase(TestCase2D test_case, 
                                    int var, int rows, int cols,
                                    double (*x)[rows][cols],
                                    double *Lx, double *Ly, double *Time) {
    
    // 根据测试算例选择相应的初始化函数
    switch(test_case) {
        case TEST_1D_SHOCKTUBE:
            initEulerpri2D_1DShocktube(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE1:
            initEulerpri2D_Shocktube1(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE2:
            initEulerpri2D_Shocktube2(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE3:
            initEulerpri2D_Shocktube3(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE4:
            initEulerpri2D_Shocktube4(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_2D_SHOCKTUBE_CASE5:
            initEulerpri2D_Shocktube5(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_TAYLOR_GREEN_VORTEX:
            initEulerpri2D_TaylorGreenVortex(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_GAUSSIAN_PULSE:
            initEulerpri2D_GaussianPulse(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_KELVIN_HELMHOLTZ:
            initEulerpri2D_KelvinHelmholtz(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_KELVIN_HELMHOLTZ_SHARP:
            initEulerpri2D_KelvinHelmholtz_Sharp(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_KELVIN_HELMHOLTZ_VP:
            initEulerpri2D_KelvinHelmholtz_VelocityPerturbation(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_RAYLEIGH_TAYLOR:
            initEulerpri2D_RayleighTaylor(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_DOUBLE_MACH_REFLECTION:
            initEulerpri2D_DoubleMachReflection(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_BACKWARD_STEP:
            initEulerpri2D_BackwardStep(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_BLAST_WAVE:
            initEulerpri2D_BlastWave(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        case TEST_NOH_PROBLEM:
            initEulerpri2D_NohProblem(var, rows, cols, x, Lx, Ly, Time);
            break;
            
        default:
            printf("Error: Unknown test case selected!\n");
            // 默认使用1D激波管
            initEulerpri2D_1DShocktube(var, rows, cols, x, Lx, Ly, Time);
            break;
    }
}


static inline void Init_Euler_2D(TestCase2D test_case,
                                int var, int rows, int cols, int GC, 
                                double (*x)[rows][cols], double (*y)[rows][cols], 
                                double (*f)[rows][cols], double (*g)[rows][cols], 
                                double *Lx, double *Ly, double *Time) {
    
    printf("========== Initializing 2D Euler Equations ==========\n");
    
    // 初始化数组（清零）
    initEuler2D(var, rows, cols, x);
    initEuler2D(var, rows, cols, y);
    initEuler2D(var, rows, cols, f);
    initEuler2D(var, rows, cols, g);
    
    // 使用统一的测试算例初始化函数
    printf("Selected test case: ");
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: printf("1D Sod Shock Tube\n"); break;
        case TEST_2D_SHOCKTUBE_CASE1: printf("2D Riemann Problem Case 1\n"); break;
        case TEST_2D_SHOCKTUBE_CASE2: printf("2D Riemann Problem Case 2\n"); break;
        case TEST_2D_SHOCKTUBE_CASE3: printf("2D Riemann Problem Case 3\n"); break;
        case TEST_2D_SHOCKTUBE_CASE4: printf("2D Riemann Problem Case 4\n"); break;
        case TEST_2D_SHOCKTUBE_CASE5: printf("2D Riemann Problem Case 5\n"); break;
        case TEST_TAYLOR_GREEN_VORTEX: printf("Taylor-Green Vortex\n"); break;
        case TEST_GAUSSIAN_PULSE: printf("Gaussian Density Pulse\n"); break;
        case TEST_KELVIN_HELMHOLTZ: printf("Kelvin-Helmholtz Instability (Smooth)\n"); break;
        case TEST_KELVIN_HELMHOLTZ_SHARP: printf("Kelvin-Helmholtz Instability (Sharp)\n"); break;
        case TEST_KELVIN_HELMHOLTZ_VP: printf("Kelvin-Helmholtz Instability (Velocity Perturbation)\n"); break;
        case TEST_RAYLEIGH_TAYLOR: printf("Rayleigh-Taylor Instability\n"); break;
        case TEST_DOUBLE_MACH_REFLECTION: printf("Double Mach Reflection\n"); break;
        case TEST_BACKWARD_STEP: printf("Backward Step Flow\n"); break;
        case TEST_BLAST_WAVE: printf("Blast Wave (Spherical Explosion)\n"); break;
        case TEST_NOH_PROBLEM: printf("Noh Problem (Converging Shock)\n"); break;
    }
    
    // 调用统一的初始化函数
    initEulerTestCase(test_case, var, rows, cols, x, Lx, Ly, Time);
    
    // 将原始变量转换为守恒变量
    initEulerconser2D(var, rows, cols, x, y);
    
    // 计算通量
    initEulerflux2D(var, rows, cols, y, f, g);
    
    printf("========== Initialization Complete ==========\n\n");
}

// 获取测试算例的描述信息
static inline const char* getTestCaseDescription(TestCase2D test_case) {
    switch(test_case) {
        case TEST_1D_SHOCKTUBE: return "1D Sod shock tube problem - classical benchmark";
        case TEST_2D_SHOCKTUBE_CASE1: return "2D Riemann problem case 1 - four interacting states";
        case TEST_2D_SHOCKTUBE_CASE2: return "2D Riemann problem case 2 - complex wave interactions";
        case TEST_2D_SHOCKTUBE_CASE3: return "2D Riemann problem case 3 - shock interactions";
        case TEST_2D_SHOCKTUBE_CASE4: return "2D Riemann problem case 4 - symmetric 4-quadrant";
        case TEST_2D_SHOCKTUBE_CASE5: return "2D Riemann problem case 5 - shock-rarefaction";
        case TEST_TAYLOR_GREEN_VORTEX: return "Taylor-Green vortex - periodic decaying turbulence";
        case TEST_GAUSSIAN_PULSE: return "Gaussian density pulse in uniform flow";
        case TEST_KELVIN_HELMHOLTZ: return "Kelvin-Helmholtz instability - shear layer with smooth interface";
        case TEST_KELVIN_HELMHOLTZ_SHARP: return "Kelvin-Helmholtz instability - sharp interface";
        case TEST_KELVIN_HELMHOLTZ_VP: return "Kelvin-Helmholtz instability - velocity perturbation";
        case TEST_RAYLEIGH_TAYLOR: return "Rayleigh-Taylor instability - heavy fluid over light fluid";
        case TEST_DOUBLE_MACH_REFLECTION: return "Double Mach reflection - strong shock interaction";
        case TEST_BACKWARD_STEP: return "Backward step flow - supersonic flow with separation";
        case TEST_BLAST_WAVE: return "Blast wave - strong spherical explosion";
        case TEST_NOH_PROBLEM: return "Noh problem - infinite strength converging shock";
        default: return "Unknown test case";
    }
}

// 打印所有可用测试算例
void printAvailableTestCases() {
    printf("\n╔══════════════════════════════════════════════════════════════════════╗\n");
    printf("║                    Available 2D Euler Test Cases             ║\n");
    printf("╠══════════════════════════════════════════════════════════════╣\n");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_1D_SHOCKTUBE, "1D Sod", "Classic 1D shock tube");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_2D_SHOCKTUBE_CASE1, "2D Riemann Case1", "Four interacting states");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_2D_SHOCKTUBE_CASE2, "2D Riemann Case2", "Complex wave interactions");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_2D_SHOCKTUBE_CASE3, "2D Riemann Case3", "Shock interactions");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_2D_SHOCKTUBE_CASE4, "2D Riemann Case4", "Symmetric 4-quadrant");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_2D_SHOCKTUBE_CASE5, "2D Riemann Case5", "Shock-rarefaction");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_TAYLOR_GREEN_VORTEX, "Taylor-Green Vortex", "Periodic decaying turbulence");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_GAUSSIAN_PULSE, "Gaussian Pulse", "Density pulse in uniform flow");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_KELVIN_HELMHOLTZ, "KH Smooth", "Shear layer with smooth interface");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_KELVIN_HELMHOLTZ_SHARP, "KH Sharp", "Shear layer with sharp interface");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_KELVIN_HELMHOLTZ_VP, "KH VP", "Velocity perturbation");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_RAYLEIGH_TAYLOR, "Rayleigh-Taylor", "Heavy over light fluid");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_DOUBLE_MACH_REFLECTION, "Double Mach", "Strong shock reflection");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_BACKWARD_STEP, "Backward Step", "Supersonic flow with separation");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_BLAST_WAVE, "Blast Wave", "Spherical explosion");
    printf("║ %2d: %-25s - %-35s ║\n", TEST_NOH_PROBLEM, "Noh Problem", "Converging shock wave");
    printf("╚═══════════════════════════════════════════════════════════════════╝\n\n");
}


#endif