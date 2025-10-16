#ifndef INITIALIZATION_H
#define INITIALIZATION_H


#include <stdio.h>
#include "Golbal.h"
#include "math.h"
#include "function.h"
#include "Boundary_Condition.h"


// 定义圆周率PI常量
#ifndef PI
    #define PI 3.14159265358979323846
#endif

static inline void initEuler2D(int var, int rows, int cols, double (*x)[rows][cols]) {
    int i, j, k;
    for (i = 0; i < var; i++) 
        for (j = 0; j < rows; j++) 
           for (k = 0; k < cols; k++) 
                x[i][j][k] = 0.0;

}

//初始条件
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
    *Ly = 1.0;    // 计算域宽度

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
static inline void initEulerpri2D_TaylorGreenVortex(int var, int rows, int cols, double (*x)[rows][cols], double *Lx, double *Ly, double *Time) {
    printf("=== Taylor-Green Vortex Test Case ===\n");

    // Set boundary conditions for periodic flow
    bc_config.left = BC_PERIODICITY;   // Left: periodic
    bc_config.right = BC_PERIODICITY;  // Right: periodic  
    bc_config.bottom = BC_PERIODICITY; // Bottom: periodic
    bc_config.top = BC_PERIODICITY;    // Top: periodic

    int i, j;
    *Lx = 2.0 * M_PI;  // x方向周期
    *Ly = 2.0 * M_PI;  // y方向周期
    double rho0 = 1.0;       // 参考密度
    double p0 = 1.0;         // 参考压力
    double U0 = 1.0;         // 参考速度
    double x_pos, y_pos;
    double nx = rows - 2* GhostCell;
    double ny = cols - 2* GhostCell;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            // 计算物理坐标位置
            x_pos = (double)i / nx * (*Lx);
            y_pos = (double)j / ny * (*Ly);
            
            
            // 泰勒格林涡速度场
            x[0][i][j] = rho0;  // 密度
            x[1][i][j] = U0 * sin(x_pos) * cos(y_pos);      // u速度
            x[2][i][j] = -U0 * cos(x_pos) * sin(y_pos);     // v速度
            x[3][i][j] = p0 + (rho0 * U0 * U0 / 4.0) * (cos(2.0 * x_pos) + cos(2.0 * y_pos)); // 压力
        }
    }

    *Time = 1.0;

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



static inline void initEulerconser2D(int var, int rows, int cols, double (*x)[rows][cols], double (*y)[rows][cols]) {
    int j,k;
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



static inline void Init_Euler_2D(int var,int rows, int cols, int GC, double (*x)[rows][cols], double (*y)[rows][cols], \
                                    double (*f)[rows][cols], double (*g)[rows][cols], double *Lx, double *Ly, double *Time){
    //初始化
    initEuler2D(var, rows, cols,x);
    initEuler2D(var, rows, cols,y);
    initEuler2D(var, rows, cols,f);
    initEuler2D(var, rows, cols,g);
   
//    initEulerpri2D_1DShocktube(var, rows,cols, x, Time);
//    initEulerpri2D_Shocktube5(var, rows, cols, x, Time);
    initEulerpri2D_KelvinHelmholtz_Sharp(var, rows, cols, x, Lx, Ly, Time);
//    initEulerpri2D_GaussianPulse(var, rows, cols, x, Lx, Ly, Time);
    initEulerconser2D(var, rows, cols, x, y);
    initEulerflux2D(var, rows, cols, y,f,g);
}



#endif