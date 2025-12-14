#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "scheme.h"
#include "function.h"

// 在头文件中定义边界类型宏
#define L "left"
#define R "right"  
#define T "top"
#define B "bottom"

//
typedef enum {
    BC_OUTFLOW = 0,     // 出口边界
    BC_REFLECTION = 1,  // 固壁反射
    BC_PERIODICITY = 2, // 周期性
    BC_INFLOW = 3,      // 入口边界
    BC_FIXED_VALUE = 4  // 固定值边界
} BC_Type;


// 边界配置结构
typedef struct {
    BC_Type left;
    BC_Type right;
    BC_Type bottom;
    BC_Type top;
} BoundaryConfig;

// 全局边界配置声明
extern BoundaryConfig bc_config;

// 全局入口状态声明（在某个头文件中）
extern double inflow_state[4]; // [密度, x动量, y动量, 压力]

// 全局固定值状态声明
extern double L_fixed_value_state[4]; 
extern double R_fixed_value_state[4]; 
extern double B_fixed_value_state[4]; 
extern double T_fixed_value_state[4]; 

// 固定值边界条件函数
static inline void BC_FixedValue_2D(int var, int rows, int cols, double (*x)[rows][cols], const char* boundary_type, int Ghost_cell) {
    double *fixed_state = NULL;
    
    // 根据边界类型选择对应的固定值状态
    if (strcmp(boundary_type, L) == 0) {
        fixed_state = L_fixed_value_state;
    }
    else if (strcmp(boundary_type, R) == 0) {
        fixed_state = R_fixed_value_state;
    }
    else if (strcmp(boundary_type, B) == 0) {
        fixed_state = B_fixed_value_state;
    }
    else if (strcmp(boundary_type, T) == 0) {
        fixed_state = T_fixed_value_state;
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
        return;
    }
    
    // 应用固定值边界条件
    if (strcmp(boundary_type, L) == 0) {
        // 左边界固定值
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][Ghost_cell-1][k] = fixed_state[i];
                x[i][Ghost_cell-2][k] = fixed_state[i];
                x[i][Ghost_cell-3][k] = fixed_state[i];
                x[i][Ghost_cell-4][k] = fixed_state[i];
            }
        }
    }
    else if (strcmp(boundary_type, R) == 0) {
        // 右边界固定值
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][rows - Ghost_cell][k] = fixed_state[i];
                x[i][rows - Ghost_cell + 1][k] = fixed_state[i];
                x[i][rows - Ghost_cell + 2][k] = fixed_state[i];
                x[i][rows - Ghost_cell + 3][k] = fixed_state[i];
            }
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // 下边界固定值
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][Ghost_cell-1] = fixed_state[i];
                x[i][j][Ghost_cell-2] = fixed_state[i];
                x[i][j][Ghost_cell-3] = fixed_state[i];
                x[i][j][Ghost_cell-4] = fixed_state[i];
            }
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // 上边界固定值
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][cols - Ghost_cell] = fixed_state[i];
                x[i][j][cols - Ghost_cell + 1] = fixed_state[i];
                x[i][j][cols - Ghost_cell + 2] = fixed_state[i];
                x[i][j][cols - Ghost_cell + 3] = fixed_state[i];
            }
        }
    }
}


// 入口边界条件函数
static inline void BC_Inflow_2D(int var, int rows, int cols, double (*x)[rows][cols], const char* boundary_type, int Ghost_cell) {
    // 使用全局的 inflow_state 数组作为入流条件
    if (strcmp(boundary_type, L) == 0) {
        // 左边界入口
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][Ghost_cell-1][k] = inflow_state[i];
                x[i][Ghost_cell-2][k] = inflow_state[i];
                x[i][Ghost_cell-3][k] = inflow_state[i];
                x[i][Ghost_cell-4][k] = inflow_state[i];
            }
        }
    }
    else if (strcmp(boundary_type, R) == 0) {
        // 右边界入口
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++) {
                x[i][rows - Ghost_cell][k] = inflow_state[i];
                x[i][rows - Ghost_cell + 1][k] = inflow_state[i];
                x[i][rows - Ghost_cell + 2][k] = inflow_state[i];
                x[i][rows - Ghost_cell + 3][k] = inflow_state[i];
            }
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // 下边界入口
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][Ghost_cell-1] = inflow_state[i];
                x[i][j][Ghost_cell-2] = inflow_state[i];
                x[i][j][Ghost_cell-3] = inflow_state[i];
                x[i][j][Ghost_cell-4] = inflow_state[i];
            }
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // 上边界入口
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++) {
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++) {
                x[i][j][cols - Ghost_cell] = inflow_state[i];
                x[i][j][cols - Ghost_cell + 1] = inflow_state[i];
                x[i][j][cols - Ghost_cell + 2] = inflow_state[i];
                x[i][j][cols - Ghost_cell + 3] = inflow_state[i];
            }
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}


static inline void BC_OutFlow_2D(int var, int rows, int cols, double (*x)[rows][cols], const char* boundary_type, int Ghost_cell){
    // 根据边界类型字符串执行相应的边界条件
    if (strcmp(boundary_type, L) == 0) {
        //左边界出口边界条件设计
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k <= cols-Ghost_cell; k++){
                x[i][Ghost_cell-1][k] = x[i][Ghost_cell][k];
                x[i][Ghost_cell-2][k] = x[i][Ghost_cell][k];
                x[i][Ghost_cell-3][k] = x[i][Ghost_cell][k];
                x[i][Ghost_cell-4][k] = x[i][Ghost_cell][k];
            }
        }   
    }
    else if (strcmp(boundary_type, R) == 0) {
        //右边界无反射边界条件设计
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k <= cols-Ghost_cell; k++){
                x[i][rows-Ghost_cell][k] = x[i][rows-Ghost_cell-1][k];
                x[i][rows-Ghost_cell+1][k] = x[i][rows-Ghost_cell-1][k];
                x[i][rows-Ghost_cell+2][k] = x[i][rows-Ghost_cell-1][k];
                x[i][rows-Ghost_cell+3][k] = x[i][rows-Ghost_cell-1][k];
            }                   
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        //下边界
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j <= rows-Ghost_cell; j++){
                x[i][j][Ghost_cell-1] = x[i][j][Ghost_cell];
                x[i][j][Ghost_cell-2] = x[i][j][Ghost_cell];
                x[i][j][Ghost_cell-3] = x[i][j][Ghost_cell];
                x[i][j][Ghost_cell-4] = x[i][j][Ghost_cell];
            }                   
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        //上边界
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j <= rows-Ghost_cell; j++){
                x[i][j][cols-Ghost_cell] = x[i][j][cols-Ghost_cell-1];
                x[i][j][cols-Ghost_cell+1] = x[i][j][cols-Ghost_cell-1];
                x[i][j][cols-Ghost_cell+2] = x[i][j][cols-Ghost_cell-1];
                x[i][j][cols-Ghost_cell+3] = x[i][j][cols-Ghost_cell-1];
            }                   
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}


//周期性边界条件设计
static inline void BC_Periodicity_2D(int var, int rows, int cols, double (*x)[rows][cols], const char* boundary_type, int Ghost_cell){
    if (strcmp(boundary_type, L) == 0) {
        // 左边界 <- 右边界内部区域
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){  
                x[i][Ghost_cell-1][k] = x[i][rows - Ghost_cell - 1][k];  // 左ghost <- 右内部
                x[i][Ghost_cell-2][k] = x[i][rows - Ghost_cell - 2][k];
                x[i][Ghost_cell-3][k] = x[i][rows - Ghost_cell - 3][k];
                x[i][Ghost_cell-4][k] = x[i][rows - Ghost_cell - 4][k];
            }
        }   
    }
    else if (strcmp(boundary_type, R) == 0) {
        // 右边界 <- 左边界内部区域
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){
                x[i][rows - Ghost_cell][k] = x[i][Ghost_cell][k];         // 右ghost <- 左内部
                x[i][rows - Ghost_cell + 1][k] = x[i][Ghost_cell + 1][k];
                x[i][rows - Ghost_cell + 2][k] = x[i][Ghost_cell + 2][k];
                x[i][rows - Ghost_cell + 3][k] = x[i][Ghost_cell + 3][k];
            }                   
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // 下边界 <- 上边界内部区域
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                x[i][j][Ghost_cell-1] = x[i][j][cols - Ghost_cell - 1];   // 下ghost <- 上内部
                x[i][j][Ghost_cell-2] = x[i][j][cols - Ghost_cell - 2];
                x[i][j][Ghost_cell-3] = x[i][j][cols - Ghost_cell - 3];
                x[i][j][Ghost_cell-4] = x[i][j][cols - Ghost_cell - 4];
            }                   
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // 上边界 <- 下边界内部区域
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                x[i][j][cols - Ghost_cell] = x[i][j][Ghost_cell];         // 上ghost <- 下内部
                x[i][j][cols - Ghost_cell + 1] = x[i][j][Ghost_cell + 1];
                x[i][j][cols - Ghost_cell + 2] = x[i][j][Ghost_cell + 2];
                x[i][j][cols - Ghost_cell + 3] = x[i][j][Ghost_cell + 3];
            }                   
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}



static inline void BC_Reflection_2D(int var, int rows, int cols, double (*x)[rows][cols], const char* boundary_type, int Ghost_cell){
    // 假设变量顺序为: [密度, x方向动量, y方向动量, 能量, ...]
    // 对于反射边界，法向速度反向，切向速度不变，其他变量不变
    if (strcmp(boundary_type, L) == 0) {
        // 左边界固壁反射：x方向速度反向，y方向速度不变
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){
                if (i == 1) { // x方向动量（速度反向）
                    x[i][Ghost_cell-1][k] = -x[i][Ghost_cell][k];
                    x[i][Ghost_cell-2][k] = -x[i][Ghost_cell+1][k];
                    x[i][Ghost_cell-3][k] = -x[i][Ghost_cell+2][k];
                    x[i][Ghost_cell-4][k] = -x[i][Ghost_cell+3][k];
                } else { // 密度、y方向动量、能量等保持不变
                    x[i][Ghost_cell-1][k] = x[i][Ghost_cell][k];
                    x[i][Ghost_cell-2][k] = x[i][Ghost_cell+1][k];
                    x[i][Ghost_cell-3][k] = x[i][Ghost_cell+2][k];
                    x[i][Ghost_cell-4][k] = x[i][Ghost_cell+3][k];
                }
            }
        }   
    }
    else if (strcmp(boundary_type, R) == 0) {
        // 右边界固壁反射：x方向速度反向，y方向速度不变
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int k = Ghost_cell; k < cols - Ghost_cell; k++){
                if (i == 1) { // x方向动量（速度反向）
                    x[i][rows - Ghost_cell][k] = -x[i][rows - Ghost_cell - 1][k];
                    x[i][rows - Ghost_cell + 1][k] = -x[i][rows - Ghost_cell - 2][k];
                    x[i][rows - Ghost_cell + 2][k] = -x[i][rows - Ghost_cell - 3][k];
                    x[i][rows - Ghost_cell + 3][k] = -x[i][rows - Ghost_cell - 4][k];
                } else { // 密度、y方向动量、能量等保持不变
                    x[i][rows - Ghost_cell][k] = x[i][rows - Ghost_cell - 1][k];
                    x[i][rows - Ghost_cell + 1][k] = x[i][rows - Ghost_cell - 2][k];
                    x[i][rows - Ghost_cell + 2][k] = x[i][rows - Ghost_cell - 3][k];
                    x[i][rows - Ghost_cell + 3][k] = x[i][rows - Ghost_cell - 4][k];
                }
            }                   
        }
    }
    else if (strcmp(boundary_type, B) == 0) {
        // 下边界固壁反射：y方向速度反向，x方向速度不变
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                if (i == 2) { // y方向动量（速度反向）
                    x[i][j][Ghost_cell-1] = -x[i][j][Ghost_cell];
                    x[i][j][Ghost_cell-2] = -x[i][j][Ghost_cell+1];
                    x[i][j][Ghost_cell-3] = -x[i][j][Ghost_cell+2];
                    x[i][j][Ghost_cell-4] = -x[i][j][Ghost_cell+3];
                } else { // 密度、x方向动量、能量等保持不变
                    x[i][j][Ghost_cell-1] = x[i][j][Ghost_cell];
                    x[i][j][Ghost_cell-2] = x[i][j][Ghost_cell+1];
                    x[i][j][Ghost_cell-3] = x[i][j][Ghost_cell+2];
                    x[i][j][Ghost_cell-4] = x[i][j][Ghost_cell+3];
                }
            }                   
        }
    }
    else if (strcmp(boundary_type, T) == 0) {
        // 上边界固壁反射：y方向速度反向，x方向速度不变
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < var; i++){
            for (int j = Ghost_cell; j < rows - Ghost_cell; j++){
                if (i == 2) { // y方向动量（速度反向）
                    x[i][j][cols - Ghost_cell] = -x[i][j][cols - Ghost_cell - 1];
                    x[i][j][cols - Ghost_cell + 1] = -x[i][j][cols - Ghost_cell - 2];
                    x[i][j][cols - Ghost_cell + 2] = -x[i][j][cols - Ghost_cell - 3];
                    x[i][j][cols - Ghost_cell + 3] = -x[i][j][cols - Ghost_cell - 4];
                } else { // 密度、x方向动量、能量等保持不变
                    x[i][j][cols - Ghost_cell] = x[i][j][cols - Ghost_cell - 1];
                    x[i][j][cols - Ghost_cell + 1] = x[i][j][cols - Ghost_cell - 2];
                    x[i][j][cols - Ghost_cell + 2] = x[i][j][cols - Ghost_cell - 3];
                    x[i][j][cols - Ghost_cell + 3] = x[i][j][cols - Ghost_cell - 4];
                }
            }                   
        }
    }
    else {
        printf("Warning: Unknown boundary type: %s\n", boundary_type);
    }
}

// 双马赫反射问题的边界条件函数 - 一次调用处理所有边界
static inline void BC_DoubleMach_2D(int var, int rows, int cols, double (*x)[rows][cols], int Ghost_cell, double current_time) {
    // 双马赫反射问题的物理参数
    double Lx = 4.0;
    double Ly = 1.0;
    double shock_slope = 1.732;  // tan(60°)
    double shock_start_x = 1.0/6.0;
    double wall_start_x = 1.0/6.0;
    double gamma = 1.4;  // 比热比
    
    // 激波前状态 - 原始变量
    double pre_shock_rho = 1.4;
    double pre_shock_u = 0.0;
    double pre_shock_v = 0.0;
    double pre_shock_p = 1.0;
    
    // 激波后状态 - 原始变量
    double post_shock_rho = 8.0;
    double post_shock_u = 7.145;
    double post_shock_v = -4.125;
    double post_shock_p = 116.83333;
    
    // 将原始变量转换为守恒变量
    // 激波前守恒变量
    double pre_rho = pre_shock_rho;
    double pre_rhou = pre_shock_rho * pre_shock_u;
    double pre_rhov = pre_shock_rho * pre_shock_v;
    double pre_E = pre_shock_p/(gamma-1.0) + 0.5*pre_shock_rho*(pre_shock_u*pre_shock_u + pre_shock_v*pre_shock_v);
    
    // 激波后守恒变量
    double post_rho = post_shock_rho;
    double post_rhou = post_shock_rho * post_shock_u;
    double post_rhov = post_shock_rho * post_shock_v;
    double post_E = post_shock_p/(gamma-1.0) + 0.5*post_shock_rho*(post_shock_u*post_shock_u + post_shock_v*post_shock_v);
    
    double nx = rows - 2 * Ghost_cell;
    
    // ========== 1. 左边界：激波后入流 ==========
    #pragma omp parallel for collapse(1)
    for (int k = Ghost_cell; k < cols-Ghost_cell; k++){
        // 密度
        x[0][Ghost_cell-1][k] = post_rho;
        x[0][Ghost_cell-2][k] = post_rho;
        x[0][Ghost_cell-3][k] = post_rho;
        x[0][Ghost_cell-4][k] = post_rho;
        
        // x动量
        x[1][Ghost_cell-1][k] = post_rhou;
        x[1][Ghost_cell-2][k] = post_rhou;
        x[1][Ghost_cell-3][k] = post_rhou;
        x[1][Ghost_cell-4][k] = post_rhou;
        
        // y动量
        x[2][Ghost_cell-1][k] = post_rhov;
        x[2][Ghost_cell-2][k] = post_rhov;
        x[2][Ghost_cell-3][k] = post_rhov;
        x[2][Ghost_cell-4][k] = post_rhov;
        
        // 总能量
        x[3][Ghost_cell-1][k] = post_E;
        x[3][Ghost_cell-2][k] = post_E;
        x[3][Ghost_cell-3][k] = post_E;
        x[3][Ghost_cell-4][k] = post_E;
    }
    
    // ========== 2. 右边界：出流边界 ==========
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < var; i++){
        for (int k = Ghost_cell; k < cols-Ghost_cell; k++){
            x[i][rows-Ghost_cell][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+1][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+2][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+3][k] = x[i][rows-Ghost_cell-1][k];
        }                   
    }
    
    // ========== 3. 下边界：部分反射壁面，部分激波后入流 ==========
    #pragma omp parallel for collapse(1)
    for (int j = Ghost_cell; j < rows-Ghost_cell; j++){
        // 计算物理坐标
        double x_pos = (double)(j+0.5-Ghost_cell)/ nx * Lx;
        if (x_pos >= wall_start_x) {
            // 反射壁面区域
            x[0][j][Ghost_cell-1] = x[0][j][Ghost_cell];
            x[0][j][Ghost_cell-2] = x[0][j][Ghost_cell+1];
            x[0][j][Ghost_cell-3] = x[0][j][Ghost_cell+2];
            x[0][j][Ghost_cell-4] = x[0][j][Ghost_cell+3];
            
            // x动量 - 保持不变
            x[1][j][Ghost_cell-1] = x[1][j][Ghost_cell];
            x[1][j][Ghost_cell-2] = x[1][j][Ghost_cell+1];
            x[1][j][Ghost_cell-3] = x[1][j][Ghost_cell+2];
            x[1][j][Ghost_cell-4] = x[1][j][Ghost_cell+3];
            
            // y动量 - 速度反向（注意：动量是ρv，所以直接取负）
            x[2][j][Ghost_cell-1] = -x[2][j][Ghost_cell];
            x[2][j][Ghost_cell-2] = -x[2][j][Ghost_cell+1];
            x[2][j][Ghost_cell-3] = -x[2][j][Ghost_cell+2];
            x[2][j][Ghost_cell-4] = -x[2][j][Ghost_cell+3];
            
            // 总能量 - 保持不变
            x[3][j][Ghost_cell-1] = x[3][j][Ghost_cell];
            x[3][j][Ghost_cell-2] = x[3][j][Ghost_cell+1];
            x[3][j][Ghost_cell-3] = x[3][j][Ghost_cell+2];
            x[3][j][Ghost_cell-4] = x[3][j][Ghost_cell+3];
            
        } else {
            // 非壁面区域 - 激波后入流
           
            // 密度
            x[0][j][Ghost_cell-1] = post_rho;
            x[0][j][Ghost_cell-2] = post_rho;
            x[0][j][Ghost_cell-3] = post_rho;
            x[0][j][Ghost_cell-4] = post_rho;
            
            // x动量
            x[1][j][Ghost_cell-1] = post_rhou;
            x[1][j][Ghost_cell-2] = post_rhou;
            x[1][j][Ghost_cell-3] = post_rhou;
            x[1][j][Ghost_cell-4] = post_rhou;
            
            // y动量
            x[2][j][Ghost_cell-1] = post_rhov;
            x[2][j][Ghost_cell-2] = post_rhov;
            x[2][j][Ghost_cell-3] = post_rhov;
            x[2][j][Ghost_cell-4] = post_rhov;
            
            // 总能量
            x[3][j][Ghost_cell-1] = post_E;
            x[3][j][Ghost_cell-2] = post_E;
            x[3][j][Ghost_cell-3] = post_E;
            x[3][j][Ghost_cell-4] = post_E;
        }
    }
    
    // ========== 4. 上边界：根据激波位置设置 ==========
// 计算激波与上边界的交点
    double shock_x_at_top = shock_start_x + (Ly / shock_slope) + (20.0 * current_time / shock_slope);
    
    #pragma omp parallel for collapse(1)
    for (int j = Ghost_cell; j < rows-Ghost_cell; j++){
        // 计算物理坐标 - 使用网格中心坐标
        double x_pos = (double)(j + 0.5 - Ghost_cell) / nx * Lx;
        if (x_pos < shock_x_at_top) {
            // 激波后区域
                // 密度
                x[0][j][cols-Ghost_cell] = post_rho;
                x[0][j][cols-Ghost_cell+1] = post_rho;
                x[0][j][cols-Ghost_cell+2] = post_rho;
                x[0][j][cols-Ghost_cell+3] = post_rho;
                
                // x动量
                x[1][j][cols-Ghost_cell] = post_rhou;
                x[1][j][cols-Ghost_cell+1] = post_rhou;
                x[1][j][cols-Ghost_cell+2] = post_rhou;
                x[1][j][cols-Ghost_cell+3] = post_rhou;
                
                // y动量
                x[2][j][cols-Ghost_cell] = post_rhov;
                x[2][j][cols-Ghost_cell+1] = post_rhov;
                x[2][j][cols-Ghost_cell+2] = post_rhov;
                x[2][j][cols-Ghost_cell+3] = post_rhov;
                
                // 总能量
                x[3][j][cols-Ghost_cell] = post_E;
                x[3][j][cols-Ghost_cell+1] = post_E;
                x[3][j][cols-Ghost_cell+2] = post_E;
                x[3][j][cols-Ghost_cell+3] = post_E;
            
        } else {
            // 激波前区域
            
                // 密度
                x[0][j][cols-Ghost_cell] = pre_rho;
                x[0][j][cols-Ghost_cell+1] = pre_rho;
                x[0][j][cols-Ghost_cell+2] = pre_rho;
                x[0][j][cols-Ghost_cell+3] = pre_rho;
                
                // x动量
                x[1][j][cols-Ghost_cell] = pre_rhou;
                x[1][j][cols-Ghost_cell+1] = pre_rhou;
                x[1][j][cols-Ghost_cell+2] = pre_rhou;
                x[1][j][cols-Ghost_cell+3] = pre_rhou;
                
                // y动量
                x[2][j][cols-Ghost_cell] = pre_rhov;
                x[2][j][cols-Ghost_cell+1] = pre_rhov;
                x[2][j][cols-Ghost_cell+2] = pre_rhov;
                x[2][j][cols-Ghost_cell+3] = pre_rhov;
                
                // 总能量
                x[3][j][cols-Ghost_cell] = pre_E;
                x[3][j][cols-Ghost_cell+1] = pre_E;
                x[3][j][cols-Ghost_cell+2] = pre_E;
                x[3][j][cols-Ghost_cell+3] = pre_E;
         
        }
    }
}

// 二维后台阶流动的边界条件函数 - 根据Woodward & Colella (1984)
static inline void BC_BackwardStep_2D(int var, int rows, int cols, double (*x)[rows][cols], int Ghost_cell, double current_time) {
    // 后台阶流动的物理参数
    double Lx = 3.0;
    double Ly = 1.0;
    double gamma = 1.4;  // 比热比
    
    // 台阶几何参数
    double step_height = 0.2;     // 台阶高度
    double step_position = 0.6;   // 台阶位置
    
    // 来流条件（马赫数3超声速流动）
    double inflow_density = 1.4;      // 来流密度
    double inflow_pressure = 1.0;     // 来流压力
    double sound_speed = sqrt(gamma * inflow_pressure / inflow_density);
    double inflow_velocity_x = 3.0 * sound_speed;  // 马赫数3
    double inflow_velocity_y = 0.0;
    
    // 计算守恒变量
    double inflow_rho = inflow_density;
    double inflow_rhou = inflow_density * inflow_velocity_x;
    double inflow_rhov = 0.0;
    double inflow_E = inflow_pressure/(gamma-1.0) + 
                      0.5 * inflow_density * inflow_velocity_x * inflow_velocity_x;
    
    double nx = rows - 2 * Ghost_cell;
    double ny = cols - 2 * Ghost_cell;
    double dx = Lx / nx;
    double dy = Ly / ny;
    
    // 计算台阶对应的网格索引
    int step_idx = (int)(step_position / dx) + Ghost_cell;
    int step_top_idx = Ghost_cell + (int)(step_height / dy);
    
    // ========== 1. 左边界：超声速入流 ==========
    // 整个左边界都是固定来流条件
    for (int k = Ghost_cell; k < cols-Ghost_cell; k++) {
        // 密度
        x[0][Ghost_cell-1][k] = inflow_rho;
        x[0][Ghost_cell-2][k] = inflow_rho;
        x[0][Ghost_cell-3][k] = inflow_rho;
        x[0][Ghost_cell-4][k] = inflow_rho;
        
        // x动量
        x[1][Ghost_cell-1][k] = inflow_rhou;
        x[1][Ghost_cell-2][k] = inflow_rhou;
        x[1][Ghost_cell-3][k] = inflow_rhou;
        x[1][Ghost_cell-4][k] = inflow_rhou;
        
        // y动量
        x[2][Ghost_cell-1][k] = inflow_rhov;
        x[2][Ghost_cell-2][k] = inflow_rhov;
        x[2][Ghost_cell-3][k] = inflow_rhov;
        x[2][Ghost_cell-4][k] = inflow_rhov;
        
        // 总能量
        x[3][Ghost_cell-1][k] = inflow_E;
        x[3][Ghost_cell-2][k] = inflow_E;
        x[3][Ghost_cell-3][k] = inflow_E;
        x[3][Ghost_cell-4][k] = inflow_E;
    }
    
    // ========== 2. 右边界：超声速出流 ==========
    // 零梯度外推
    for (int i = 0; i < var; i++) {
        for (int k = Ghost_cell; k < cols-Ghost_cell; k++) {
            x[i][rows-Ghost_cell][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+1][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+2][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+3][k] = x[i][rows-Ghost_cell-1][k];
        }
    }
    
    // ========== 3. 下边界：反射壁面 ==========
    // 整个下边界都是反射壁面（包括台阶前和台阶后）
    for (int j = Ghost_cell; j < rows-Ghost_cell; j++) {
        // 计算物理坐标
        double x_pos = (double)(j + 0.5 - Ghost_cell) * dx;
        
        // 标准反射壁面条件
        // 密度
        x[0][j][Ghost_cell-1] = x[0][j][Ghost_cell];
        x[0][j][Ghost_cell-2] = x[0][j][Ghost_cell+1];
        x[0][j][Ghost_cell-3] = x[0][j][Ghost_cell+2];
        x[0][j][Ghost_cell-4] = x[0][j][Ghost_cell+3];
        
        // x动量 - 切向速度不变
        x[1][j][Ghost_cell-1] = x[1][j][Ghost_cell];
        x[1][j][Ghost_cell-2] = x[1][j][Ghost_cell+1];
        x[1][j][Ghost_cell-3] = x[1][j][Ghost_cell+2];
        x[1][j][Ghost_cell-4] = x[1][j][Ghost_cell+3];
        
        // y动量 - 法向速度反向
        x[2][j][Ghost_cell-1] = -x[2][j][Ghost_cell];
        x[2][j][Ghost_cell-2] = -x[2][j][Ghost_cell+1];
        x[2][j][Ghost_cell-3] = -x[2][j][Ghost_cell+2];
        x[2][j][Ghost_cell-4] = -x[2][j][Ghost_cell+3];
        
        // 总能量
        x[3][j][Ghost_cell-1] = x[3][j][Ghost_cell];
        x[3][j][Ghost_cell-2] = x[3][j][Ghost_cell+1];
        x[3][j][Ghost_cell-3] = x[3][j][Ghost_cell+2];
        x[3][j][Ghost_cell-4] = x[3][j][Ghost_cell+3];
    }
    
    // ========== 4. 上边界：反射壁面 ==========
    // 整个上边界都是反射壁面
    for (int j = Ghost_cell; j < rows-Ghost_cell; j++) {
        // 密度
        x[0][j][cols-Ghost_cell] = x[0][j][cols-Ghost_cell-1];
        x[0][j][cols-Ghost_cell+1] = x[0][j][cols-Ghost_cell-2];
        x[0][j][cols-Ghost_cell+2] = x[0][j][cols-Ghost_cell-3];
        x[0][j][cols-Ghost_cell+3] = x[0][j][cols-Ghost_cell-4];
        
        // x动量 - 切向速度不变
        x[1][j][cols-Ghost_cell] = x[1][j][cols-Ghost_cell-1];
        x[1][j][cols-Ghost_cell+1] = x[1][j][cols-Ghost_cell-2];
        x[1][j][cols-Ghost_cell+2] = x[1][j][cols-Ghost_cell-3];
        x[1][j][cols-Ghost_cell+3] = x[1][j][cols-Ghost_cell-4];
        
        // y动量 - 法向速度反向
        x[2][j][cols-Ghost_cell] = -x[2][j][cols-Ghost_cell-1];
        x[2][j][cols-Ghost_cell+1] = -x[2][j][cols-Ghost_cell-2];
        x[2][j][cols-Ghost_cell+2] = -x[2][j][cols-Ghost_cell-3];
        x[2][j][cols-Ghost_cell+3] = -x[2][j][cols-Ghost_cell-4];
        
        // 总能量
        x[3][j][cols-Ghost_cell] = x[3][j][cols-Ghost_cell-1];
        x[3][j][cols-Ghost_cell+1] = x[3][j][cols-Ghost_cell-2];
        x[3][j][cols-Ghost_cell+2] = x[3][j][cols-Ghost_cell-3];
        x[3][j][cols-Ghost_cell+3] = x[3][j][cols-Ghost_cell-4];
    }
    
    // ========== 5. 台阶垂直面边界条件 ==========
    // 处理x = step_position处的垂直壁面
    // 注意：这里处理的是台阶右侧的幽灵单元格
    if (step_idx >= Ghost_cell && step_idx < rows) {
        for (int k = Ghost_cell; k < step_top_idx; k++) {
            if (k >= Ghost_cell && k < cols) {
                // 垂直壁面：法向为x方向
                // 使用镜像反射边界条件
                
                // 第1层幽灵单元格
                x[0][step_idx][k] = x[0][step_idx-1][k];      // 密度外推
                x[1][step_idx][k] = -x[1][step_idx-1][k];     // x动量反射
                x[2][step_idx][k] = x[2][step_idx-1][k];      // y动量外推
                x[3][step_idx][k] = x[3][step_idx-1][k];      // 能量外推
                
                // 第2层幽灵单元格
                if (step_idx+1 < rows) {
                    x[0][step_idx+1][k] = x[0][step_idx-2][k];
                    x[1][step_idx+1][k] = -x[1][step_idx-2][k];
                    x[2][step_idx+1][k] = x[2][step_idx-2][k];
                    x[3][step_idx+1][k] = x[3][step_idx-2][k];
                }
                
                // 第3层幽灵单元格
                if (step_idx+2 < rows) {
                    x[0][step_idx+2][k] = x[0][step_idx-3][k];
                    x[1][step_idx+2][k] = -x[1][step_idx-3][k];
                    x[2][step_idx+2][k] = x[2][step_idx-3][k];
                    x[3][step_idx+2][k] = x[3][step_idx-3][k];
                }
                
                // 第4层幽灵单元格
                if (step_idx+3 < rows) {
                    x[0][step_idx+3][k] = x[0][step_idx-4][k];
                    x[1][step_idx+3][k] = -x[1][step_idx-4][k];
                    x[2][step_idx+3][k] = x[2][step_idx-4][k];
                    x[3][step_idx+3][k] = x[3][step_idx-4][k];
                }
            }
        }
    }
    
    // ========== 6. 台阶上缘边界条件 ==========
    // 处理y = step_height处的水平壁面
    // 注意：这里处理的是台阶上方的幽灵单元格
    if (step_top_idx >= Ghost_cell && step_top_idx < cols) {
        for (int j = step_idx; j < rows-Ghost_cell; j++) {
            if (j >= Ghost_cell && j < rows) {
                // 水平壁面：法向为y方向
                
                // 第1层幽灵单元格
                x[0][j][step_top_idx] = x[0][j][step_top_idx+1];      // 密度外推
                x[1][j][step_top_idx] = x[1][j][step_top_idx+1];      // x动量外推
                x[2][j][step_top_idx] = -x[2][j][step_top_idx+1];     // y动量反射
                x[3][j][step_top_idx] = x[3][j][step_top_idx+1];      // 能量外推
                
                // 第2层幽灵单元格
                if (step_top_idx-1 >= 0) {
                    x[0][j][step_top_idx-1] = x[0][j][step_top_idx+2];
                    x[1][j][step_top_idx-1] = x[1][j][step_top_idx+2];
                    x[2][j][step_top_idx-1] = -x[2][j][step_top_idx+2];
                    x[3][j][step_top_idx-1] = x[3][j][step_top_idx+2];
                }
                
                // 第3层幽灵单元格
                if (step_top_idx-2 >= 0) {
                    x[0][j][step_top_idx-2] = x[0][j][step_top_idx+3];
                    x[1][j][step_top_idx-2] = x[1][j][step_top_idx+3];
                    x[2][j][step_top_idx-2] = -x[2][j][step_top_idx+3];
                    x[3][j][step_top_idx-2] = x[3][j][step_top_idx+3];
                }
                
                // 第4层幽灵单元格
                if (step_top_idx-3 >= 0) {
                    x[0][j][step_top_idx-3] = x[0][j][step_top_idx+4];
                    x[1][j][step_top_idx-3] = x[1][j][step_top_idx+4];
                    x[2][j][step_top_idx-3] = -x[2][j][step_top_idx+4];
                    x[3][j][step_top_idx-3] = x[3][j][step_top_idx+4];
                }
            }
        }
    }
    
    // ========== 7. 台阶拐角特殊处理 ==========
    // 处理x=step_position, y=step_height拐角处的幽灵单元格
    if (step_idx >= Ghost_cell && step_idx < rows && 
        step_top_idx >= Ghost_cell && step_top_idx < cols) {
        
        // 拐角点本身（如果是幽灵单元格）
        if (step_idx < rows-Ghost_cell || step_top_idx < cols-Ghost_cell) {
            // 使用对角线反射：两个方向都反射
            x[0][step_idx][step_top_idx] = x[0][step_idx-1][step_top_idx+1];
            x[1][step_idx][step_top_idx] = -x[1][step_idx-1][step_top_idx+1];  // x反射
            x[2][step_idx][step_top_idx] = -x[2][step_idx-1][step_top_idx+1];  // y反射
            x[3][step_idx][step_top_idx] = x[3][step_idx-1][step_top_idx+1];
        }
        
        // 拐角附近的幽灵单元格
        for (int layer = 1; layer < 4; layer++) {
            // 垂直方向
            if (step_idx+layer < rows && step_top_idx < cols) {
                x[0][step_idx+layer][step_top_idx] = x[0][step_idx-1][step_top_idx+1];
                x[1][step_idx+layer][step_top_idx] = -x[1][step_idx-1][step_top_idx+1];
                x[2][step_idx+layer][step_top_idx] = x[2][step_idx-1][step_top_idx+1];
                x[3][step_idx+layer][step_top_idx] = x[3][step_idx-1][step_top_idx+1];
            }
            
            // 水平方向
            if (step_idx < rows && step_top_idx-layer >= 0) {
                x[0][step_idx][step_top_idx-layer] = x[0][step_idx-1][step_top_idx+1];
                x[1][step_idx][step_top_idx-layer] = x[1][step_idx-1][step_top_idx+1];
                x[2][step_idx][step_top_idx-layer] = -x[2][step_idx-1][step_top_idx+1];
                x[3][step_idx][step_top_idx-layer] = x[3][step_idx-1][step_top_idx+1];
            }
        }
    }
}



// 更新边界条件函数指针数组
static inline void Boundary_Conditions(int var, int rows, int cols, double (*y)[rows][cols], int GC) {
    // 边界条件函数指针数组（添加固定值边界条件）
    void (*bc_funcs[])(int, int, int, double (*)[*][*], const char*, int) = {
        BC_OutFlow_2D,      // 0: BC_OUTFLOW
        BC_Reflection_2D,   // 1: BC_REFLECTION
        BC_Periodicity_2D,  // 2: BC_PERIODICITY
        BC_Inflow_2D,       // 3: BC_INFLOW
        BC_FixedValue_2D    // 4: BC_FIXED_VALUE
    };
    
    const char* sides[] = {L, R, B, T};
    BC_Type bc_types[] = {bc_config.left, bc_config.right, bc_config.bottom, bc_config.top};
    
    for (int i = 0; i < 4; i++) {
        bc_funcs[bc_types[i]](var, rows, cols, y, sides[i], GC);
    }
}


#endif // BOUNDARY_CONDITION_H