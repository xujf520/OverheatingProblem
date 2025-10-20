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


// 边界配置类型
typedef enum {
    BC_OUTFLOW = 0,     // 出口边界
    BC_REFLECTION = 1,  // 固壁反射
    BC_PERIODICITY = 2, // 周期性
    BC_INFLOW = 3       // 入口边界
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

// 入口边界条件函数
static inline void BC_Inflow_2D(int var, int rows, int cols, double (*x)[rows][cols], const char* boundary_type, int Ghost_cell) {
    // 使用全局的 inflow_state 数组作为入流条件
    
    if (strcmp(boundary_type, L) == 0) {
        // 左边界入口
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
    for (int i = 0; i < var; i++){
        for (int k = Ghost_cell; k < cols-Ghost_cell; k++){
            x[i][rows-Ghost_cell][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+1][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+2][k] = x[i][rows-Ghost_cell-1][k];
            x[i][rows-Ghost_cell+3][k] = x[i][rows-Ghost_cell-1][k];
        }                   
    }
    
    // ========== 3. 下边界：部分反射壁面，部分激波后入流 ==========
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




static inline void Boundary_Conditions(int var, int rows, int cols, double (*y)[rows][cols], int GC) {
    // 边界条件函数指针数组（添加入口边界条件）
    void (*bc_funcs[])(int, int, int, double (*)[*][*], const char*, int) = {
        BC_OutFlow_2D,      // 0: BC_OUTFLOW
        BC_Reflection_2D,   // 1: BC_REFLECTION
        BC_Periodicity_2D,  // 2: BC_PERIODICITY
        BC_Inflow_2D        // 3: BC_INFLOW
    };
    
    const char* sides[] = {L, R, B, T};
    BC_Type bc_types[] = {bc_config.left, bc_config.right, bc_config.bottom, bc_config.top};
    
    for (int i = 0; i < 4; i++) {
        bc_funcs[bc_types[i]](var, rows, cols, y, sides[i], GC);
    }
}


#endif // BOUNDARY_CONDITION_H