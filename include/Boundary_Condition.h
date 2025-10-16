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
    BC_PERIODICITY = 2  // 周期性
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



//指定固定数值给边界条件
static inline void BC_In_rho(int var, int rows,double (*x)[rows], int Control_Pos, int Ghost_cell){
    switch (Control_Pos){

        //左边界反射边界条件设计！！！！！
        case 0:
            switch (Ghost_cell){
                case 1:                                      //此时填充虚拟网格为1
                    for (int i = 0; i <= 2; i++){
                        if (i == 1)//条件意思是速度反射
                            x[i][Ghost_cell-1] = -x[i][Ghost_cell];
                        else
                            x[i][Ghost_cell-1] = x[i][Ghost_cell];  
                    }
                    break;
                case 2:                                      //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        if (i == 1) {                        //条件意思是速度反射
                            x[i][Ghost_cell-1] = -x[i][Ghost_cell];
                            x[i][Ghost_cell-2] = -x[i][Ghost_cell+1];
                        }
                        else{
                            x[i][Ghost_cell-1] = x[i][Ghost_cell];
                            x[i][Ghost_cell-2] = x[i][Ghost_cell+1];
                        } 
                    }
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;
        //右边界反射边界条件设计！！！！！
        case 1:
            switch (Ghost_cell){
                case 1:                                        //此时填充虚拟网格为1
                    for (int i = 0; i <= 2; i++){
                        if (i == 1)//条件意思是速度反射
                            x[i][rows-Ghost_cell] = -x[i][rows-Ghost_cell-1];
                        else
                            x[i][rows-Ghost_cell] = x[i][rows-Ghost_cell-1];
                    }
                    break;

                case 2:                                      //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[0][rows-Ghost_cell] = 2.079156;
                        x[0][rows-Ghost_cell+1] = 2.079156;

                        x[1][rows-Ghost_cell] = 0.0;
                        x[1][rows-Ghost_cell+1] = 0.0;

                        x[2][rows-Ghost_cell] =  7.316625;
                        x[2][rows-Ghost_cell+1] =  7.316625;
                    }
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;
    }
}


static inline void Boundary_Conditions(int var, int rows, int cols, double (*y)[rows][cols], int GC) {
    // 边界条件函数指针数组
    void (*bc_funcs[])(int, int, int, double (*)[*][*], const char*, int) = {
        BC_OutFlow_2D,
        BC_Reflection_2D, 
        BC_Periodicity_2D
    };
    
    const char* sides[] = {L, R, B, T};
    BC_Type bc_types[] = {bc_config.left, bc_config.right, bc_config.bottom, bc_config.top};
    
    for (int i = 0; i < 4; i++) {
        bc_funcs[bc_types[i]](var, rows, cols, y, sides[i], GC);
    }
}


#endif // BOUNDARY_CONDITION_H