#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scheme.h"
#include "function.h"

static inline void BC_Reflect(int rows, int cols,double (*x)[cols], int Control_Pos, int Ghost_cell){
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
                            x[i][cols-Ghost_cell] = -x[i][cols-Ghost_cell-1];
                        else
                            x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];
                    }
                    break;

                case 2:                                      //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        if (i == 1) {                        //条件意思是速度反射
                            x[i][cols-Ghost_cell] = -x[i][cols-Ghost_cell-1];
                            x[i][cols-Ghost_cell+1] = -x[i][cols-Ghost_cell-2];
                        }
                        else{
                            x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];
                            x[i][cols-Ghost_cell+1] = x[i][cols-Ghost_cell-2];
                        } 
                    }
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;
    }
}


//未完成：入口边界需要指定
static inline void BC_InFlow(int rows, int cols,double (*x)[cols], int Control_Pos, int Ghost_cell){
    switch (Control_Pos){

        //左边界反射边界条件设计！！！！！
        case 0:
            switch (Ghost_cell){
                //此时填充虚拟网格为1
                case 1:
                    for (int i = 0; i <= 2; i++){
                        if (i == 1)//条件意思是速度反射
                            x[i][0] = -x[i][1];
                        else
                            x[i][0] = x[i][1];  
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
                //此时填充虚拟网格为1
                case 1:
                    for (int i = 0; i <= 2; i++){
                        if (i == 1)//条件意思是速度反射
                            x[i][cols-Ghost_cell] = -x[i][cols-Ghost_cell-1];
                        else
                            x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];
                    }
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;
    }
}


//出口边界条件设计
static inline void BC_OutFlow(int rows, int cols,double (*x)[cols], int Control_Pos, int Ghost_cell){
    switch (Control_Pos){

        //左边界无反射边界条件设计！！！！！
        case 0:
            switch (Ghost_cell){
                case 1:                                                                 //此时填充虚拟网格为1
                    for (int i = 0; i <= 2; i++)
                        x[i][0] = x[i][1];                                              //出口流动的情况
                    break;
                case 2:                                                                 //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[i][Ghost_cell-1] = x[i][Ghost_cell];                           //出口流动的情况
                        x[i][Ghost_cell-2] = x[i][Ghost_cell];                           //出口流动的情况
                    }
                       
                    break;
                case 3:                                                                 //此时填充虚拟网格为3
                    for (int i = 0; i <= 2; i++){
                        x[i][Ghost_cell-1] = x[i][Ghost_cell];                           //出口流动的情况
                        x[i][Ghost_cell-2] = x[i][Ghost_cell];                           //出口流动的情况
                        x[i][Ghost_cell-3] = x[i][Ghost_cell];                           //出口流动的情况
                    }
                       
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;
        //右边界无反射边界条件设计！！！！！
        case 1:
            switch (Ghost_cell){
                case 1:                                                                 //此时填充虚拟网格为1
                    for (int i = 0; i <= 2; i++)
                        x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];                //出口流动的情况
                    break;
                case 2:                                                                 //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];                           //出口流动的情况
                        x[i][cols-Ghost_cell+1] = x[i][cols-Ghost_cell-1];                           //出口流动的情况
                    }
                    break;
                case 3:                                                                 //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];                           //出口流动的情况
                        x[i][cols-Ghost_cell+1] = x[i][cols-Ghost_cell-1];                           //出口流动的情况
                        x[i][cols-Ghost_cell+2] = x[i][cols-Ghost_cell-1];                           //出口流动的情况
                    }
                    break;
                default:
                    break;
            }
            break;
    }
}


//出口边界条件设计
static inline void BC_Periodicity(int rows, int cols,double (*x)[cols], int Control_Pos, int Ghost_cell){
    switch (Control_Pos){

        //左边界周期性！！！！！
        case 0:
            switch (Ghost_cell){
                case 1:                                                                 //此时填充虚拟网格为1
                    for (int i = 0; i <= 2; i++)
                        x[i][Ghost_cell-1] = x[i][cols-Ghost_cell-1];                                      
                    break;
                case 2:                                                                 //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[i][Ghost_cell-1] = x[i][cols-Ghost_cell-1];             
                        x[i][Ghost_cell-2] = x[i][cols-Ghost_cell-2];             
                    }
                       
                    break;
                case 3:                                                                 //此时填充虚拟网格为3
                    for (int i = 0; i <= 2; i++){
                        x[i][Ghost_cell-1] = x[i][cols-Ghost_cell-1];             
                        x[i][Ghost_cell-2] = x[i][cols-Ghost_cell-2];             
                        x[i][Ghost_cell-3] = x[i][cols-Ghost_cell-3];             
                    }
                       
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;

        
        case 1:
            switch (Ghost_cell){
                case 1:                                                                 //此时填充虚拟网格为1
                    for (int i = 0; i <= 2; i++)
                        x[i][cols-Ghost_cell] = x[i][Ghost_cell];       
                    break;
                case 2:                                                                 //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[i][cols-Ghost_cell] = x[i][Ghost_cell];                  
                        x[i][cols-Ghost_cell+1] = x[i][Ghost_cell+1];                  
                    }
                    break;
                case 3:                                                                 //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[i][cols-Ghost_cell] = x[i][Ghost_cell];                  
                        x[i][cols-Ghost_cell+1] = x[i][Ghost_cell+1];                  
                        x[i][cols-Ghost_cell+2] = x[i][Ghost_cell+2];                  
                    }
                    break;
                default:
                    break;
            }
            break;
    }
}



//指定固定数值给边界条件
static inline void BC_In_rho(int rows, int cols,double (*x)[cols], int Control_Pos, int Ghost_cell){
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
                            x[i][cols-Ghost_cell] = -x[i][cols-Ghost_cell-1];
                        else
                            x[i][cols-Ghost_cell] = x[i][cols-Ghost_cell-1];
                    }
                    break;

                case 2:                                      //此时填充虚拟网格为2
                    for (int i = 0; i <= 2; i++){
                        x[0][cols-Ghost_cell] = 2.079156;
                        x[0][cols-Ghost_cell+1] = 2.079156;

                        x[1][cols-Ghost_cell] = 0.0;
                        x[1][cols-Ghost_cell+1] = 0.0;

                        x[2][cols-Ghost_cell] =  7.316625;
                        x[2][cols-Ghost_cell+1] =  7.316625;
                    }
                    break;
                //填充更多的虚拟网格，程序待完成
                default:
                    break;
            }
            break;
    }
}


#endif // BOUNDARY_CONDITION_H