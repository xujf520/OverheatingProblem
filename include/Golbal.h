// global.h
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdbool.h>

//计算控制变量
static int Control_Compution = 0;
static int Control_output = 100;


//重构方式控制变量{0是0阶，2是2阶段TVD格式；3是3阶weno，5是5阶weno重构}
//TVD包括Vanleer Limter，Minmod limter等，具体在CFD_convection.h中修改： 
static int Recon_Accur = 2;
//是否特征重构
static bool Characteriz = false;


//边界条件控制：
static bool Periodicity = false;
static bool Reflect = true;

//程序计算参数
static double CFL = 0.4; 

//虚拟网格参数
static int GhostCell = 4;

//物质属性参数：
static double M_gamma = 1.4;

#endif