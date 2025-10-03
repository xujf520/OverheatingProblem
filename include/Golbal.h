// global.h
#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdbool.h>

//部分逻辑变量
static bool Characteriz = false;
//重构方式控制变量{0是0阶，2是2阶段TVD格式；3是3阶weno，5是5阶weno重构}
//TVD包括Vanleer Limter，Minmod limter等，具体在CFD_convection.h中修改： 
static int Recon_Accur = 5;



//虚拟网格参数
static int GhostCell = 3;



//物质属性参数：
static double M_gamma = 1.4;



#endif