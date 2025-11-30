#ifndef CFD_CONVECTION_H
#define CFD_CONVECTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Golbal.h"
#include "Reconstruction.h"
#include "scheme.h"
#include "function.h"

/*                               ************************************                               */
/*                               ************************************                               */
/*                                      重构格式：TVD类格式                                          */
/*                               ************************************                               */
/*                               ************************************                               */


/*                                      ******************                                          */
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */

                                    /*……………………………………………………*/
                                        /*近似黎曼求解*/
                                    /*……………………………………………………*/


static inline void Flux_Reconstruction_RP(int AR_scheme, int var, int rows, int cols, int GC, double (*y)[rows][cols], \
                                                double (*f)[rows][cols],double (*g)[rows][cols], double dt, double dx, double dy) {
    int i,j,k;

    double (*Flux)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conserl)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    double (*Conserr)[rows][cols] = malloc(var * sizeof(double[rows][cols]));
    // 检查内存分配是否成功
    if (Flux == NULL || Conserl == NULL || Conserr == NULL) {
        fprintf(stderr, "Memory allocation failed in Flux_Reconstruction_RP\n");
        // 释放已分配的内存
        free(Flux);
        free(Conserl);
        free(Conserr);
        return;
    }

    //维度分裂
    //x方向重构
    int space_dir = 1;
    switch (Recon_Accur){
        case 1:
            Reconstruction_Godunov(space_dir,var,rows,cols,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
            TVD_Reconstruction(space_dir,var,rows,cols,GC,y,Conserl,Conserr,dx,dy);
            break;
        case 3:
            WENO3_Reconstruction(space_dir,var,rows,cols,GC,y,Conserl,Conserr);
            break;
        case 5:
            WENO5_Reconstruction(space_dir,var,rows,cols,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }

    
    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 1:
            HLL_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 2:
            HLLC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 3:
            Roe_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 11:
            HLLHC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 22:
            HLLCHC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 33:
            RoeHC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 5:
//            RS_Marquina(3,rows,y,Flux,dt,dx);
            break;
        default:
//            ER_Flux(var, rows, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }


    for ( i = 0; i < var; i++)
        for ( j = GC-1; j <= rows-GC; j++)
            for ( k = GC-1; k <= cols-GC; k++)
                f[i][j][k] = Flux[i][j][k];

    
    //y方向重构

    space_dir = 2;
    switch (Recon_Accur){
        case 1:
            Reconstruction_Godunov(space_dir,var,rows,cols,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
            TVD_Reconstruction(space_dir,var,rows,cols,GC,y,Conserl,Conserr,dx,dy);
            break;
        case 3:
            WENO3_Reconstruction(space_dir,var,rows,cols,GC,y,Conserl,Conserr);
            break;
        case 5:
            WENO5_Reconstruction(space_dir,var,rows,cols,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }


    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 1:
            HLL_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 2:
            HLLC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 3:
            Roe_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 11:
            HLLHC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 22:
            HLLCHC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 33:
            RoeHC_Flux(space_dir,var, rows,cols, GC, Conserl,Conserr,Flux);
            break;
        case 5:
//            RS_Marquina(3,rows,y,Flux,dt,dx);
            break;
        default:
//            ER_Flux(var, rows, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }



    for ( i = 0; i < var; i++)
        for ( j = GC-1; j <= rows-GC; j++)
            for ( k = GC-1; k <= cols-GC; k++)
                g[i][j][k] = Flux[i][j][k];


    free(Conserl);
    free(Conserr);
    free(Flux);

}


/*                                      ******************                                          */
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */

                                    /*……………………………………………………*/
                                /*具有人工热传导的近似黎曼求解数值方法*/
                                    /*……………………………………………………*/
/*……………………………………………………………………………………………………*/


static inline void Flux_Reconstruction_RP_Heat(int AR_scheme, int var, int rows, int GC, double (*y)[rows],double (*z)[rows]\
                                ,double dt,double dx, double u_Refer) {
    int i,j;
    double Conserl[var][rows],Conserr[var][rows];
    double prileft[var][rows], priright[var][rows];
    double Flux[var][rows];
    
     switch (Recon_Accur){
        case 1:
//            Reconstruction_Godunov(var,rows,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
//            TVD_Reconstruction(var,rows,GC,y,Conserl,Conserr,dx);
            break;
        case 3:
//            WENO3_Reconstruction(var,rows,GC,y,Conserl,Conserr);
            break;
        case 5:
//            WENO5_Reconstruction(var,rows,GC,y,Conserl,Conserr);
//            WENO5_Reconstruction_C(var,rows,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }
    
//    Con_to_Pri_1D(3,rows,prileft,Conserl);
//   Con_to_Pri_1D(3,rows,priright,Conserr);


    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 2:
            HLL_Flux_HeatConduction(var, rows, GC, prileft,priright,Flux,u_Refer);
            break;
        case 3:
            HLLC_Flux_HeatConduction(var, rows, GC, prileft,priright,Flux,u_Refer);
            break;
        case 4:
            Roe_Flux_HeatConduction(var, rows, GC, prileft,priright,Flux,u_Refer);
            break;
        case 5:
            RS_Marquina(3,rows,y,z,dt,dx);
            break;
        default:
            ER_Flux_PlusHeat(var, rows, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }


    for ( i = 0; i < var; i++){
        for ( j = 1; j < rows-1; j++){
            z[i][j] = Flux[i][j];
        }
    }

}



/*……………………………………………………………………………………………………*/

/*                                      ******************                                          */
/*                                      重构步：MUSCL类格式                                           */
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */
                                    /*……………………………………………………*/
                                        /*近似黎曼求解*/
                                    /*……………………………………………………*/




#endif 
