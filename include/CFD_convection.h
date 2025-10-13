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


static inline void Flux_Reconstruction_RP(int AR_scheme, int rows, int cols, int depth, int GC, \
                                            double (*y)[cols][depth],double (*f)[cols][depth],double (*g)[cols][depth], double dt, double dx) {
    int i,j,k;

    double (*Flux)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
    double (*Conserl)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
    double (*Conserr)[cols][depth] = malloc(rows * sizeof(double[cols][depth]));
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
            Reconstruction_Godunov(space_dir,rows,cols,depth,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
            TVD_Reconstruction(space_dir,rows,cols,depth,GC,y,Conserl,Conserr,dx);
            break;
        case 3:
            WENO3_Reconstruction(space_dir,rows,cols,depth,GC,y,Conserl,Conserr);
            break;
        case 5:
            WENO5_Reconstruction(space_dir,rows,cols,depth,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }

    
    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 2:
            HLL_Flux(space_dir,rows, cols,depth, GC, Conserl,Conserr,Flux);
            break;
        case 3:
//            HLLC_Flux(rows, cols, GC, prileft,priright,Flux);
            break;
        case 4:
//            Roe_Flux(rows, cols, GC, prileft,priright,Flux);
            break;
        case 5:
//            RS_Marquina(3,cols,y,Flux,dt,dx);
            break;
        default:
//            ER_Flux(rows, cols, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }


    for ( i = 0; i < rows; i++)
        for ( j = GC-1; j <= cols-GC; j++)
            for ( k = GC-1; k <= cols-GC; k++)
                f[i][j][k] = Flux[i][j][k];

    
    //y方向重构

    space_dir = 2;
    switch (Recon_Accur){
        case 1:
            Reconstruction_Godunov(space_dir,rows,cols,depth,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
            TVD_Reconstruction(space_dir,rows,cols,depth,GC,y,Conserl,Conserr,dx);
            break;
        case 3:
            WENO3_Reconstruction(space_dir,rows,cols,depth,GC,y,Conserl,Conserr);
            break;
        case 5:
            WENO5_Reconstruction(space_dir,rows,cols,depth,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }


    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 2:
            HLL_Flux(space_dir,rows, cols,depth, GC, Conserl,Conserr,Flux);
            break;
        case 3:
//            HLLC_Flux(rows, cols, GC, prileft,priright,Flux);
            break;
        case 4:
//            Roe_Flux(rows, cols, GC, prileft,priright,Flux);
            break;
        case 5:
//            RS_Marquina(3,cols,y,Flux,dt,dx);
            break;
        default:
//            ER_Flux(rows, cols, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }


    for ( i = 0; i < rows; i++)
        for ( j = GC-1; j <= cols-GC; j++)
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


static inline void Flux_Reconstruction_RP_Heat(int AR_scheme, int rows, int cols, int GC, double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double u_Refer) {
    int i,j;
    double Conserl[rows][cols],Conserr[rows][cols];
    double prileft[rows][cols], priright[rows][cols];
    double Flux[rows][cols];
    
     switch (Recon_Accur){
        case 1:
//            Reconstruction_Godunov(rows,cols,GC,y,Conserl,Conserr,dx);
            break;
        case 2:
//            TVD_Reconstruction(rows,cols,GC,y,Conserl,Conserr,dx);
            break;
        case 3:
//            WENO3_Reconstruction(rows,cols,GC,y,Conserl,Conserr);
            break;
        case 5:
//            WENO5_Reconstruction(rows,cols,GC,y,Conserl,Conserr);
//            WENO5_Reconstruction_C(rows,cols,GC,y,Conserl,Conserr);
            break;

        default:
            printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
            exit(1);
    }
    
    Con_to_Pri_1D(3,cols,prileft,Conserl);
    Con_to_Pri_1D(3,cols,priright,Conserr);


    //演化过程：
    //AR_scheme is Approximate Riemann Solver
    switch (AR_scheme) {
        case 2:
            HLL_Flux_HeatConduction(rows, cols, GC, prileft,priright,Flux,u_Refer);
            break;
        case 3:
            HLLC_Flux_HeatConduction(rows, cols, GC, prileft,priright,Flux,u_Refer);
            break;
        case 4:
            Roe_Flux_HeatConduction(rows, cols, GC, prileft,priright,Flux,u_Refer);
            break;
        case 5:
            RS_Marquina(3,cols,y,z,dt,dx);
            break;
        default:
            ER_Flux_PlusHeat(rows, cols, GC, prileft,priright,Flux);
            // 你可以根据实际需求添加相应的处理逻辑
            break;
    }


    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            z[i][j] = Flux[i][j];
        }
    }

}



/*……………………………………………………………………………………………………*/
//“The algorithmic description of Marquina’s flux formula is as follows:” ([Donat 和 Marquina, 1996, p. 44]
static inline void RS_Marquina_XRela(int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx,double u_r) {
    double Conserl[rows][cols],Conserr[rows][cols];
    double fluxl[rows][cols],fluxr[rows][cols];
    double prileft[rows][cols], priright[rows][cols];
    double pri[rows][cols];
    double slope[rows][cols-1],a[rows][cols-1],b[rows][cols-1];
    
    double eigen_l[rows][rows][cols],eigen_r[rows][rows][cols];
    double w_l[rows][cols], w_r[rows][cols];
    double phi_fl[rows][cols], phi_fr[rows][cols];
    double phi_fp[rows][cols], phi_fm[rows][cols];
    double flux[rows][cols];
    double lamda[rows][cols], lamda_plus[rows][cols];
    double alpha[rows][cols];

    //初始化数组
     for (int k = 0; k < rows; k++){
        for (int j = 0; j < cols; j++){
            w_l[k][j] = 0.0;
            w_r[k][j] = 0.0;
            phi_fl[k][j] = 0.0;
            phi_fr[k][j] = 0.0;
            phi_fp[k][j] = 0.0;
            phi_fm[k][j] = 0.0;
            flux[k][j] = 0.0;
        }
        
    }


    Con_to_Pri_1D(3,cols,pri,y);
    //计算特征矩阵
	for(int j=0;  j < cols; j++) {
			
        double  q2, c2, b1, b2;
        double _u, _H, _c;
        //preparing some interval value
        _u = pri[1][j];
        _H = 0.5 * pow(pri[1][j],2) + M_gamma * pri[2][j] /((M_gamma-1) * pri[0][j]);
        q2 = _u*_u ;
        c2 = (M_gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
        _c = sqrt(c2);
        b1 = (M_gamma - 1.0)/c2;
        b2 = 1.0 + b1*q2 - b1*_H;
        // left eigen vectors 
        eigen_l[0][0][j] = 0.5*(b2 + _u/_c);
        eigen_l[0][1][j] = -0.5*(b1*_u + 1/_c);
        eigen_l[0][2][j] = 0.5*b1;
            
        eigen_l[1][0][j] = -q2 + _H;
        eigen_l[1][1][j] = _u;;
        eigen_l[1][2][j] = -1.0;

        eigen_l[2][0][j] = 0.5*(b2 - _u/_c);
        eigen_l[2][1][j] = 0.5*(-b1*_u + 1/_c);
        eigen_l[2][2][j] = 0.5*b1;

        //right eigen vectors
        eigen_r[0][0][j] = 1.0;
        eigen_r[0][1][j] = b1;
        eigen_r[0][2][j] = 1.0;
            
        eigen_r[1][0][j] = _u - _c;
        eigen_r[1][1][j] = _u*b1;
        eigen_r[1][2][j] = _u + _c;

        eigen_r[2][0][j] = _H - _u*_c;
        eigen_r[2][1][j] = _H*b1 - 1.0;
        eigen_r[2][2][j] = _H + _u*_c;
    }


    for(int j=0;  j < cols; j++) {
			
        double  q2, c2;
        double _u, _H, _c;
        //preparing some interval value
        _u = pri[1][j] - u_r;
        _H = 0.5 * pow(pri[1][j] - u_r,2) + M_gamma * pri[2][j] /((M_gamma-1) * pri[0][j]);
        q2 = _u*_u ;
        c2 = (M_gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
        _c = sqrt(c2);

        //计算对应特征值
        lamda[0][j] = _u - _c;
        lamda[1][j] = _u;
        lamda[2][j] = _u + _c;

        lamda_plus[0][j] = fabs(_u) + _c ;
        lamda_plus[1][j] = _u;
        lamda_plus[2][j] = fabs(_u) + _c ;
    }
    

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-1; j++)
            a[i][j]=(y[i][j] - y[i][j-1])/dx;

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-1; j++)
            b[i][j]=(y[i][j+1] - y[i][j])/dx;
    
    //重构
    //重构
    for (int i = 0; i < rows; i++){   
        for (int j = 1; j < cols-1; j++)
        {
            switch (Recon_Accur)
            {
                case 0:
//                    Reconstruction_Godunov(rows,cols,y,Conserl,Conserr,dx);
                    break;
                case 1:
                    //vanleer
                    slope[i][j]= ((sgn(a[i][j])+sgn(b[i][j]))*a[i][j]*b[i][j])/(fabs(a[i][j])+fabs(b[i][j])+1e-15);
                    //minibee
                    //slope[i][j] = 0.5 * (sgn(a[i][j])+sgn(b[i][j])) * min_of_two(fabs(a[i][j]), fabs(b[i][j]));
                    //vanalbada
                    //slope[i][j]=(fmax(a[i][j]*b[i][j],0) * (a[i][j]+b[i][j]))/(pow(a[i][j],2)+pow(b[i][j],2)+10e-6);
                    break;
                
                default:
                    printf("The reconstruction program with the %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-2; j++)
            Conserl[i][j]=y[i][j] + 0.5  * slope[i][j]*dx;

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-2; j++)
            Conserr[i][j]=y[i][j+1] - 0.5 * slope[i][j+1]*dx;
    
            

    //边界条件，降价为0阶重构
    for (int i = 0; i < rows; i++){
        Conserr[i][0] = y[i][1];
        Conserl[i][0] = y[i][0];
        Conserr[i][cols-2] = y[i][cols-1];
        Conserl[i][cols-2] = y[i][cols-2];
    }
    
    //执行计算Marquina flux
    Con_to_Pri_1D(3,cols,prileft,Conserl);
    Con_to_Pri_1D(3,cols,priright,Conserr);

    initEulerflux1D(rows, cols, Conserl, fluxl);
    initEulerflux1D(rows, cols, Conserr, fluxr);


    // 投影到特征空间
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_l[k][j] += Conserl[m][j]*eigen_l[k][m][j];
                phi_fl[k][j] += fluxl[m][j]*eigen_l[k][m][j];
            }
        }
        
    }
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_r[k][j] += Conserr[m][j]*eigen_l[k][m][j+1];
                phi_fr[k][j] += fluxr[m][j]*eigen_l[k][m][j+1];
            }
        }
        
    }

    //计算通量分量
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            if (lamda[k][j] * lamda[k][j+1] > 0.0){
                if (lamda[k][j] > 0){
                    phi_fp[k][j] = phi_fl[k][j];
                    phi_fm[k][j] = 0.0;
                }
                else{
                    phi_fp[k][j] = 0.0;
                    phi_fm[k][j] = phi_fr[k][j];
                }
            }
            else{
                alpha[k][j] = max_of_two(fabs(lamda[k][j])+fabs(u_r),fabs(lamda[k][j+1])+fabs(u_r));
                phi_fp[k][j] = 0.5*(phi_fl[k][j] + alpha[k][j] * w_l[k][j]);
                phi_fm[k][j] = 0.5*(phi_fr[k][j] - alpha[k][j] * w_r[k][j]);
            }
        }
    }

    //投影到物理空间计算Flux
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
               flux[k][j] += phi_fp[m][j]*eigen_r[k][m][j] + phi_fm[m][j] * eigen_r[k][m][j+1];
            }
        }  
    }

    //边界条件，降价为0阶重构
    for (int i = 0; i < rows; i++){
        flux[i][0] = flux[i][2];
        flux[i][1] = flux[i][2];
        flux[i][cols-1] = flux[i][cols-3];
        flux[i][cols-2] = flux[i][cols-3];
    }


    for (int i = 0; i < rows; i++){
        for (int  j = 1; j < cols-1; j++){
            z[i][j] = flux[i][j];
        }
    }

}


/*                                      ******************                                          */
/*                                      重构步：MUSCL类格式                                           */
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */
                                    /*……………………………………………………*/
                                        /*近似黎曼求解*/
                                    /*……………………………………………………*/




#endif 
