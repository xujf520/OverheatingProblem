#ifndef CFD_CONVECTION_H
#define CFD_CONVECTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
/*……………………………………………………………………………………………………*/
//双激波近似黎曼求解
static inline void RS_HLL(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma, double u_r) {
    int i,j;
    double coverl[3][cols],coverr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double Flux[3][cols];
    double slope[3][cols-1],a[3][cols-1],b[3][cols-1];


    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            a[i][j]=(y[i][j] - y[i][j-1])/dx;
        }
    }

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            b[i][j]=(y[i][j+1] - y[i][j])/dx;
        }
    }
    
    //重构
    for (int i = 0; i < rows; i++){   
        for (int j = 1; j < cols-1; j++)
        {
            switch (Recon_Accur)
            {
                case 0:
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }
    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverl[i][j]=y[i][j]+0.5*slope[i][j]*dx;
        }
    }
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverr[i][j]=y[i][j+1]-0.5*slope[i][j+1]*dx;
        }
    }
    //边界条件，降价为0阶重构
    for ( i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }

    

    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);
    
    HLL_Flux(rows, cols, prileft,priright,Flux,gamma);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            z[i][j] = Flux[i][j];
        }
    }

}


/*……………………………………………………………………………………………………*/
//三波近似黎曼求解
static inline void RS_HLLC(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma) {
    int i,j;
    double coverl[3][cols],coverr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double Flux[3][cols];
    double slope[3][cols-1],a[3][cols-1],b[3][cols-1];

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            a[i][j]=(y[i][j] - y[i][j-1])/dx;
        }
    }

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            b[i][j]=(y[i][j+1] - y[i][j])/dx;
        }
    }
    
    //重构
    for (int i = 0; i < rows; i++){   
        for (int j = 1; j < cols-1; j++)
        {
            switch (Recon_Accur)
            {
                case 0:
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }
    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverl[i][j]=y[i][j]+0.5*slope[i][j]*dx;
        }
    }
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverr[i][j]=y[i][j+1]-0.5*slope[i][j+1]*dx;
        }
    }
    //边界条件，降价为0阶重构
    for ( i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }
    
    
    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);

    HLLC_Flux(rows, cols, prileft,priright,Flux,gamma);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            z[i][j] = Flux[i][j];
        }
    }

}

/*……………………………………………………………………………………………………*/
//Roe近似黎曼求解并且存在两种重构&&作用于守恒变量
static inline void RS_Roe(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma) {
    //直接作用守恒变量
    int i,j;
    double coverl[rows][cols],coverr[rows][cols];
    double slope[rows][cols-1],a[rows][cols-1],b[rows][cols-1];
    double prileft[rows][cols], priright[rows][cols];
    double Flux[rows][cols];
    
    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            a[i][j]=(y[i][j] - y[i][j-1])/dx;
        }
    }

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            b[i][j]=(y[i][j+1] - y[i][j])/dx;
        }
    }
    
    //重构

    //重构
    for (int i = 0; i < rows; i++){   
        for (int j = 1; j < cols-1; j++)
        {
            switch (Recon_Accur)
            {
                case 0:
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }

    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverl[i][j]=y[i][j]+0.5*slope[i][j]*dx;
        }
    }
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverr[i][j]=y[i][j+1]-0.5*slope[i][j+1]*dx;
        }
    }
    //边界条件，降价为0阶重构
    for ( i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }
    
    //计算Roe平均
    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);

    Roe_Flux(rows, cols, prileft,priright,Flux,gamma);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            z[i][j] = Flux[i][j];
        }
    }
}




/*……………………………………………………………………………………………………*/
//精确黎曼求解
static inline void RS_ER(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma) {
    int i,j;
    double coverl[3][cols],coverr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double fluxl[3][cols],fluxr[3][cols];
    double flux[3][cols];
    double slope[3][cols-1],a[3][cols-1],b[3][cols-1];

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            a[i][j]=(y[i][j] - y[i][j-1])/dx;
        }
    }

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            b[i][j]=(y[i][j+1] - y[i][j])/dx;
        }
    }
    
    //重构
    //重构
    for (int i = 0; i < rows; i++){   
        for (int j = 1; j < cols-1; j++)
        {
            switch (Recon_Accur)
            {
                case 0:
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }
    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverl[i][j] = y[i][j];
        }
    }
    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverl[i][j]=y[i][j]+ 0.5  * slope[i][j]*dx;
        }
    }
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverr[i][j]=y[i][j+1] - 0.5 * slope[i][j+1]*dx;
        }
    }
    //边界条件，降价为0阶重构
    for ( i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }
    
    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);
    initEulerflux1D(rows, cols, coverl, fluxl,gamma);
    initEulerflux1D(rows, cols, coverr, fluxr,gamma);
    
    ER_Flux(rows, cols, prileft,priright,flux,gamma);
    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            z[i][j] = flux[i][j];
        }
    }


}


/*……………………………………………………………………………………………………*/
//“The algorithmic description of Marquina’s flux formula is as follows:” ([Donat 和 Marquina, 1996, p. 44]
static inline void RS_Marquina(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma) {
    double coverl[3][cols],coverr[3][cols];
    double fluxl[3][cols],fluxr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double pri[rows][cols];
    double slope[3][cols-1],a[3][cols-1],b[3][cols-1];
    
    double eigen_l[3][3][cols],eigen_r[3][3][cols];
    double w_l[rows][cols], w_r[rows][cols];
    double phi_fl[rows][cols], phi_fr[rows][cols];
    double phi_fp[rows][cols], phi_fm[rows][cols];
    double flux[3][cols];
    double lamda[3][cols];
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


    Con_to_Pri_1D(3,cols,pri,y,gamma);
    //计算特征矩阵
	for(int j=0;  j < cols; j++) {
			
        double  q2, c2, b1, b2;
        double _u, _H, _c;
        //preparing some interval value
        _u = pri[1][j];
        _H = 0.5 * pow(pri[1][j],2) + gamma * pri[2][j] /((gamma-1) * pri[0][j]);
        q2 = _u*_u ;
        c2 = (gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
        _c = sqrt(c2);
        b1 = (gamma - 1.0)/c2;
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


        //计算对应特征值
        lamda[0][j] = _u - _c ;
        lamda[1][j] = _u;
        lamda[2][j] = _u + _c;
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
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-2; j++)
            coverl[i][j]=y[i][j] + 0.5  * slope[i][j]*dx;

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-2; j++)
            coverr[i][j]=y[i][j+1] - 0.5 * slope[i][j+1]*dx;
    
            

    //边界条件，降价为0阶重构
    for (int i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }
    
    //执行计算Marquina flux
    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);

    initEulerflux1D(rows, cols, coverl, fluxl,gamma);
    initEulerflux1D(rows, cols, coverr, fluxr,gamma);

    // 投影到特征空间
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_l[k][j] += coverl[m][j]*eigen_l[k][m][j];
                phi_fl[k][j] += fluxl[m][j]*eigen_l[k][m][j];
            }
        }
        
    }
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_r[k][j] += coverr[m][j]*eigen_l[k][m][j+1];
                phi_fr[k][j] += fluxr[m][j]*eigen_l[k][m][j+1];
            }
        }
        
    }

    //计算通量分量
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            if (lamda[k][j]  * lamda[k][j+1] > 0){
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
                alpha[k][j] = max_of_two(fabs(lamda[k][j]),fabs(lamda[k][j+1]));
                //alpha[k][j] = dx/dt;
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
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */

                                    /*……………………………………………………*/
                                /*具有相对速度的近似黎曼求解数值方法*/
                                    /*……………………………………………………*/
/*……………………………………………………………………………………………………*/
//双激波近似黎曼求解
static inline void RS_HLL_XRela(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma, double u_r) {
    int i,j;
    double coverl[3][cols],coverr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double Flux[3][cols];
    double slope[3][cols-1],a[3][cols-1],b[3][cols-1];


    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            a[i][j]=(y[i][j] - y[i][j-1])/dx;
        }
    }

    for ( i = 0; i < rows; i++)
    {
        for ( j = 1; j < cols-1; j++){
            b[i][j]=(y[i][j+1] - y[i][j])/dx;
        }
    }
    
    //重构
    for (int i = 0; i < rows; i++){   
        for (int j = 1; j < cols-1; j++)
        {
            switch (Recon_Accur)
            {
                case 0:
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }
    
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverl[i][j]=y[i][j]+0.5*slope[i][j]*dx;
        }
    }
    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-2; j++){
            coverr[i][j]=y[i][j+1]-0.5*slope[i][j+1]*dx;
        }
    }
    //边界条件，降价为0阶重构
    for ( i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }

    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);
    
    HLL_Flux_XRela(rows, cols, prileft,priright,Flux,gamma,u_r);

    for ( i = 0; i < rows; i++){
        for ( j = 1; j < cols-1; j++){
            z[i][j] = Flux[i][j];
        }
    }


}



/*……………………………………………………………………………………………………*/
//“The algorithmic description of Marquina’s flux formula is as follows:” ([Donat 和 Marquina, 1996, p. 44]
static inline void RS_Marquina_XRela(int Recon_Accur,int rows, int cols,double (*y)[cols],double (*z)[cols]\
                                ,double dt,double dx, double cfl,double gamma,double u_r) {
    double coverl[3][cols],coverr[3][cols];
    double fluxl[3][cols],fluxr[3][cols];
    double prileft[3][cols], priright[3][cols];
    double pri[rows][cols];
    double slope[3][cols-1],a[3][cols-1],b[3][cols-1];
    
    double eigen_l[3][3][cols],eigen_r[3][3][cols];
    double w_l[rows][cols], w_r[rows][cols];
    double phi_fl[rows][cols], phi_fr[rows][cols];
    double phi_fp[rows][cols], phi_fm[rows][cols];
    double flux[3][cols];
    double lamda[3][cols];
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


    Con_to_Pri_1D(3,cols,pri,y,gamma);
    //计算特征矩阵
	for(int j=0;  j < cols; j++) {
			
        double  q2, c2, b1, b2;
        double _u, _H, _c;
        //preparing some interval value
        _u = pri[1][j];
        _H = 0.5 * pow(pri[1][j],2) + gamma * pri[2][j] /((gamma-1) * pri[0][j]);
        q2 = _u*_u ;
        c2 = (gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
        _c = sqrt(c2);
        b1 = (gamma - 1.0)/c2;
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
        _H = 0.5 * pow(pri[1][j] - u_r,2) + gamma * pri[2][j] /((gamma-1) * pri[0][j]);
        q2 = _u*_u ;
        c2 = (gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
        _c = sqrt(c2);

        //计算对应特征值
        lamda[0][j] = _u - _c ;
        lamda[1][j] = _u;
        lamda[2][j] = _u + _c;
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
                    slope[i][j]=0;
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
                    printf(" The reconstruction program with a precision greater than the  %d-th order has not been implemented.\n",Recon_Accur);
                    exit(1);
            }
        }
    }

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-2; j++)
            coverl[i][j]=y[i][j] + 0.5  * slope[i][j]*dx;

    for (int i = 0; i < rows; i++)
        for (int j = 1; j < cols-2; j++)
            coverr[i][j]=y[i][j+1] - 0.5 * slope[i][j+1]*dx;
    
            

    //边界条件，降价为0阶重构
    for (int i = 0; i < rows; i++){
        coverr[i][0] = y[i][1];
        coverl[i][0] = y[i][0];
        coverr[i][cols-2] = y[i][cols-1];
        coverl[i][cols-2] = y[i][cols-2];
    }
    
    //执行计算Marquina flux
    Con_to_Pri_1D(3,cols,prileft,coverl,gamma);
    Con_to_Pri_1D(3,cols,priright,coverr,gamma);

    initEulerflux1D(rows, cols, coverl, fluxl,gamma);
    initEulerflux1D(rows, cols, coverr, fluxr,gamma);


    // 投影到特征空间
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_l[k][j] += coverl[m][j]*eigen_l[k][m][j];
                phi_fl[k][j] += fluxl[m][j]*eigen_l[k][m][j];
            }
        }
        
    }
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            for (int m = 0; m < 3; m++){
                w_r[k][j] += coverr[m][j]*eigen_l[k][m][j+1];
                phi_fr[k][j] += fluxr[m][j]*eigen_l[k][m][j+1];
            }
        }
        
    }

    //计算通量分量
    for (int k = 0; k < 3; k++){
        for (int j = 1; j < cols-1; j++){
            if (lamda[k][j] * lamda[k][j+1] >= 0){
                if (lamda[k][j] >= 0){
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
/*                                      重构步：基于原始变量                                          */
/*                                      ******************                                          */




/*                                      ******************                                          */
/*                                      重构步：基于特征变量                                          */
/*                                      ******************                                          */

/*……………………………………………………………………………………………………*/
//Roe近似黎曼求解并且存在两种重构（TVD）&&作用特征变量


/*                                      ******************                                          */
/*                                      重构步：MUSCL类格式                                           */
/*                                      重构步：基于守恒变量                                          */
/*                                      ******************                                          */
                                    /*……………………………………………………*/
                                        /*近似黎曼求解*/
                                    /*……………………………………………………*/






#endif 
