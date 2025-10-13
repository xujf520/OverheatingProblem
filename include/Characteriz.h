#ifndef CHARACTERIZ_H
#define CHARACTERIZ_H

#include <stdio.h>
#include <math.h>
#include "Golbal.h"


static inline void Compute_Eigen(int rowss, int rows, int cols, double (*Pri)[cols], double (*eigen_l)[rows][cols], double (*eigen_r)[rows][cols]) {
    
    for(int j = 0; j < cols; j++) {
        double q2, c2, b1, b2;
        double _u, _H, _c, _uc;
        
        // 计算中间变量
        _u = Pri[1][j];
        _H = 0.5 * pow(Pri[1][j], 2) + M_gamma * Pri[2][j] / ((M_gamma - 1) * Pri[0][j]);
        q2 = _u * _u;
        c2 = (M_gamma - 1.0) * (_H - 0.5 * q2);  // 从焓计算声速
        _c = sqrt(c2);
        _uc = _u * _c;
        b1 = (M_gamma - 1.0) / (2.0 * c2);
        b2 = 1.0 + b1 * q2 - b1 * _H;
        
        // 左特征向量
        eigen_l[0][0][j] = 0.5 * (b2 + _u / _c);
        eigen_l[0][1][j] = -0.5 * (b1 * _u + 1 / _c);
        eigen_l[0][2][j] = 0.5 * b1;
            
        eigen_l[1][0][j] = -q2 + _H;
        eigen_l[1][1][j] = _u;
        eigen_l[1][2][j] = -1.0;

        eigen_l[2][0][j] = 0.5 * (b2 - _u / _c);
        eigen_l[2][1][j] = 0.5 * (-b1 * _u + 1 / _c);
        eigen_l[2][2][j] = 0.5 * b1;

        // 右特征向量
        eigen_r[0][0][j] = 1.0;
        eigen_r[0][1][j] = b1;
        eigen_r[0][2][j] = 1.0;
            
        eigen_r[1][0][j] = _u - _c;
        eigen_r[1][1][j] = _u * b1;
        eigen_r[1][2][j] = _u + _c;

        eigen_r[2][0][j] = _H - _u * _c;
        eigen_r[2][1][j] = _H * b1 - 1.0;
        eigen_r[2][2][j] = _H + _u * _c;
    }
}



static inline void Compute_Eigen_2D(double AA, double BB, int rows, int cols, int depth,\
                                    double (*Pri)[cols][depth], double (*eigen_l)[rows][cols][depth],\
                                    double (*eigen_r)[rows][cols][depth]){

    
    double (*H)[depth] = malloc(cols * sizeof(double[depth]));
    
    // 检查内存分配是否成功
    if (H == NULL) {
        fprintf(stderr, "Memory allocation failed in Cumpute_Eigen_2D\n");
        // 释放已分配的内存
        free(H);
        return;
    }

    for (int j = 1; j < cols-1; j++){ 
        for (int k = 1; k < depth-1; k++){
            H[j][k] = 0.5 * (pow(Pri[1][j][k],2) + pow(Pri[2][j][k],2)) + (M_gamma/(M_gamma-1))*(Pri[3][j][k]/Pri[0][j][k]);
        }
    }


	for(int j = GhostCell-1; j <= cols-GhostCell; j++) {
		for(int k = GhostCell-1; k <= depth-GhostCell; k++) {
			double D, D1, q2, uu, vv, c2, b1, b2;
			double _u, _v, _H, _c;
            //preparing some interval value
            D = sqrt((AA*Pri[0][j+1][k] + BB*Pri[0][j][k+1])/Pri[0][j][k]);
            D1 = D + 1.0;
            _u = (Pri[1][j][k] + D*(AA*Pri[1][j+1][k] + BB*Pri[1][j][k+1]))/D1;
            _v = (Pri[2][j][k] + D*(AA*Pri[2][j+1][k] + BB*Pri[2][j][k+1]))/D1;
            _H = (H[j][k] + D*(AA*H[j+1][k] + BB*H[j][k+1]))/D1;
            q2 = _u*_u + _v*_v;
            uu = AA*_u + BB*_v;
            vv = AA*_v - BB*_u;
            c2 = (M_gamma - 1.0)*(_H - 0.5*q2); 						//sound speed form H
            _c = sqrt(c2);
            b1 = (M_gamma - 1.0)/c2;
            b2 = 1.0 + b1*q2 - b1*_H;

            // left eigen vectors 
            eigen_l[0][0][j][k] = 0.5*(b2 + uu/_c);
            eigen_l[0][1][j][k] = -0.5*(b1*_u + AA/_c);
            eigen_l[0][2][j][k] = -0.5*(b1*_v + BB/_c);
            eigen_l[0][3][j][k] = 0.5*b1;
            
            eigen_l[1][0][j][k] = -q2 + _H;
            eigen_l[1][1][j][k] = _u;
            eigen_l[1][2][j][k] = _v;
            eigen_l[1][3][j][k] = -1.0;

            eigen_l[2][0][j][k] = vv;
            eigen_l[2][1][j][k] = BB;
            eigen_l[2][2][j][k] = -AA;
            eigen_l[2][3][j][k] = 0.0;

            eigen_l[3][0][j][k] = 0.5*(b2 - uu/_c);
            eigen_l[3][1][j][k] = 0.5*(-b1*_u + AA/_c);
            eigen_l[3][2][j][k] = 0.5*(-b1*_v + BB/_c);
            eigen_l[3][3][j][k] = 0.5*b1;

            //right eigen vectors
            eigen_r[0][0][j][k] = 1.0;
            eigen_r[0][1][j][k] = b1;
            eigen_r[0][2][j][k] = 0.0;
            eigen_r[0][3][j][k] = 1.0;
            
            eigen_r[1][0][j][k] = _u - AA*_c;
            eigen_r[1][1][j][k] = _u*b1;
            eigen_r[1][2][j][k] = BB;
            eigen_r[1][3][j][k] = _u + AA*_c;

            eigen_r[2][0][j][k] = _v - BB*_c;
            eigen_r[2][1][j][k] = _v*b1;
            eigen_r[2][2][j][k] = -AA;
            eigen_r[2][3][j][k] = _v + BB*_c;

            eigen_r[3][0][j][k] = _H - uu*_c;
            eigen_r[3][1][j][k] = _H*b1 - 1.0;
            eigen_r[3][2][j][k] = -vv;
            eigen_r[3][3][j][k] = _H + uu*_c;


		}
	}
    
    free(H);
}


#endif 