#ifndef CHARACTERIZ_H
#define CHARACTERIZ_H

#include <stdio.h>
#include <math.h>
#include "Golbal.h"


static inline void Compute_Eigen_Matrix(int rowss, int rows, int cols, double (*Pri)[cols], double (*eigen_l)[rows][cols], double (*eigen_r)[rows][cols]) {
    
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


#endif  