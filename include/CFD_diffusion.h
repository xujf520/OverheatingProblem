#ifndef CFD_DIFFUSION_H
#define CFD_DIFFUSION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static inline void CFD_diffusion1up(float*alongx, float*grad_diff, float n, float u, float dx) {
    int i;
    for(i=0;i<n;i++){                    //迎风差分
        if(i==0){
            grad_diff[i]=(alongx[i]-u)/(2*pow(dx,2));
        }
        if(i==1){
            grad_diff[i]=(alongx[i]-2*alongx[i-1]+u)/(2*pow(dx,2));
        }
        if(i > 1){
            grad_diff[i]=(alongx[i]-2*alongx[i-1]+alongx[i-2])/(2*pow(dx,2));
        }
    }
}



static inline void CFD_diffusion1CD(float*alongx, float*grad_diff, float n, float u, float dx){
    int i = 0;
    for ( i = 0; i < n; i++){                    //一阶中心差分
        if (i == 0){
            grad_diff[i] = (alongx[i+1] - 2*alongx[i] + u) / pow(dx, 2);
        }
        if (i > 0){
            grad_diff[i] = (alongx[i+1] - 2*alongx[i] + alongx[i-1]) / pow(dx, 2);
        }
    }
}

#endif 