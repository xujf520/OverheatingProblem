#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "initialize.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "time_advance.h"
#include "scheme.h"


/*                                            *********                                          */
/*                                            声明子程序                                          */
/*                                            *********                                          */
// 控制计算步数子程序
void StepLoop();
void OutputData_file();

//初始条件子程序
void Init_Sod(double pri_Ver1[3], double pri_Ver2[3]);
void Init_Sod_Rare(double pri_Ver1[3], double pri_Ver2[3]);
void Init_Sod_Shock(double pri_Ver1[3], double pri_Ver2[3]);
void Init_DRare(double pri_Ver1[3], double pri_Ver2[3]);
void Init_DRare_S(double pri_Ver1[3], double pri_Ver2[3]);
void Init_Shock_Impact(double pri_Ver1[3], double pri_Ver2[3],double u_r);
void Init_Euler(int rows, int cols ,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double pri_Ver1[3], double pri_Ver2[3],double gamma);
void Init_Euler_S(int rows, int cols ,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double pri_Ver1[3], double pri_Ver2[3]);

//简单的网格代码
void Mesh(int n, double deltax, double *x);


/*                                            *********                                          */
/*                                              主程序                                            */
/*                                            *********                                          */
const int n = 200;
const int vis = 1;                      //重构方式控制变量{0是0阶，1是TVD；高精度数值方法则0是3阶weno，1是5阶weno}
const int scheme = 2;                   //近似黎曼问题求解器的控制变量            &&& 控制是否使用相对参考系的Riemann Solver

int main(){
    //定义初始参数
    
    double t = 0;
    double l = 1.,       tmax = 0.15;       //计算域参数
    double cfl = 0.4; 
    double gamma = 1.4;                     //物性参数
    double deltax;
    double * deltat = malloc(sizeof(double));

    // 添加声明，2个虚拟网格
    double x[n];
    double pri[3][n + 4];
    double U[3][n + 4];
    double FU[3][n + 4];
    double U_S[3][n + 4];
    double FU_S[3][n + 4];
    double pri_Ver1[3], pri_Ver2[3];

    deltax = l / n; 
    *deltat= 0.0;
    double  dt = 1e-4;
    int Ite = 0;


//初始条件
    double u_r = 0.0;

    Init_Shock_Impact(pri_Ver1,pri_Ver2,u_r);                        //激波对撞
    printf("Read initial conditions successfully!\n");
  
    //mesh
    Mesh(n,deltax,x);
    printf("Mesh successfully!\n");
   
    Init_Euler(3,n+4,pri,U,FU,pri_Ver1,pri_Ver2,gamma);            //初始化欧拉方程
   
    

    //时间推进：时间一阶和时间二阶格式
    for (t = 0; t < tmax; t = t+dt) {
        RK1_TVD(vis,scheme,3,n+4,x,U,FU,dt,deltax,cfl,gamma,1,u_r);
        //RK1_TVD_FluxM(vis,scheme,3,n+4,x,U,FU,dt,deltax,cfl,gamma,1,u_r);                                   //Riemann问题是否设置相对运动
        Ite++;
        if (Ite % 100 == 0 ){
              printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t, Ite);
        }
        
    }
    printf("end of calculation!\n");

    //实现输出最后的结果
   //ConS_to_Pri_1D(3,n+4,pri,U_S);
    Con_to_Pri_1D(3,n+4,pri,U,gamma);
    //打开文件并输出结果
    OutputData_file(0,n+4,x,U_S,FU_S,pri,t);
    
    free(deltat);

    printf("The program has completed its execution.\n");
    return 0;

}


//初始条件
void Init_Sod(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 1.0;
    double rho2 = 0.125;
    double u1 = -0;
    double u2 = 0;
    double p1 = 1.0;
    double p2 = 0.1;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

void Init_Sod_Rare(double pri_Ver1[3], double pri_Ver2[3]) {


    double rho1 = 1.0;
    double rho2 =  0.43757818061324344;
    //double rho2 =  0.5;
    double u1 = 0.0;
    double u2 =  0.9013775087441291 ;
    double p1 = 1.0;
    double p2 =  0.31439665844271514;
    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

void Init_Sod_Shock(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 0.23750813461769094; 
    double rho2 =  0.125;
    double u1 = 0.9013775087441291;
    double u2 =  0.0;
    double p1 = 0.3143966584427151;
    double p2 = 0.1;
    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


void Init_DRare(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = -2.0;
    double u2 = 2.0;
    double p1 = 0.4;
    double p2 = 0.4;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

void Init_DRare_S(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = -2.5;
    double u2 = 2.5;
    double p1 = 1.0;
    double p2 = 0.5;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

void Init_Shock_Impact(double pri_Ver1[3], double pri_Ver2[3],double u_r) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = 1.0 + u_r;
    double u2 = -1.0 + u_r;
    double p1 = 1.0;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


void Init_Euler(int rows, int cols ,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double pri_Ver1[3], double pri_Ver2[3] ,double gamma){
    //初始化
    initEuler1D(3, cols, x);
    initEuler1D(3, cols, y);
    initEuler1D(3, cols, z);
    //进一步初始化
    initEulerpri1D_Shocktube(3, cols, x,pri_Ver1,pri_Ver2);
    //initEulerpri1D_Osher(3, cols, pri, x);
    initEulerconser1D(3, cols,x, y, gamma);
    initEulerflux1D(3, cols, y,z, gamma);
}

void Init_Euler_S(int rows, int cols ,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double pri_Ver1[3], double pri_Ver2[3]){
    //初始化
    initEuler1D(3, cols, x);
    initEuler1D(3, cols, y);
    initEuler1D(3, cols, z);
    //进一步初始化
    initEulerpri1D_Shocktube(3, cols, x,pri_Ver1,pri_Ver2);
    init_ConS(3, cols,x, y);
    init_FluxS(3, cols, x,z);
}


//划分网格代码
void Mesh(int n, double deltax, double *x) {
    // 初始化网格节点坐标
    for (int i = 0; i < n; i++) {
        x[i] = deltax * i + 0.5 * deltax;
    }
}



// 控制计算步数子程序
void StepLoop(double* deltat, double deltax, double cfl, double t, int l, int cols, int vis, int scheme, 
                                            double* x, double (*U)[cols], double (*FU)[cols],double (*pri)[cols],double gamma) {
    for (int k = 0; k < l; k++) {
        //计算格式
        //RK1_TVD(vis, scheme, 3, cols, x, U, FU, 1e-10, deltax, cfl, gamma, k+1);
        //RK_2Roe(vis, scheme, 3, n+4, U, FU, deltat, deltax, cfl);

        t += *deltat;
        Con_to_Pri_1D(3, cols, pri, U,gamma);
        OutputData_file(k+1,cols,x,U,FU,pri,t);
    }
}


void OutputData_file(int k, int cols, double* x, double (*U)[cols], double (*FU)[cols],double (*pri)[cols] , double t) {
    // 文件输出

    if (k == 0)
    {
        FILE* file = fopen("/mnt/d/Desktop/RP_FVM/data/output_data.dat","w");
        if (file == NULL) {
            
            printf("File opening failed\n");
            printf("---------------Error----------------\n");
            return; // 返回错误代码
        }
        fprintf(file, "variables=x \t rho \t u\t p\t T\n");
        for (int i = 2; i < cols-2 ; i++) {
            fprintf(file,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",x[i-2], pri[0][i], pri[1][i], pri[2][i],pri[2][i]/pri[0][i]);
        }
        fclose(file);
        printf("Output calculation result successful\n");
       // printf("%f\n",t);
    }
    else
    {
        char filename[100];
        sprintf(filename, "/mnt/d/Desktop/RP_FVM/data/output_%d.dat", k);
        FILE* file = fopen(filename, "w");
        if (file == NULL) {
            printf("无法打开文件 %s\n", filename);
            return;
        }
        fprintf(file, "variables =  'x' 'rho' 'u' 'p' 'rhou' \n");
        for (int i = 2; i < cols-2; i++) {
            fprintf(file, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", 
            x[i-2], pri[0][i], pri[1][i], pri[2][i], U[1][i]);
        }
        fclose(file);
        printf("The %d th calculation ended at %f\n", k, t);
    }
}
