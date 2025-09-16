#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "initialize.h"
#include "CFD_convection.h"
#include "CFD_diffusion.h"
#include "Boundary_Condition.h"
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
void Init_DRare(double pri_Ver1[3], double pri_Ver2[3],double u_r);
void Init_Shock_Impact(double pri_Ver1[3], double pri_Ver2[3],double u_r);
void Init_Shock1(double pri_Ver1[3], double pri_Ver2[3]);
void Init_Shock2(double pri_Ver1[3], double pri_Ver2[3]);
void Init_Shock3(double pri_Ver1[3], double pri_Ver2[3]);
void Init_Euler(int rows, int cols ,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double pri_Ver1[3], double pri_Ver2[3],double gamma);

// 读取黎曼问题
int Get_RP_input(int argc, char *argv[]);

//简单的网格代码
void Mesh(int n, double deltax, double *x);



/*                                            *********                                          */
/*                                              主程序                                            */
/*                                            *********                                          */

//网格参数
 

const int n = 200;                                                      //网格数量 
const int RP_Method = 0;                                                //Riemann Solver的具体方法  
const int Recon_Accur = 1;                                              //重构方式控制变量{0是0阶，1是2阶段TVD格式；3是3阶weno，5是5阶weno重构}
                                                                        //TVD包括Vanleer Limter，Minmod limter等，具体在CFD_convection.h中修改：                                              
const int Control_Compution = 0;
const int Control_output = 100;

//Riemann Solver
int scheme;                                                   
int main(int argc, char *argv[]) {

    scheme = Get_RP_input(argc, argv);

    double t = 0;
    double L = 1.,       Tmax = 0.14;       //计算域参数
    double CFL = 0.6; 
    double gamma = 1.4;                     //物性参数
    double Delta_x;
    double Delta_T;

    // 添加声明，2个虚拟网格
    double x[n];
    double pri[3][n+4];
    double U[3][n+4];
    double FU[3][n+4];
    double pri_Ver1[3], pri_Ver2[3];

    Delta_x = L / n; 
    int Ite = 0;

//初始条件
    double u_r = 0.0;

//    Init_Sod(pri_Ver1,pri_Ver2);                        //激波对撞
    Init_Shock_Impact(pri_Ver1,pri_Ver2,u_r); 
//    Init_DRare(pri_Ver1,pri_Ver2,u_r);  
//    Init_Shock3(pri_Ver1,pri_Ver2);  
    printf("Read initial conditions successfully!\n");
  
    //mesh
    Mesh(n,Delta_x,x);
    printf("Mesh successfully!\n");

    Init_Euler(3,n+4,pri,U,FU,pri_Ver1,pri_Ver2,gamma);            //初始化欧拉方程
    printf("Euler equation initialization successful!\n");


    switch (Control_Compution){
        case 0:
            //时间推进：时间一阶和时间二阶格式
            for (t = 0; t < Tmax; t = t+Delta_T) {
                Delta_T = Get_Delta_T(3,n+4,U,Delta_x,CFL,gamma);
                if (t + Delta_T > Tmax)
                    Delta_T = Tmax - t ;

                switch (RP_Method) {
                    case 0:
                        RK1_TVD(Recon_Accur,scheme,3,n+4,x,U,FU,Delta_T,Delta_x,gamma,1,u_r);
                        break;
                    case 1:
                        RK1_TVD_RP_HeatConduction(Recon_Accur,scheme,3,n+4,x,U,FU,Delta_T,Delta_x,gamma,1,u_r);
                        break;
                    default:
                        RK1_TVD_FluxRela(Recon_Accur,scheme,3,n+4,x,U,FU,Delta_T,Delta_x,gamma,1,u_r);     
                        break;
                }
                Ite++;
                if (Ite % Control_output == 0 )
                    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t, Ite);
            }
            break;
        case 1:
            //迭代步数控制
            t = 0;
            for (int m = 0; m <= 100; m++) {
                Delta_T = 0.2*Delta_x;
                switch (RP_Method) {
                    case 0:
                        RK1_TVD(Recon_Accur,scheme,3,n+4,x,U,FU,Delta_T,Delta_x,gamma,1,u_r);
                        break;
                    case 1:
                        RK1_TVD_RP_HeatConduction(Recon_Accur,scheme,3,n+4,x,U,FU,Delta_T,Delta_x,gamma,1,u_r);
                        break;
                    default:
                        RK1_TVD_FluxRela(Recon_Accur,scheme,3,n+4,x,U,FU,Delta_T,Delta_x,gamma,1,u_r);     
                        break;
                }
                Ite++;
                t = t+Delta_T;
                if (Ite % Control_output == 0 ){
                    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t, Ite);
                }
            }
            break;
        default:
            break;
    }
    printf("Step = %d     Time = %f  \nCalculation of step %d is completed \n", Ite, t, Ite);
    printf("end of calculation!\n");

    //实现输出最后的结果
    Con_to_Pri_1D(3,n+4,pri,U,gamma);
    //打开文件并输出结果
    OutputData_file(0,n+4,x,U,FU,pri,t);
    
    printf("The program has completed its execution.\n");
    return 0;

}


// 函数定义：从命令行参数或用户输入获取黎曼求解器方案
int Get_RP_input(int argc, char *argv[]) {
    int selected_scheme;
    
    if (argc > 1) {
        selected_scheme = atoi(argv[1]);
        printf("Using scheme %d from command line argument\n", selected_scheme);
    } else {
        // 显示可用的黎曼求解器选项
        printf("Available Riemann Solvers:\n");
        printf("0: Lax\n");
        printf("1: Rusanov\n");
        printf("2: HLL\n");
        printf("3: HLLC\n");
        printf("4: Roe\n");
        printf("5: Marquina\n");
        printf("6: StegerWarming\n");
        printf("7: VanLeer\n");
        printf("8: LiouSteffen\n");
        printf("Other: Exact Riemann\n");
        
        printf("Please enter the scheme value (0-8): ");
        scanf("%d", &selected_scheme);
        printf("Using scheme %d from user input\n", selected_scheme);
    }

    // 验证输入的有效性
    if (selected_scheme < 0 || selected_scheme >= 10) {
        printf("Warning: Scheme value %d is outside recommended range (0-8)\n", selected_scheme);
        printf("Using Exact Riemann solver as default\n");
    }
    return selected_scheme;
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
    double rho1 = 0.2655737117053071;
    double rho2 =  0.125;
    double u1 = 0.9274526200489499;
    double u2 =  0.0;
    double p1 =  0.30313017805064685;
    double p2 = 0.1;
    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


void Init_DRare(double pri_Ver1[3], double pri_Ver2[3],double u_r) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = -1.0 + u_r;
    double u2 = 1.0 + u_r;
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

void Init_Shock_Impact(double pri_Ver1[3], double pri_Ver2[3],double u_r) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = 4.0 + u_r;
    double u2 = -4.0 + u_r;
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


void Init_Shock1(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 24.0/11;
    double rho2 = 1.0;
    double u1 = 13.0/12;
    double u2 = 0;
    double p1 = 19.0/6;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}


void Init_Shock2(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 600.0/107.0;
    double rho2 = 1.0;
    double u1 = 493.0/60.0;
    double u2 = 0.0;
    double p1 = 499.0/6.0;
    double p2 = 1.0;

    // 将变量赋值给数组
    pri_Ver1[0] = rho1;
    pri_Ver1[1] = u1;
    pri_Ver1[2] = p1;

    pri_Ver2[0] = rho2;
    pri_Ver2[1] = u2;
    pri_Ver2[2] = p2;
}

void Init_Shock3(double pri_Ver1[3], double pri_Ver2[3]) {
    double rho1 = 1.0;
    double rho2 = 1.0;
    double u1 = 4.0;
    double u2 = 4.0;
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


//划分网格代码
void Mesh(int n, double deltax, double *x) {
    // 初始化网格节点坐标
    for (int i = 0; i < n; i++) {
        x[i] = deltax * i + 0.5 * deltax;
    }
}



// 控制计算步数子程序
void StepLoop(double* deltat, double deltax, double CFL, double t, int l, int cols, int vis, int scheme, 
                                            double* x, double (*U)[cols], double (*FU)[cols],double (*pri)[cols],double gamma) {
    for (int k = 0; k < l; k++) {
        //计算格式
        //RK1_TVD(vis, scheme, 3, cols, x, U, FU, 1e-10, deltax, CFL, gamma, k+1);
        //RK_2Roe(vis, scheme, 3, n+4, U, FU, deltat, deltax, CFL);

        t += *deltat;
        Con_to_Pri_1D(3, cols, pri, U,gamma);
        OutputData_file(k+1,cols,x,U,FU,pri,t);
    }
}


void OutputData_file(int k, int cols, double* x, double (*U)[cols], double (*FU)[cols],double (*pri)[cols] , double t) {
    // 文件输出
    printf("Output result(rho, u, p, T, U_M, U_E)\n");
    if (k == 0)
    {
        FILE* file = fopen("/mnt/d/Desktop/RP_FVM/data/output_data.dat","w");
        if (file == NULL) {
            
            printf("File opening failed\n");
            printf("---------------Error----------------\n");
            return; // 返回错误代码
        }
        fprintf(file, "variables=x \t rho \t u\t p\t T\t rhou\t rhoE\t rhoE_K\t rhoE_I\t F_rho\t F_rhou\t F_rhoE\n");
        for (int i = 2; i < cols-2 ; i++) {
            fprintf(file,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",\
                        x[i-2], pri[0][i], pri[1][i], pri[2][i], pri[2][i]/pri[0][i],\
                        U[1][i], U[2][i], 0.5*pri[1][i]*U[1][i], pri[2][i]/(1.4-1) , \
                        U[1][i], U[1][i] * pri[1][i] + pri[2][i] ,pri[1][i] * (U[2][i] + pri[2][i]));
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
