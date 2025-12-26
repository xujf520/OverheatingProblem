#ifndef FUNCTION_H
#define FUNCTION_H

#include <math.h>
#include <omp.h>

#ifndef PI
    #define PI 3.14159265358979323846
#endif

// Function Pointer Type Used To Pass The Integrand
//Modify The Function Pointer Type And Add A U Parameter
typedef double (*FuncPtrWithU)(double x, double* u);

// Example Test Function These Can Also Be Set As Static Inline If They Need To Be Used In Multiple Files
static inline double test_function1(double x) {
    return x * x;  // f(x) = x²，积分结果应为 b³/3 - a³/3
}

static inline double test_function2(double x) {
    return sin(x);  // f(x) = sin(x)，积分结果应为 -cos(b) + cos(a)
}

static inline double test_function3(double x) {
    return exp(x);  // f(x) = e^x，积分结果应为 e^b - e^a
}


static inline int sgn(double num) {
    if (num > 0) {
        return 1;
    } else if (num < 0) {
        return -1;
    } else {
        return 0;
    }
}


static inline double max_of_two(double a, double b) {
    if (a > b)
        return a;
    else 
        return b;
}



static inline double max_of_three(double a, double b, double c) {
    a = max_of_two(a,b);
    return max_of_two(a,c);
}


static inline double max_of_four(double a, double b, double c, double d) {
    a= max_of_three(a,b,c);
    return max_of_two(a,d);
}

static inline double min_of_two(double a, double b) {
    if (a > b)
        return b;
    else 
        return a;
}



static inline double min_mod(double a, double b){
	return 0.5 * (sgn(a) + sgn(b)) * min_of_two(fabs(a), fabs(b));
}

static inline double van_leer(double a, double b){
    double epsilo = 1e-10;
	return ((sgn(a)+sgn(b)) * a * b)/(fabs(a)+fabs(b) + epsilo);
}


static inline double van_albada(double a, double b){
    double epsilo = 1e-6;
	return (fmax(a*b,0) * (a+b))/(pow(a,2)+pow(b,2)+epsilo);
}


static inline double Get_Delta_T(int rows, int cols,double (*x)[cols], double dx) {
    double S_plus = 0;
    for (int j = 1; j < cols-1; j++) {
        //Read The Known Left And Right Original Variables
        double rho = x[0][j];
        double rhou = x[1][j];
        double rhoe = x[2][j];
        double u = rhou / rho;
        double p = (rhoe - 0.5 * rho * pow(u, 2))*(M_gamma-1);
        //Calculate The Speed Of Sound
        double a= sqrt(M_gamma * p / rho);
        //Calculate The Global Maximum Wave Speed
        S_plus = max_of_two (fabs(u) + a,S_plus);
    
    }   

//    return CFL * dx / S_plus;
    return CFL * dx / S_plus;
//    return  0.5 * pow(dx, 5.0/3.0);
}

static inline double Get_Delta_T_2D(int rows, int cols, int depth, double (*x)[cols][depth], double dx, double dy) {
    double S_plus_x = 0.0, S_plus_y = 0.0;

    #pragma omp parallel for reduction(max: S_plus_x, S_plus_y) collapse(2)
    for (int j = GhostCell; j < cols - GhostCell; j++) {
        for (int k = GhostCell; k < depth - GhostCell; k++) {
            // Read The Known Left And Right Original Variables
            double rho = x[0][j][k];
            double rhou = x[1][j][k];
            double rhov = x[2][j][k];
            double rhoe = x[3][j][k];
            double u = rhou / rho;
            double v = rhov / rho;
            double p = (rhoe - 0.5 * rho * (u*u + v*v)) * (M_gamma - 1);
            // Calculate The Speed Of Sound
            double a = sqrt(M_gamma * p / rho);
            //Calculate Local Maximum Wave Speed
            double local_S_plus_x = fabs(u) + a;
            double local_S_plus_y = fabs(v) + a;
            
            // Reduction Operation Updates The Global Maximum
            if (local_S_plus_x > S_plus_x) S_plus_x = local_S_plus_x;
            if (local_S_plus_y > S_plus_y) S_plus_y = local_S_plus_y;
        }
    }
    
    return CFL * min_of_two(dx, dy) / (S_plus_x + S_plus_y);
}

//计算总守恒量
static inline void Total_Conser(int rows, int cols, double Ghost_Cell , double (*x)[cols], double Conser[rows],double delta_x) {
    
    for (int i = 0; i < rows; i++){
        for (int j = Ghost_Cell ; j < cols - Ghost_Cell; j++){
            Conser[i] += x[i][j] * delta_x;  
        } 
        
    }
}

//数值积分程序

// 梯形法则数值积分（修改版本）
// f: 被积函数; x_up: 积分上限; x_down: 积分下限; k: 积分精度; u: 额外参数
static inline double trapezoid_rule(FuncPtrWithU f, double x_down, double x_up, int k, double* u) {
    if (k <= 0) {
        printf("Error: n must be positive\n");
        return 0.0;
    }
    
    double h = (x_up - x_down) / k;                             // 步长
    double sum = 0.5 * (f(x_down, u) + f(x_up, u));            // 端点值，传递 u
    
    for (int i = 1; i < k; i++) {
        double x = x_down + i * h;
        sum += f(x, u);  // 传递 u
    }
    
    return sum * h;
}

// 辛普森法则数值积分（修改版本）
// f: 被积函数; x_up: 积分上限; x_down: 积分下限; k: 积分精度; u: 额外参数
static inline double simpson_rule(FuncPtrWithU f, double x_down, double x_up, int k, double* u) {
    if (k <= 0 || k % 2 != 0) {
        printf("Error: n must be positive and even\n");
        return 0.0;
    }
    
    double h = (x_up - x_down) / k;                             // 步长
    double sum = f(x_down, u) + f(x_up, u);                    // 端点值，传递 u
    
    // 奇数项系数为4，偶数项系数为2
    for (int i = 1; i < k; i++) {
        double x = x_down + i * h;
        if (i % 2 == 0) {
            sum += 2.0 * f(x, u);  // 传递 u
        } else {
            sum += 4.0 * f(x, u);  // 传递 u
        }
    }
    
    return sum * h / 3.0;
}

// 矩形法则数值积分（中点法则，修改版本）
// f: 被积函数; x_up: 积分上限; x_down: 积分下限; k: 积分精度; u: 额外参数
static inline double rectangle_rule(FuncPtrWithU f, double x_down, double x_up, int k, double* u) {
    if (k <= 0) {
        printf("Error: n must be positive\n");
        return 0.0;
    }
    
    double h = (x_up - x_down) / k;                             // 步长
    double sum = 0.0;
    
    for (int i = 0; i < k; i++) {
        double x_mid = x_down + (i + 0.5) * h;                  // 中点
        sum += f(x_mid, u);  // 传递 u
    }
    
    return sum * h;
}

// 数值积分主函数（修改版本）
/* method: 积分方法
 *         1 - 梯形法则
 *         2 - 辛普森法则  
 *         3 - 矩形法则
 */

 
static inline double numerical_integration(FuncPtrWithU f, double a, double b, int n, int method, double* u) {
    switch (method) {
        case 1:
            return trapezoid_rule(f, a, b, n, u);  // 传递 u
        case 2:
            return simpson_rule(f, a, b, n, u);    // 传递 u
        case 3:
            return rectangle_rule(f, a, b, n, u);  // 传递 u
        default:
            printf("Error: Unknown method. Using trapezoid rule.\n");
            return trapezoid_rule(f, a, b, n, u);  // 传递 u
    }
}

// 计算五阶格式的数值误差（保持不变）
static inline double Smooth_function(double x, double* u) {
    // 积分平均值
    double v1 = *(u-2);
    double v2 = *(u-1);
    double v3 = *(u);
    double v4 = *(u+1);
    double v5 = *(u+2);

    // 计算多项式系数：
    double a4 = (v1+v5 - 4.0*(v2+v4) + 6.0*v3)/24.0;
    double a3 = (v5-v1 - 2.0*(v4-v2))/12.0;
    double a2 = (v4+v2)/2.0 - v3 -3.0*a4/2.0;
    double a1 = (v4-v2)/2.0 -5.0*a3/4.0;
    double a0 = v3 - a2/12.0 - a4/80.0;

//    return sin(4.0*PI*x) + 2 - a4*pow(x,4) - a3*pow(x,3) - a2*pow(x,2) - a1*x - a0;  // 简化了 pow(x,1) 和 pow(x,0)

    return sin(4.0*PI*x) + 2 -v3;  // 简化了 pow(x,1) 和 pow(x,0)
}



#endif