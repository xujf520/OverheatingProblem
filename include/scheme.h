#ifndef WENO_FLUX_H
#define WENO_FLUX_H

#include <stdio.h>
#include <math.h>

//提前声明部分函数
double ExactRieamnna_pstar(double rhol, double rhor, double ul,double ur, double pl, \
							double pr, double gammal,double gammar ,double tol, double maxit);
double ExactRieamnna_ustar(double p_star, double rhol, double rhor, double ul,double ur, double pl,double pr, double gammal,double gammar);

//将中间变量第二次还原
static inline void restore2(double*u1,double*u2,double*u3,int n,double*roarr,double*uarr,double*parr){
    int i;
    for(i = 0;i < n;i++){
        roarr[i] = u1[i];
    }
    for(i = 0;i < n;i++){
        uarr[i] = u2[i]/u1[i];
    }
    for(i = 0;i < n;i++){
        parr[i] = (u3[i]-(0.5*pow(u2[i],2))/u1[i])*0.4;
    }
}


static inline void Con_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols],double gamma){
    int j;
    for (j = 0; j < cols-1; j++) {
        x[0][j] = y[0][j];
        x[1][j] = y[1][j]/y[0][j];
        x[2][j] = (gamma-1)*(y[2][j] - 0.5*y[1][j]*y[1][j]/y[0][j]);
    }
}

static inline void Pri_to_Con_1D(int rows, int cols, double (*x)[cols], double (*y)[cols],double gamma){
    int j;
    for (j = 0; j < cols-1; j++) {
        y[0][j] = x[0][j];
        y[1][j] = x[0][j] * x[1][j];
        y[2][j] = 0.5 * x[0][j] * pow(x[1][j],2) + x[2][j]/(gamma-1);
    }
}

static inline void Pri_to_S_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j < cols-1; j++) {
        y[0][j] = x[2][j] / pow(x[0][j], 1.4);
        y[1][j] = x[1][j];
        y[2][j] = x[2][j];
    }
}

static inline void S_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j < cols-1; j++) {
        x[0][j] = pow(y[2][j]/y[0][j],1/1.4);
        x[1][j] = y[1][j];
        x[2][j] = y[2][j];
    }
}

static inline void ConS_to_Pri_1D(int rows, int cols, double (*x)[cols], double (*y)[cols]){
    int j;
    for (j = 0; j < cols-1; j++) {
        x[0][j] = y[0][j];
        x[1][j] = y[1][j] / y[0][j];
        x[2][j] = y[2][j] * pow(y[0][j],0.4);
    }
}


//Flux计算方法
static inline void HLL_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {


        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数
        //计算声速

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;

        //确定HLL数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;
        if (sleft >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && sright >0){
            rho_F = (sright*rho_FL - sleft*rho_FR + sleft*sright * (rho_R - rho_L))/(sright-sleft);
            rhou_F = (sright*rhou_FL - sleft*rhou_FR + sleft*sright * (rhou_R - rhou_L))/(sright-sleft);
            rhoe_F = (sright*rhoe_FL - sleft*rhoe_FR + sleft*sright * (rhoe_R - rhoe_L))/(sright-sleft);
        }
        else if (sright <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void HLLC_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_L = rho_L * u_L;
        double rhou_R = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_L = 0.5 * rho_L * pow(u_L, 2) + p_L/(gamma-1);
        double rhoe_R = 0.5 * rho_R * pow(u_R, 2) + p_R/(gamma-1);
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));

        //利用Roe平均的变量计算近似波速
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));
        double sleft = ubar - cbar;
        double sright = ubar + cbar;
        double s_star = (p_R - p_L + rho_L*u_L*(sleft - u_L) - rho_R*u_R*(sright - u_R))\
                                 /(rho_L*(sleft - u_L) - rho_R*(sright - u_R));

        double u_stat_L = rho_L * (sleft-u_L)/(sleft-s_star);
        double u_stat_R = rho_R * (sright-u_R)/(sright-s_star);
        double fe_star_L = rhoe_L/rho_L + (s_star-u_L)*(s_star + p_L/(rho_L *(sleft-u_L)));
        double fe_star_R = rhoe_R/rho_R + (s_star-u_R)*(s_star + p_R/(rho_R *(sright-u_R)));


        //HLLC数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        if (sleft >= 0 ){
            rho_F = rho_FL;
            rhou_F = rhou_FL;
            rhoe_F = rhoe_FL;
        }
        else if(sleft < 0 && s_star >=0 ){
            rho_F = rho_FL + sleft * (u_stat_L - rho_L);
            rhou_F = rhou_FL + sleft *(u_stat_L * s_star  - rhou_L);
            rhoe_F = rhoe_FL + sleft * (u_stat_L * fe_star_L -  rhoe_L);

        }
        else if(s_star < 0 && sright > 0){
            rho_F = rho_FR + sright * (u_stat_R - rho_R);
            rhou_F = rhou_FR + sright *(u_stat_R * s_star  - rhou_R);
            rhoe_F = rhoe_FR + sright * (u_stat_R * fe_star_R -  rhoe_R);
        }
        else if (sright <= 0){
            rho_F = rho_FR;
            rhou_F = rhou_FR;
            rhoe_F = rhoe_FR;
        }

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void Roe_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    double epsilon = 1e-10;
    
    for (int j = 1; j < cols-1; j++) {
        //读取已知的左右原始变量
        double rho_L = x[0][j];
        double rho_R = y[0][j];
        double u_L = x[1][j];
        double u_R = y[1][j];
        double p_L = x[2][j];
        double p_R = y[2][j];

        //计算需要使用的参数

        //计算总焓H
        double H_L = 0.5 * pow(u_L,2) + (gamma/(gamma-1) ) * (p_L/rho_L);
        double H_R = 0.5 * pow(u_R,2) + (gamma/(gamma-1) ) * (p_R/rho_R);

        //计算左右守恒变量和通量
        double rho_FL = rho_L * u_L;
        double rho_FR = rho_R * u_R;
        double rhou_FL = rho_L * u_L*u_L + p_L;
        double rhou_FR = rho_R * u_R*u_R + p_R;
        double rhoe_FL = rho_L * H_L * u_L;
        double rhoe_FR = rho_R * H_R * u_R;
    
        //计算Roe平均物理量
        double rhobar = sqrt(rho_L * rho_R); 
        double ubar = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R)/(sqrt(rho_L)+sqrt(rho_R));
        double Hbar = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R)/(sqrt(rho_L)+sqrt(rho_R));
        double cbar = sqrt((gamma-1) * (Hbar - 0.5 * pow(ubar,2)));


        //利用Roe平均计算稳定项

        double lambda1, lambda2, lambda3;
        double alpha1, alpha2, alpha3;
        lambda1 = fabs(ubar - cbar);
        if (lambda1 < epsilon){
            lambda1 = (fabs(ubar - cbar) + pow(epsilon,2)) / (2 * epsilon);
        }
        lambda2 = fabs(ubar);
        if (lambda2 < epsilon){
            lambda2 = (fabs(ubar) + pow(epsilon,2)) / (2 * epsilon);
        }
        lambda3 = fabs(ubar + cbar);
        if (lambda3 < epsilon){
            lambda3 = (fabs(ubar + cbar) + pow(epsilon,2)) / (2 * epsilon);
        }

        alpha1 = ((p_R - p_L) - rhobar * cbar * (u_R - u_L)) / (2 * pow(cbar,2));
        alpha2 = (rho_R - rho_L) - (p_R - p_L) / pow(cbar,2);
        alpha3 = ((p_R - p_L) + rhobar * cbar * (u_R - u_L)) / (2 * pow(cbar,2));


        //Roe数值通量
        double rho_F = 0, rhou_F = 0, rhoe_F = 0;

        rho_F = 0.5*(rho_FL + rho_FR) - 0.5 *(lambda1*alpha1 + lambda2*alpha2 + lambda3*alpha3);
        rhou_F = 0.5 * (rhou_FL + rhou_FR)\
                    - 0.5*(lambda1 * alpha1 * (ubar- cbar)\
                    + lambda2 * alpha2 * ubar\
                    + lambda3 * alpha3 * (ubar + cbar));
        rhoe_F = 0.5*(rhoe_FL + rhoe_FR)\
                    -0.5*(lambda1 * alpha1 * (Hbar - ubar*cbar)\
                    +lambda2 * alpha2 * 0.5 * pow(ubar,2)\
                    +lambda3 * alpha3 * (Hbar + ubar*cbar));

        z[0][j] = rho_F; 
        z[1][j] = rhou_F; 
        z[2][j] = rhoe_F; 
    }   
}



static inline void ER_Flux(int rows, int cols,double (*x)[cols], double (*y)[cols] ,double (*z)[cols],double gamma) {
    
    for (int j = 1; j < cols-1; j++) {
        double rho_F = 0,u_F=0,p_F=0;
        double rho_l = x[0][j];
        double rho_r = y[0][j];
        double u_l = x[1][j];
        double u_r = y[1][j];
        double p_l = x[2][j];
        double p_r = y[2][j];
        double c_l = sqrt(gamma*p_l/rho_l);
        double c_r = sqrt(gamma*p_r/rho_r);
        double p_star = ExactRieamnna_pstar(rho_l, rho_r, u_l, u_r, p_l, p_r, gamma, gamma ,1e-10, 1e6);
		double u_star = ExactRieamnna_ustar(p_star, rho_l, rho_r, u_l,u_r, p_l,p_r, gamma,gamma);
        
         
        //分别考虑左激波，左稀疏波的情况
        if (p_star >=  p_l){
            double part1 = (gamma + 1)* p_star + (gamma - 1) * p_l ;
            double part2 = (gamma - 1) * p_star + (gamma + 1) * p_l;
            double rho_star_L = rho_l * part1 / part2;
            double A_L = 2.0 / ((gamma + 1) * rho_l);
            //计算激波的速度
            double B_L = (gamma - 1) * p_l  / (gamma + 1);
            double S_L = u_l - (1 / rho_l) * sqrt((B_L + p_star) / A_L); 
            double si = 0;
            if (si <= S_L){
                rho_F = rho_l;
                u_F = u_l;
                p_F = p_l;
            }
            else if (S_L < si && si <= u_star)
            {
                rho_F = rho_star_L;
                u_F = u_star;
                p_F = p_star;
            }
        }
        else {
            //计算接触间断左侧密度
            double part1 = p_star;
            double part2 = p_l;
            double rho_star_L = rho_l * pow(part1/part2, 1/gamma);
            //计算左稀疏波波头和波尾的速度
            double aL_star = c_l * pow((p_star / p_l), ((gamma - 1) / (2 *gamma)));
            double S_HL = u_l - c_l;
            double S_TL = u_star - aL_star;
            double si=0;
            if (si <= S_HL){
                rho_F = rho_l;
                u_F = u_l;
                p_F = p_l;
            }
            else if (si <= S_TL && si > S_HL)
            {
                rho_F = rho_l * pow(2/(gamma + 1) + (gamma - 1) * (u_l - si) / ((gamma + 1) * c_l), (2 / (gamma - 1))); 
                u_F = 2 * (c_l + (gamma - 1) * u_l/2 + si) / (gamma + 1);
                p_F = p_l * pow(2 / (gamma + 1) + (gamma - 1) * (u_l - si) / ((gamma + 1) * c_l),  2 * gamma / (gamma - 1));
            }
            else if (si <= u_star && si > S_TL)
            {
                rho_F = rho_star_L;
                u_F = u_star;
                p_F = p_star;
            }
        }

        //分别考虑右稀疏波和右激波的情况
        if ( p_star >=  p_r){
            double part1 =  (gamma + 1) * p_star + (gamma - 1) * p_r;
            double part2 =  (gamma - 1) * p_star + (gamma + 1) * p_r;
            double rho_star_r = rho_r * part1 / part2;
            double A_r = 2 / ((gamma + 1) * rho_r);
            //计算激波的速度
            double B_r = (gamma - 1) * p_r / (gamma + 1);
            double S_r = u_r + (1 / rho_r) * sqrt((B_r + p_star) / A_r); 
            double si = 0;
            if (si <= S_r &&  u_star < si)
            {
                rho_F = rho_star_r;
                u_F = u_star;
                p_F = p_star;
            }
            else if (S_r < 0){
                rho_F = rho_r;
                u_F = u_r;
                p_F = p_r;
            }
        }
        else{
            //计算接触间断左侧密度
            double part1 = p_star;
            double part2 = p_r;
            double rho_star_r = rho_r * pow(part1 / part2, 1/gamma);
            //计算左稀疏波波头和波尾的速度
            double aR_star = c_r * pow((p_star / p_r), ((gamma - 1) / (2 *gamma)));
            double S_HR = u_r + c_r;
            double S_TR = u_star + aR_star;   
            double si=0;
            if (si <= S_TR && si > u_star)
            {
                rho_F = rho_star_r;
                u_F = u_star;
                p_F = p_star;
            }    
            else if (si<= S_HR && si >S_TR)
            {
                rho_F = rho_r * pow(2/(gamma + 1) - (gamma - 1) * (u_r - si) / ((gamma + 1) * c_r), (2 / (gamma - 1))); 
                u_F = 2 * (-c_r + (gamma - 1) * u_r / 2 + si) / (gamma + 1);
                p_F = p_r * pow(2 / (gamma + 1) - (gamma - 1) * (u_r - si) / ((gamma + 1) * c_r),  2 * gamma / (gamma - 1));
            }
            else if (si > S_HR){
                rho_F = rho_r;
                u_F = u_r;
                p_F = p_r;
            }
        }

        z[0][j] = rho_F * u_F; 
        z[1][j] = rho_F * u_F *u_F + p_F; 
        z[2][j] = ((0.5*rho_F * u_F *u_F + p_F/(gamma-1)) + p_F) * u_F; 
    }   
}



//这里是计算精确Riemann解

double F_K_P(double p, double p_refer, double rho_refer, double gamma,double c_refer){
	//计算f_K(p) K=L/R
    //:param p: 变量 求解过程中的压强
    //:param _refer: 固定量要么是L 要么是R
    //:return: f_K 函数值
	if (p >=  p_refer)
	{
		double A_K = 2 / ((gamma + 1) * rho_refer);
        double B_K = (gamma - 1) * (p_refer) / (gamma + 1);
        return (p - p_refer) * sqrt(A_K / (B_K + p ));
	}
	else
		return (2 * c_refer / (gamma - 1)) * (pow(p/p_refer,(gamma - 1) / (2 *gamma)) - 1);
}

double dF_K_P(double p, double p_refer, double rho_refer, double gamma,double c_refer){
    //计算f'_K(p) K=L/R 它是 f_K(p) 关于 p 的导数
    //:param p: 变量 求解过程中的压强
    //:param _refer: 固定量要么是L 要么是R
	//:return:f'_K(p) 函数值
	if (p >=  p_refer)
	{
		double A_K = 2 / ((gamma + 1) * rho_refer);
        double B_K = (gamma - 1) * (p_refer) / (gamma + 1);
		double part1 = sqrt((A_K) / (B_K + p ));
        double part2 = 1 - (p - p_refer) / (2 * (B_K + p));
        return part1 * part2;

	}
	else{
		double part1 = pow(p/p_refer, -(gamma + 1) / (2 * gamma));
        double part2 = 1 / (rho_refer * c_refer);
        return part1 * part2;
	}
}


// 定义一个函数来计算f(x)
double ExactRieamnna_pstar(double rhol, double rhor, double ul,double ur, double pl, \
							double pr, double gammal,double gammar ,double tol, double maxit){
	
	//利用状态方方程计算声速度
	double cl = sqrt(gammal * pl/rhol);
	double cr = sqrt(gammar * pr/rhor);
	//计算迭代初始值
    double p = max_of_two(tol, 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rhor + rhol) * (cl+cr));
    double p_new = 0.0;
	double f_p = 0; 
	double df_p = 0;

	for ( int i = 0; i < maxit; i++)
	{
		//迭代初始值问题
		if (p < 0.0)
		{
			printf("negative pressure\n");
			exit(1); 
		}

		//计算Fpp的函数
        f_p = F_K_P(p, pl,rhol,gammal,cl) + F_K_P(p, pr,rhor,gammar,cr) + ur - ul;
        df_p = dF_K_P(p, pl,rhol,gammal,cl) + dF_K_P(p, pr,rhor,gammar,cr);
        p_new = p - f_p / df_p;

		if(p_new < 0.0 ){ 
            printf("negative pressure\n");
			exit(1); 
		}
		if (2 * fabs(p_new - p) / (p + p_new) < tol){
			p = p_new;
            break;
		}
		else{
			p = p_new;
		}
	}


	double p_star = p;
	return p_star;
}

double ExactRieamnna_ustar(double p_star, double rhol, double rhor, double ul,double ur, double pl, \
							double pr, double gammal,double gammar){
	
	//利用状态方方程计算声速度
	double cl = sqrt(gammal * pl/rhol);
	double cr = sqrt(gammar * pr/rhor);
	//计算迭代初始值
	double u_star = 0.5 * (ul + ur) + 0.5 * (F_K_P(p_star, pr,rhor,gammar,cr) - F_K_P(p_star, pl,rhol,gammal,cl));
    return u_star;
}


#endif  