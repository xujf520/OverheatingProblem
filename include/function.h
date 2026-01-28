#ifndef FUNCTION_H
#define FUNCTION_H

#include <math.h>
#include <omp.h>

#ifndef PI
    #define PI 3.14159265358979323846
#endif

// Function pointer type for integrand functions with additional parameter u
typedef double (*FuncPtrWithU)(double x, double* u);

// Test function examples (can be made static inline if needed in multiple files)
static inline double test_function1(double x) {
    return x * x;  // f(x) = x², integral result should be b³/3 - a³/3
}

static inline double test_function2(double x) {
    return sin(x);  // f(x) = sin(x), integral result should be -cos(b) + cos(a)
}

static inline double test_function3(double x) {
    return exp(x);  // f(x) = e^x, integral result should be e^b - e^a
}

// Sign function: returns 1 for positive numbers, -1 for negative, 0 for zero
static inline int sgn(double num) {
    if (num > 0) 
        return 1;
    else if (num < 0)
        return -1;
    else 
        return 0;
}

// Returns the maximum of two numbers
static inline double max_of_two(double a, double b) {
    if (a > b)
        return a;
    else 
        return b;
}

// Returns the maximum of three numbers
static inline double max_of_three(double a, double b, double c) {
    a = max_of_two(a,b);
    return max_of_two(a,c);
}

// Returns the maximum of four numbers
static inline double max_of_four(double a, double b, double c, double d) {
    a = max_of_three(a,b,c);
    return max_of_two(a,d);
}

// Returns the minimum of two numbers
static inline double min_of_two(double a, double b) {
    if (a > b)
        return b;
    else 
        return a;
}

// Minmod limiter function commonly used in flux limiters
static inline double min_mod(double a, double b){
    return 0.5 * (sgn(a) + sgn(b)) * min_of_two(fabs(a), fabs(b));
}

// Van Leer limiter function for MUSCL schemes
static inline double van_leer(double a, double b){
    double epsilo = 1.0e-6;
    return ((sgn(a)+sgn(b)) * a * b) / (fabs(a)+fabs(b)+epsilo);
}

// Van Albada limiter function for smooth solutions
static inline double van_albada(double a, double b){
    double epsilo = 1e-6;
    return (fmax(a*b,0) * (a+b))/(pow(a,2)+pow(b,2) + epsilo);
}

// SuperBee limiter function for sharp discontinuity capturing
static inline double SuperBee(double a, double b){
    double epsilo = 1e-6;
    return 0.5 * (sgn(a) + sgn(b)) * max_of_two(
        min_of_two(2. * fabs(a), fabs(b)), 
        min_of_two(fabs(a), 2. * fabs(b))
    );
}

/* UltraBee limiter function (commented out for reference)
static inline double UltraBee(double a, double b){
    double epsilo = 1e-6;
    return (sgn(a) + sgn(b)) * min_of_two(
        fabs(a) / (fabs(c) + epsm), 
        abs(y) / (1 - abs(c) + epsm)
    );
}*/

/**
 * Calculate time step Δt for 1D simulations based on CFL condition
 * Δt = CFL * Δx / max(|u| + a) where a is speed of sound
 * 
 * @param rows Number of rows (unused but kept for interface consistency)
 * @param cols Number of columns in the grid
 * @param x Array of conservative variables [3][cols]
 * @param dx Grid spacing in x-direction
 * @return Time step Δt satisfying CFL condition
 */
static inline double Get_Delta_T(int rows, int cols, double (*x)[cols], double dx) {
    double S_plus = 0;
    for (int j = 1; j < cols-1; j++) {
        // Read conservative variables
        double rho = x[0][j];
        double rhou = x[1][j];
        double rhoe = x[2][j];
        
        // Calculate primitive variables
        double u = rhou / rho;
        double p = (rhoe - 0.5 * rho * pow(u, 2)) * (M_gamma - 1);
        
        // Calculate speed of sound
        double a = sqrt(M_gamma * p / rho);
        
        // Update maximum wave speed
        S_plus = max_of_two(fabs(u) + a, S_plus);
    }   

    return CFL * dx / S_plus;
}

/**
 * Calculate time step Δt for 2D simulations based on CFL condition
 * Δt = CFL * min(Δx, Δy) / (max(|u|+a) + max(|v|+a))
 * 
 * @param rows Number of variables (unused but kept for interface consistency)
 * @param cols Number of grid points in x-direction
 * @param depth Number of grid points in y-direction
 * @param x Array of conservative variables [4][cols][depth]
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @return Time step Δt satisfying CFL condition
 */
static inline double Get_Delta_T_2D(int rows, int cols, int depth, double (*x)[cols][depth], double dx, double dy) {
    double S_plus_x = 0.0, S_plus_y = 0.0;

    #pragma omp parallel for reduction(max: S_plus_x, S_plus_y) collapse(2)
    for (int j = GhostCell; j < cols - GhostCell; j++) {
        for (int k = GhostCell; k < depth - GhostCell; k++) {
            // Read conservative variables
            double rho = x[0][j][k];
            double rhou = x[1][j][k];
            double rhov = x[2][j][k];
            double rhoe = x[3][j][k];
            
            // Calculate primitive variables
            double u = rhou / rho;
            double v = rhov / rho;
            double p = (rhoe - 0.5 * rho * (u*u + v*v)) * (M_gamma - 1);
            
            // Calculate speed of sound
            double a = sqrt(M_gamma * p / rho);
            
            // Calculate local maximum wave speeds
            double local_S_plus_x = fabs(u) + a;
            double local_S_plus_y = fabs(v) + a;
            
            // Update global maximum wave speeds through reduction
            if (local_S_plus_x > S_plus_x) S_plus_x = local_S_plus_x;
            if (local_S_plus_y > S_plus_y) S_plus_y = local_S_plus_y;
        }
    }
    
    return CFL * min_of_two(dx, dy) / (S_plus_x + S_plus_y);
}

/**
 * Alternative method for calculating time step Δt for 2D simulations
 * Uses harmonic mean of Δtx and Δty: Δt = CFL * (Δtx * Δty) / (Δtx + Δty)
 * 
 * @param rows Number of variables (unused but kept for interface consistency)
 * @param cols Number of grid points in x-direction
 * @param depth Number of grid points in y-direction
 * @param x Array of conservative variables [4][cols][depth]
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @return Time step Δt satisfying CFL condition
 */
static inline double Get_Delta_T_2D_X(int rows, int cols, int depth, double (*x)[cols][depth], double dx, double dy) {
    double S_plus_x = 0.0, S_plus_y = 0.0;

    #pragma omp parallel for reduction(max: S_plus_x, S_plus_y) collapse(2)
    for (int j = GhostCell; j < cols - GhostCell; j++) {
        for (int k = GhostCell; k < depth - GhostCell; k++) {
            // Read conservative variables
            double rho = x[0][j][k];
            double rhou = x[1][j][k];
            double rhov = x[2][j][k];
            double rhoe = x[3][j][k];
            
            // Calculate primitive variables
            double u = rhou / rho;
            double v = rhov / rho;
            double p = (rhoe - 0.5 * rho * (u*u + v*v)) * (M_gamma - 1);
            
            // Calculate speed of sound
            double a = sqrt(M_gamma * p / rho);
            
            // Calculate local maximum wave speeds
            double local_S_plus_x = fabs(u) + a;
            double local_S_plus_y = fabs(v) + a;
            
            // Update global maximum wave speeds through reduction
            if (local_S_plus_x > S_plus_x) S_plus_x = local_S_plus_x;
            if (local_S_plus_y > S_plus_y) S_plus_y = local_S_plus_y;
        }
    }
    
    // Calculate Δt for each direction
    double delta_tx = dx / S_plus_x;
    double delta_ty = dy / S_plus_y;
    
    // Calculate final Δt using harmonic mean
    return CFL * (delta_tx * delta_ty) / (delta_tx + delta_ty);
}

/**
 * Calculate time step for precision test simulations
 * Uses fixed time step based on grid spacing for stability in convergence tests
 * 
 * @param rows Number of variables (unused but kept for interface consistency)
 * @param cols Number of grid points in x-direction
 * @param depth Number of grid points in y-direction
 * @param x Array of conservative variables [4][cols][depth]
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @return Fixed time step Δt = 0.5 * dx^(5/3)
 */
static inline double Precision_Get_Delta_T_2D(int rows, int cols, int depth, double (*x)[cols][depth], double dx, double dy) {
    return 0.5 * pow(dx, 5.0/3);
}

/**
 * Calculate total conserved quantities in the domain
 * Sums up cell values multiplied by cell volume (Δx)
 * 
 * @param rows Number of conserved variables
 * @param cols Number of grid points including ghost cells
 * @param Ghost_Cell Number of ghost cells on each side
 * @param x Array of conserved variables [rows][cols]
 * @param Conser Output array for total conserved quantities [rows]
 * @param delta_x Grid spacing in x-direction
 */
static inline void Total_Conser(int rows, int cols, double Ghost_Cell, double (*x)[cols], double Conser[rows], double delta_x) {
    for (int i = 0; i < rows; i++){
        Conser[i] = 0.0;  // Initialize to zero
        for (int j = Ghost_Cell; j < cols - Ghost_Cell; j++){
            Conser[i] += x[i][j] * delta_x;  
        } 
    }
}



/**
 * Trapezoidal rule for numerical integration
 * Approximates ∫f(x)dx from x_down to x_up using trapezoidal rule
 * 
 * @param f Function to integrate (takes x and additional parameter u)
 * @param x_down Lower integration limit
 * @param x_up Upper integration limit
 * @param k Number of subdivisions
 * @param u Additional parameter passed to function f
 * @return Numerical approximation of the integral
 */
static inline double trapezoid_rule(FuncPtrWithU f, double x_down, double x_up, int k, double* u) {
    if (k <= 0) {
        printf("Error: n must be positive\n");
        return 0.0;
    }
    
    double h = (x_up - x_down) / k;                    // Step size
    double sum = 0.5 * (f(x_down, u) + f(x_up, u));    // Endpoints
    
    for (int i = 1; i < k; i++) {
        double x = x_down + i * h;
        sum += f(x, u);
    }
    
    return sum * h;
}

/**
 * Simpson's rule for numerical integration
 * Approximates ∫f(x)dx from x_down to x_up using Simpson's 1/3 rule
 * 
 * @param f Function to integrate (takes x and additional parameter u)
 * @param x_down Lower integration limit
 * @param x_up Upper integration limit
 * @param k Number of subdivisions (must be even)
 * @param u Additional parameter passed to function f
 * @return Numerical approximation of the integral
 */
static inline double simpson_rule(FuncPtrWithU f, double x_down, double x_up, int k, double* u) {
    if (k <= 0 || k % 2 != 0) {
        printf("Error: n must be positive and even\n");
        return 0.0;
    }
    
    double h = (x_up - x_down) / k;                    // Step size
    double sum = f(x_down, u) + f(x_up, u);           // Endpoints
    
    // Apply Simpson's rule coefficients: 4 for odd indices, 2 for even indices
    for (int i = 1; i < k; i++) {
        double x = x_down + i * h;
        if (i % 2 == 0) {
            sum += 2.0 * f(x, u);
        } else {
            sum += 4.0 * f(x, u);
        }
    }
    
    return sum * h / 3.0;
}

/**
 * Rectangle rule (midpoint rule) for numerical integration
 * Approximates ∫f(x)dx from x_down to x_up using midpoint rule
 * 
 * @param f Function to integrate (takes x and additional parameter u)
 * @param x_down Lower integration limit
 * @param x_up Upper integration limit
 * @param k Number of subdivisions
 * @param u Additional parameter passed to function f
 * @return Numerical approximation of the integral
 */
static inline double rectangle_rule(FuncPtrWithU f, double x_down, double x_up, int k, double* u) {
    if (k <= 0) {
        printf("Error: n must be positive\n");
        return 0.0;
    }
    
    double h = (x_up - x_down) / k;                    // Step size
    double sum = 0.0;
    
    for (int i = 0; i < k; i++) {
        double x_mid = x_down + (i + 0.5) * h;        // Midpoint of each subinterval
        sum += f(x_mid, u);
    }
    
    return sum * h;
}

/**
 * Main numerical integration function
 * Selects integration method based on parameter
 * 
 * @param f Function to integrate (takes x and additional parameter u)
 * @param a Lower integration limit
 * @param b Upper integration limit
 * @param n Number of subdivisions
 * @param method Integration method: 1=Trapezoidal, 2=Simpson, 3=Rectangle
 * @param u Additional parameter passed to function f
 * @return Numerical approximation of the integral
 */
static inline double numerical_integration(FuncPtrWithU f, double a, double b, int n, int method, double* u) {
    switch (method) {
        case 1:
            return trapezoid_rule(f, a, b, n, u);
        case 2:
            return simpson_rule(f, a, b, n, u);
        case 3:
            return rectangle_rule(f, a, b, n, u);
        default:
            printf("Error: Unknown method. Using trapezoid rule.\n");
            return trapezoid_rule(f, a, b, n, u);
    }
}

/**
 * Smooth function for calculating fifth-order scheme numerical error
 * Used in convergence analysis to compare exact and numerical solutions
 * 
 * @param x Coordinate at which to evaluate function
 * @param u Array of cell averages for polynomial reconstruction
 * @return Difference between exact solution and reconstructed polynomial value
 */
static inline double Smooth_function(double x, double* u) {
    // Extract neighboring cell averages for polynomial reconstruction
    double v1 = *(u-2);
    double v2 = *(u-1);
    double v3 = *(u);
    double v4 = *(u+1);
    double v5 = *(u+2);

    // Calculate polynomial coefficients for fifth-order reconstruction
    double a4 = (v1+v5 - 4.0*(v2+v4) + 6.0*v3)/24.0;
    double a3 = (v5-v1 - 2.0*(v4-v2))/12.0;
    double a2 = (v4+v2)/2.0 - v3 - 3.0*a4/2.0;
    double a1 = (v4-v2)/2.0 - 5.0*a3/4.0;
    double a0 = v3 - a2/12.0 - a4/80.0;

    // Return difference between exact solution and reconstructed polynomial
    // For testing with smooth exact solution: sin(4πx) + 2
    return sin(4.0*PI*x) + 2 - a4*pow(x,4) - a3*pow(x,3) - a2*pow(x,2) - a1*x - a0;
}


#endif