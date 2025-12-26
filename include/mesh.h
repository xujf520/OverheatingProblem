#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "initialize.h"

// Define block structure
typedef struct {
    int num_blocks;      // Number of blocks
    int *LNX;           // Number of interior points in x-direction for each block (array, length = num_blocks)
    int *LNY;           // Number of interior points in y-direction for each block
    int *LNX_NGC;       // Number of points including ghost cells in x-direction
    int *LNY_NGC;       // Number of points including ghost cells in y-direction
    int *start_i;       // Starting i-index in global grid for each block (interior points, starting from 0)
    int *start_j;       // Starting j-index in global grid for each block
} BlockGeometricParameters;

// Function to set up block decomposition parameters
static inline BlockGeometricParameters setup_block_decomposition(
    int test_case, int LNX, int LNY, int GC) {
    
    BlockGeometricParameters BGP;
    BGP.num_blocks = 1;         // Default: no blocking
    BGP.LNX = NULL;
    BGP.LNY = NULL;
    BGP.LNX_NGC = NULL;
    BGP.LNY_NGC = NULL;
    BGP.start_i = NULL;
    BGP.start_j = NULL;
    
    // Set block parameters based on test case
    switch (test_case) {
        case TEST_BACKWARD_STEP:
            // Backward step: divide into main region and step region (two blocks)
            BGP.num_blocks = 2;
            BGP.LNX = (int*)malloc(2 * sizeof(int));
            BGP.LNY = (int*)malloc(2 * sizeof(int));
            BGP.LNX_NGC = (int*)malloc(2 * sizeof(int));
            BGP.LNY_NGC = (int*)malloc(2 * sizeof(int));
            BGP.start_i = (int*)malloc(2 * sizeof(int));
            BGP.start_j = (int*)malloc(2 * sizeof(int));

            BGP.LNX[0] = LNX;                                                             // Main region occupies full length
            BGP.LNY[0] = LNY * 0.8;                                                       // Main region width: 80%
            BGP.LNX_NGC[0] = LNX + 2 * GC;     // Step height: 20%
            BGP.LNY_NGC[0] = LNY * 0.8 + 2 * GC;     // Step width: 20%
            BGP.start_i[0] = 0;
            BGP.start_j[0] = LNY * 0.2; // Start from right side of step
            
            // Block1
            BGP.LNX[1] = LNX * 0.2;     // Step height: 20%
            BGP.LNY[1] = LNY * 0.2;     // Step width: 20%
            BGP.LNX_NGC[1] = LNX * 0.2 + 2 * GC;     // Step height: 20%
            BGP.LNY_NGC[1] = LNY * 0.2 + 2 * GC;     // Step width: 20%
            BGP.start_i[1] = 0;
            BGP.start_j[1] = 0;
            
            break;
            
        default:
            // Case with no blocking
            BGP.num_blocks = 1;
            BGP.LNX = (int*)malloc(sizeof(int));
            BGP.LNY = (int*)malloc(sizeof(int));
            BGP.LNX_NGC = (int*)malloc(sizeof(int));
            BGP.LNY_NGC = (int*)malloc(sizeof(int));
            BGP.start_i = (int*)malloc(sizeof(int));
            BGP.start_j = (int*)malloc(sizeof(int));
            
            BGP.LNX[0] = LNX;
            BGP.LNY[0] = LNY;
            BGP.LNX_NGC[0] = LNX + 2 * GC;     // Step height: 20%
            BGP.LNY_NGC[0] = LNY + 2 * GC;     // Step width: 20%
            BGP.start_i[0] = 0;
            BGP.start_j[0] = 0;
            break;
    }
    
    return BGP;
}

// Function to free memory of block structure
static inline void free_block_decomposition(BlockGeometricParameters *BGP) {
    if (BGP->LNX) free(BGP->LNX);
    if (BGP->LNY) free(BGP->LNY);
    if (BGP->LNX_NGC) free(BGP->LNX_NGC);
    if (BGP->LNY_NGC) free(BGP->LNY_NGC);
    if (BGP->start_i) free(BGP->start_i);
    if (BGP->start_j) free(BGP->start_j);
    BGP->num_blocks = 0;
}

// Define data exchange area structure
typedef struct {
    int Overlap;      
    int *OverLap_X;           // Overlap size in x-direction
    int *OverLap_Y;           // Overlap size in y-direction
    int *start_OLX;           // Starting x-index of overlap region
    int *start_OLY;           // Starting y-index of overlap region
} DataExchangeArea;

// Function to set up overlap parameters
static inline DataExchangeArea Get_OverLap(int test_case, int LNX, int LNY, int GC) {
    
    DataExchangeArea DEA;
    DEA.Overlap = 0;  // Default: no blocking
    DEA.OverLap_X = NULL;
    DEA.OverLap_Y = NULL;
    DEA.start_OLX = NULL;
    DEA.start_OLY = NULL;
    
    // Set overlap parameters based on test case
    switch (test_case) {
        case TEST_BACKWARD_STEP:
            DEA.Overlap = 1;
            DEA.OverLap_X = (int*)malloc(sizeof(int));
            DEA.OverLap_Y = (int*)malloc(sizeof(int));
            DEA.start_OLX = (int*)malloc(sizeof(int));
            DEA.start_OLY = (int*)malloc(sizeof(int));

            DEA.OverLap_X[0] = LNX * 0.2;                                                            
            DEA.OverLap_Y[0] = 2 * GC; 
            DEA.start_OLX[0] = GC;                                                            
            DEA.start_OLY[0] = LNY * 0.2;                
            break;
            
        default:
            DEA.Overlap = 0;  // Default: no blocking
            DEA.OverLap_X = NULL;
            DEA.OverLap_Y = NULL;
            DEA.start_OLX = NULL;
            DEA.start_OLY = NULL;
    }   
    return DEA;
}

// Function to free memory of data exchange area structure
static inline void free_Get_OverLap(DataExchangeArea *DEA) {
    if (DEA->OverLap_X) free(DEA->OverLap_X);
    if (DEA->OverLap_Y) free(DEA->OverLap_Y);
    if (DEA->start_OLX) free(DEA->start_OLX);
    if (DEA->start_OLY) free(DEA->start_OLY);
    DEA->Overlap = 0;
}


//Mesh 
//geometrical parameter
static inline void Mesh_2D(int rows, int cols, double *mesh_x, double *mesh_y, double deltax, double deltay) {

    for (int j = 0; j < rows; j++) 
            mesh_x[j] = deltax * j + 0.5 * deltax;
    
    
    for (int k = 0; k < cols; k++) 
        mesh_y[k] = deltay * k + 0.5 * deltay;
        
}


//Block geometrical parameter
static inline void Block_mesh(int rows, int cols, double *mesh_x, double *mesh_y, double deltax, double deltay) {

    for (int j = 0; j < rows; j++) 
            mesh_x[j] = deltax * j + 0.5 * deltax;
    
    
    for (int k = 0; k < cols; k++) 
        mesh_y[k] = deltay * k + 0.5 * deltay;
        
}


#endif