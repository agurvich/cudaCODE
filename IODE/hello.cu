// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
//#include <cutil.h>

// includes, kernels
#include "hello.h"

////////////////////////////////////////////////////////////////////////////////
// declarations, forward

//extern "C"

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    // number of ode systems ("elements"), e.g. 10 million
    int numODE = 1;

    // number of equations, e.g. 157
    int NEQN = 2;

    // the actual equations
    //double ** y = 
    //double * g = 

    double y[2];
    double g[2];
    y[0] = 2.4;
    y[1] = 0.2; // cycles per second, matches spring constant

    g[0] = 0;
    g[1] = 0.2;

    double tEnd = 10;//seconds

    double t0 = 0;
    double h = 0.1;// seconds

    // Format host matrix into 1-d array
    double * yHost ;
    yHost = ( double *) malloc ( numODE * NEQN * sizeof ( double ));
    yHost[0] = y[0];
    yHost[1] = y[1];

    /*

    for (int i = 0; i < numODE ; ++i) {
        for (int j = 0; j < NEQN ; ++j) {
            yHost [i + numODE * j] = y[i][j];
        }
    }
    */

    // allocate memory on the device
    double * yDevice ;
    cudaMalloc (( void **) & yDevice , numODE * NEQN * sizeof ( double ));
    
    double * gDevice ;
    cudaMalloc (( void **) & gDevice , numODE * sizeof ( double ));

    // setup grid dimensions
    int blockSize ;
    if ( numODE < 4194304) {
        blockSize = 64;
    } 
    else if ( numODE < 8388608) {
        blockSize = 128;
    } 
    else if ( numODE < 16777216) {
        blockSize = 256;
    }
    else {
        blockSize = 512;
    }
    dim3 dimBlock ( blockSize , 1);
    dim3 dimGrid ( numODE / dimBlock .x, 1);

    // set initial time
    double t = t0;
    double tNext = t + h;
    
    while (t < tEnd ) {
        // transfer memory to GPU
        cudaMemcpy ( yDevice , yHost , numODE * NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
        
        intDriver <<<dimGrid , dimBlock >>> (t, tNext , numODE , NEQN, gDevice , yDevice );
        
         // transfer memory back to CPU
        cudaMemcpy (yHost , yDevice , numODE * NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );
        
        for (int j=0; j<NEQN; j++){
            printf("r = %.2f",y[0]);
            printf("phi = %.2f",y[1]);
        }
        
        t = tNext ;
        tNext += h;
    }
    
     cudaFree ( gDevice );
     cudaFree ( yDevice );
}

