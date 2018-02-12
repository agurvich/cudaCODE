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
    int numODE = 

    // number of equations, e.g. 157
    int NEQN = 

    // the actual equations
    double ** y = 
    double * g = 

    // Format host matrix into 1-d array
    double * yHost ;
    yHost = ( double *) malloc ( numODE * NEQN * sizeof ( double ));

    for (int i = 0; i < numODE ; ++i) {
        for (int j = 0; j < NEQN ; ++j) {
            yHost [i + numODE * j] = y[i][j];
        }
    }

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
        
        t = tNext ;
        tNext += h;
    }
    
     cudaFree ( gDevice );
     cudaFree ( yDevice );




    int size = 32;
    int memsize = size*sizeof(int);

    int* A;
    int* B;
    int* C;
    A = (int*) malloc(size*sizeof(int));
    B = (int*) malloc(size*sizeof(int));
    C = (int*) malloc(size*sizeof(int));

    int* d_A;
    int* d_B;
    int* d_C;

    cudaMalloc((void**)&d_A, memsize);
    cudaMalloc((void**)&d_B, memsize);
    cudaMalloc((void**)&d_C, memsize);

    // set values of host arrays
    for (int j = 0; j<size; j++){
        A[j]=j;
        B[j]=2*j;
    }

    // copy host arrays to device
    cudaMemcpy((void *)d_A,A, memsize,cudaMemcpyHostToDevice);
    cudaMemcpy((void *)d_B,B, memsize,cudaMemcpyHostToDevice);
    cudaMemcpy((void *)d_C,C, memsize,cudaMemcpyHostToDevice);


    // invoke kernel
    dim3 threadsPerBlock(32,1);
    printf("Starting GPU evaluation\n");
    HelloKernel<<<1,threadsPerBlock>>>( d_A, d_B, d_C);

    // wait for threads to finish
    cudaThreadSynchronize();
    printf("\nFinished GPU evaluation\n");

    // grab output
    cudaMemcpy((void *)C,(void *) d_C, memsize, 
                    cudaMemcpyDeviceToHost);

    for (int i =0; i<size; i++){
        printf("%d ",C[i]);
    }
    printf("\n");

    cudaFree((void *)d_A);
    cudaFree((void *)d_B);
    cudaFree((void *)d_C);
}

