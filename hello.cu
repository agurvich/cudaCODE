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
    printf("Hello world\n");

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

