#include <stdio.h>
#include "hello.h"

__global__ void HelloKernel(int * A, int * B, int * C){
    printf("%d %d\t",A[threadIdx.x],B[threadIdx.x]);
    C[threadIdx.x]=A[threadIdx.x]+B[threadIdx.x];
}

