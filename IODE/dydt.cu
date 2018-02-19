#include <stdio.h>
#include "hello.h"

__global__ void dydt(double t, double * y, double * g, double * F, int NEQN){
    for (int i=0;i<NEQN;i++){
        F[i]=5;
    }
}
