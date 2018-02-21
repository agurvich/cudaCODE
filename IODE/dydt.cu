#include <stdio.h>
#include <math.h>
#include "hello.h"

#define EQUIL_LENGTH = 2
#define SPRING_CONSTANT = .04 // 1/25 -> w = 0.2 cycles per second

__device__ void dydt(double t, double * y, double * g, double * F, int NEQN){
    // r equation
    //F[0] = g[1]*SPRING_CONSTANT*(y[0]-EQUIL_LENGTH) + g[0];
    F[0] = 0.2*2*(y[0]-.04);

    // phi equation
    F[1] = g[1];
}
