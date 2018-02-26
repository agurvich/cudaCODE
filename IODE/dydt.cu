#include <stdio.h>
#include <math.h>
#include "hello.h"

#define EQUIL_LENGTH = 2
#define SPRING_CONSTANT = .04 // 1/25 -> w = 0.2 cycles per second
#define MASS = 1

__device__ void dydt(double t, double * y, double * g, double * F, int NEQN){
    // r equation
    //F[0] = g[1]*SPRING_CONSTANT*(y[0]-EQUIL_LENGTH) + g[0];
    double dx = (y[0]-EQUIL_LENGTH);
    double dx0 = (g[0]-EQUIL_LENGTH);
    F[0] = sqrt(SPRING_CONSTANT/MASS*(dx0*dx0 - dx*dx);

    // phi equation
    F[1] = g[1];
}
