#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NUM_OF_EQUATIONS_PER_ELEMENT 32

void cpu_dydt(double t, double * y, double * g, double * F,int NEQN){
    F[0] = g[0];
    F[1] = -cosf(g[1]*t);
    F[2] = (y[0] + y[1])*g[2];

}

void riemann_step(double h, double * y, double * F, int NEQN){
    for (int i=0; i<NEQN;i++){
    y[i]+=F[i]*h;
    }

}

void rk4_step(double t, double * y, double * F, double h, double* g,  double * yTemp, double *  yErr,int NEQN){
    //USING THIS: https://rosettacode.org/wiki/Runge-Kutta_method#C
    double tempF[NUM_OF_EQUATIONS_PER_ELEMENT];
    double tempY[NUM_OF_EQUATIONS_PER_ELEMENT];
    for (int i=0; i < NEQN; i++){
        
        //calculate k0
        double k0 = h * F[i];

        //calculate k1
        tempY[i] = y[i] + k0 / 2;
        cpu_dydt(t+h/2,tempY, g, tempF, NEQN);
        double k1 = h * tempF[i];

        //calculate k2
        tempY[i] = y[i] + k1 / 2;
        cpu_dydt(t+h/2,tempY, g, tempF, NEQN);
        double k2 = h * tempF[i];

        //calculate k3
        tempY[i] = y[i] + k2;
        cpu_dydt(t+h, tempY, g, tempF, NEQN);
        double k3 = h * tempF[i];

        // combine using RK4 weighting
        yTemp[i]=y[i]+(k0 + 2*k1 +2*k2 + k3)/6;

        double timeTemp = h + t;
        //testing for error
        double y2 = ((timeTemp * timeTemp) / 4) + 1.0;
        y2 = y2 * y2;
        //printf("y[i](%f)/y2(%f) - 1 = %f\n", y[i], y2, y[i]/y2 - 1);


        ////printf("About to calculate error using h = %f and t = %f\n meaning my timeTemp should be %f, but it is %f", h, t, h + t, timeTemp);
        
        //taking absolute difference between rk4 and riemann to determine eror
        yErr[i] = fabs(yTemp[i] - (y[i] + F[i]*h));
        //yErr[i]= yTemp[i]/(pow(((timeTemp)*(timeTemp)) / 4 + 1, 2)) - 1;
        //yErr[i] = (y[i]/y2) - 1;
        //printf("yErr[%d] = %f (abs(yTemp[%d](%f) - (y[%d](%f) + F[%d](%f)*h(%f)))\n", i, yErr[i], i, yTemp[i], i, y[i], i, F[i],h);
        ////printf("I calculated the error to be:\n %f(yErr[%d]) = %f(yTemp[%d]) / (pow(%d * %d (tempTime) / 4 + 1, 2)) - 1\n", yErr[i], i, yTemp[i], i, timeTemp, timeTemp);
    }
}


void CPU_intDriver(double t, double tEnd,  int numODE, int NEQN, double * g,double * y, int method_flag){
    double * F;
    double h = 0.5*(tEnd-t);
    while (t < tEnd){
        // don't overstep the boundary
        h = min(h,tEnd-t);

        // allocate derivate array
        F = (double *) malloc(sizeof(double)*NEQN);
        double yTemp [ NUM_OF_EQUATIONS_PER_ELEMENT ], yErr [ NUM_OF_EQUATIONS_PER_ELEMENT ];

        for (int i = 0; i<numODE; i++){
            // step through the list and have the pointer start at this system
            cpu_dydt(t,(y+i*NEQN),(g+i*NEQN),F,NEQN);
            // integrate 
            if (method_flag==3){
                riemann_step(h,(y+i*NEQN),F,NEQN);
            }
            if (method_flag==4){
                rk4_step(t,(y+i*NEQN), F, h, (g+i*NEQN),yTemp, yErr,NEQN);
            }
            for (int j =0; j<NEQN; j++){
                (y+i*NEQN)[j] = yTemp[j];
            }
            
        }
    
        t+=h;
    
    }

}
