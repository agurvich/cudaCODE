#include <stdio.h>
#include <math.h>
#include "harness.h"

# define UROUND (2.22e-16)
# define SAFETY 0.9
# define PGROW ( -0.2)
# define PSHRNK ( -0.25)
# define ERRCON (1.89e-4)
# define TINY (1.0e-30)
# define STEP_INCREASE 1.0
# define P1 0.9


#define EQUIL_LENGTH  2
#define SPRING_CONSTANT  .04 // 1/25 -> w = 0.2 cycles per second
#define MASS  1
#define NUM_OF_EQUATIONS_PER_ELEMENT 3
const float eps = 1.0e-10;

__device__ void dydt(double t, double * y, double * g, double * F, int NEQN){

    // test set of equations that have closed forms
    F[0] = g[0];//y[0] = g[0]t+y00
    F[1] = -cosf((float) t* g[1]);//y[1] = -1/g[1]sin(g[1] t) + y10
    F[2] = (y[0] + y[1])*g[2]; //y[2] = (g[0]t^2 + (y00+y10)t + 1/g[1]^2cos(g[1]t))*g[2] + y20
}

 __device__ void riemannStep(double * y, double * F, double h, double * yTemp, double *  yErr,int NEQN){
    for (int i=0; i < NEQN; i++){
       yTemp[i]=y[i]+F[i]*h;
       // assume riemann is super awesome
       yErr[i]=0;
    }
 }

__device__ void rk4Step(double t, double * y, double * F, double h, double* g,  double * yTemp, double *  yErr,int NEQN){
    //USING THIS: https://rosettacode.org/wiki/Runge-Kutta_method#C
    double tempF[NUM_OF_EQUATIONS_PER_ELEMENT];
    double tempY[NUM_OF_EQUATIONS_PER_ELEMENT];
    for (int i=0; i < NEQN; i++){
        
        //calculate k0
        double k0 = h * F[i];

        //calculate k1
        tempY[i] = y[i] + k0 / 2;
        dydt(t+h/2,tempY, g, tempF, NEQN);
        double k1 = h * tempF[i];

        //calculate k2
        tempY[i] = y[i] + k1 / 2;
        dydt(t+h/2,tempY, g, tempF, NEQN);
        double k2 = h * tempF[i];

        //calculate k3
        tempY[i] = y[i] + k2;
        dydt(t+h, tempY, g, tempF, NEQN);
        double k3 = h * tempF[i];

        //These K-values are for error checking
        tempY[i] = y[i] + (7*k0 + 10*k1 + k3)/27;
        dydt(t + 2*h/3, tempY, g, tempF, NEQN);
        double k4 = h * tempF[i];

        tempY[i] = y[i] + (28*k0 - 125*k1 + 546*k2 + 54*k3 - 378*k4)/625;
        dydt(t + 2*h/10, tempY, g, tempF, NEQN);
        double k5 = h * tempF[i];

        double rk5Solution = y[i] + (14*k0 + 35*k3 + 162*k4 + 125*k5)/336;

        // combine using RK4 weighting
        yTemp[i]=y[i]+(k0 + 2*k1 +2*k2 + k3)/6;

        
        //printf("y[i](%f)/y2(%f) - 1 = %f\n", y[i], y2, y[i]/y2 - 1);


        ////printf("About to calculate error using h = %f and t = %f\n meaning my timeTemp should be %f, but it is %f", h, t, h + t, timeTemp);
        
        //taking absolute difference between rk4 and riemann to determine eror
        yErr[i] = fabs(yTemp[i] - rk5Solution); //fabs(yTemp[i] - (y[i] + F[i]*h));
        //printf("Error is %f\n", yErr[i]);
        //yErr[i]= yTemp[i]/(pow(((timeTemp)*(timeTemp)) / 4 + 1, 2)) - 1;
        //yErr[i] = (y[i]/y2) - 1;
        //printf("yErr[%d] = %f (abs(yTemp[%d](%f) - (y[%d](%f) + F[%d](%f)*h(%f)))\n", i, yErr[i], i, yTemp[i], i, y[i], i, F[i],h);
        ////printf("I calculated the error to be:\n %f(yErr[%d]) = %f(yTemp[%d]) / (pow(%d * %d (tempTime) / 4 + 1, 2)) - 1\n", yErr[i], i, yTemp[i], i, timeTemp, timeTemp);
    }
}

__device__ void
 innerstep ( double t, const double tEnd , double * g,
 double * y, int NEQN, int method_flag) {
    int tid = threadIdx.x + ( blockDim.x * blockIdx.x);
    // maximum and minimum allowable step sizes
    const double hMax = fabs ( tEnd - t);
    const double hMin = 1.0e-20;

    // initial step size
    double h = 0.5 * fabs ( tEnd - t);
    // integrate until specified end time
        while (t < tEnd ) {
        // for implicit solver, there would be a step size estimate here
        // maybe every 25 steps or something so overhead is invested so that
        // you don't reject steps as often. 
    
        // limit step size based on remaining time
        h = fmin ( tEnd - t, h);
        
        // y and error vectors temporary until error determined
        double yTemp [ NUM_OF_EQUATIONS_PER_ELEMENT ], yErr [ NUM_OF_EQUATIONS_PER_ELEMENT ];
        
        // evaluate derivative
        double F[ NUM_OF_EQUATIONS_PER_ELEMENT ];


        //Ok so basically the y and g pointer is being incremented 
        //to the start of where that thread is going to start looking.
        //This can be changed, but isn't necessary right now.
        dydt (t, (y + (tid * NEQN)), (g + (tid * NEQN)), F, NEQN);

        // take a trial step
        if (method_flag == 0) riemannStep ((y + (tid * NEQN)), F, h, yTemp , yErr, NEQN);
        else if (method_flag ==1) rk4Step(t,
            (y + (tid * NEQN)), F, h, (g + (tid * NEQN)), yTemp , yErr, NEQN);

        // calculate error
        double err = 0.0;
        int nanFlag = 0;

        for (int i = 0; i < (NEQN) ; ++i) {
            if ( isnan ( yErr [i])) nanFlag = 1;
            // when we take the max, we necessarily limit each set of equations to proceed
            // at the rate of the *slowest* equation. If we relax this constraint (and maybe 
            // give each equation its own thread, we could see some serious optimization if 
            // shared memory was used correctly

            //Brian changed this
            //err = fmax(err, yErr[i]);
            ////printf("yErr = %f, y = %f, h = %f, f = %f\n", yErr[i], y[i], h, F[i] );
            //err = fmax(err, fabs(yErr[i]));
            err = fmax (err , yErr[i]);

        }
        //printf("Error is calculated to be %f\n", err);

        //This might be needed, but I am killing it for now
        //err /= eps ;
        
        // check if error too large
        if (( err > 0.000001) || isnan (err ) || ( nanFlag == 1)) {
            //printf("Failed error check, error is %f\n", err);
            //printf("Error is too large / wrong\n");
        // step failed , error too large
            if ( isnan (err ) || ( nanFlag == 1)) {
                h *= P1;
            } 
            else {
                //h = fmax ( SAFETY * h * pow(err , PSHRNK ), P1 * h);
                h *= P1;
            }
            if(h < 0.0000000001){
                //printf("h has bottomed out\n");
                h = .00001;
            }
            //printf("h is now %.16f\n", h);
        } 
        else {
            // step accepted
            //printf("Step accepted, t = %f will now be %f\n", t, t + h);
            t += h;
            if ( err > ERRCON ) {
                //h = SAFETY * h * pow(err , PGROW );
            }
            else {
                h *= STEP_INCREASE;
            }
            // ensure step size is bounded
            h = fmax (hMin , fmin (hMax , h));
            for (int i = (NEQN * tid); i < (NEQN*tid + NEQN) ; ++i)
                y[i] = yTemp [i - NEQN*tid];
        }
    }
    __syncthreads();
 }

__global__ void
intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            double * gGlobal , double * yGlobal, int method_flag ) {

    // unique thread ID , based on local ID in block and block ID
    int tid = threadIdx.x + ( blockDim.x * blockIdx.x);

    // ensure thread within limit
    if (tid < numODE ) {
        innerstep(t, tEnd , gGlobal, yGlobal, NEQN, method_flag);
     }
}



