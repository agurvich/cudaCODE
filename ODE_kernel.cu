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


// flag to use RK5 to find RK4 error
//#define RK4_ERROR

#define EQUIL_LENGTH  2
#define SPRING_CONSTANT  .04 // 1/25 -> w = 0.2 cycles per second
#define MASS  1
#define NUM_OF_EQUATIONS_PER_ELEMENT 32
const float eps = 1.0e-10;

__device__ void dummy_dydt(double t, double * y, double * g, double * F, int NEQN, int numODE){
}

__device__ void dydt(double t, double * y, double * g, double * F, int NEQN, int numODE){
    // data pointers will be iterated so this thread's 0 is appropriate,
    // then need to move a stride away for the next element

#ifdef THREADTOELEMENTBINDING
    for (int i=0; i<NEQN; i++){
        // cycle through the 3 equations as necessary
        if (i%3 ==0){
            F[i*numODE] = g[i*numODE];
        }
        else if (i%3 ==1){
            F[i*numODE] = -cosf((float) t* g[i*numODE]);
        }
        else if (i%3 ==2){
            F[i*numODE] = (y[(i-2)*numODE] + y[(i-1)*numODE])*g[i*numODE]; 
        }
    }
#else
    if (threadIdx.x < NEQN){
        if (threadIdx.x%3 ==0){
            F[threadIdx.x] = g[threadIdx.x];
        }
        else if (threadIdx.x%3 ==1){
            F[threadIdx.x] = -cosf((float) t* g[threadIdx.x]);
        }
        // there are guaranteed to be 2 threads before this one
        else if (threadIdx.x%3 ==2){
            F[threadIdx.x] = (y[threadIdx.x-2] + y[threadIdx.x-1])*g[threadIdx.x]; 
        }
    }

    // make sure all the threads in the block have set the dydt
    //__syncthreads();
#endif

}

 __device__ void riemannStep(double * y, double * F, double h, double * yTemp, double *  yErr,int NEQN, int numODE){

// pointer will point to start of element
#ifdef THREADTOELEMENTBINDING
    for (int i=0; i < NEQN; i++){
        // look for next element a stride numODE away
        yTemp[i*numODE]=y[i*numODE]+F[i*numODE]*h;
        // assume riemann is super awesome
        yErr[i*numODE]=0;
    }
    //printf("Made it into riemann step! %.2f \n",h);
#else
    // each thread will work on its equation
    if (threadIdx.x < NEQN){
        yTemp[threadIdx.x]=y[threadIdx.x]+F[threadIdx.x]*h;
        // assume riemann is super awesome
        yErr[threadIdx.x]=0;
    }
    // make sure all the threads in the block have completed the step
    __syncthreads();
#endif

 }

__device__ void rk4Step(
    double t, double * y, double * F, double h, double* g,
    double * yTemp, double *  yErr, int NEQN, int numODE){
    //USING THIS: https://rosettacode.org/wiki/Runge-Kutta_method#C

#ifdef THREADTOELEMENTBINDING
    double tempF[NUM_OF_EQUATIONS_PER_ELEMENT];
    double tempY[NUM_OF_EQUATIONS_PER_ELEMENT];
    for (int i=0; i < NEQN; i++){
        //calculate k0
        double k0 = h * F[i*numODE];

        //calculate k1
        tempY[i] = y[i*numODE] + k0 / 2;
        dydt(t+h/2,tempY, g, tempF, NEQN, numODE);
        double k1 = h * tempF[i];

        //calculate k2
        tempY[i] = y[i*numODE] + k1 / 2;
        dydt(t+h/2,tempY, g, tempF, NEQN, numODE);
        double k2 = h * tempF[i];

        //calculate k3
        tempY[i] = y[i*numODE] + k2;
        dydt(t+h, tempY, g, tempF, NEQN, numODE);
        double k3 = h * tempF[i];

        //These K-values are for error checking
        tempY[i] = y[i*numODE] + (7*k0 + 10*k1 + k3)/27;
        dydt(t + 2*h/3, tempY, g, tempF, NEQN, numODE);
        double k4 = h * tempF[i];

        tempY[i*numODE] = y[i*numODE] + (28*k0 - 125*k1 + 546*k2 + 54*k3 - 378*k4)/625;
        dydt(t + 2*h/10, tempY, g, tempF, NEQN, numODE);
        double k5 = h * tempF[i];

        // combine using RK4 weighting
        yTemp[i]=y[i*numODE]+(k0 + 2*k1 +2*k2 + k3)/6;
#ifdef RK4_ERROR
        //taking absolute difference between rk4 and riemann to determine eror
        double rk5Solution = y[i*numODE] + (14*k0 + 35*k3 + 162*k4 + 125*k5)/336;
        yErr[i*numODE] = fabs(yTemp[i] - rk5Solution); 
#else
        yErr[i*numODE] = 0; 
#endif
    }

#else
    double tempF[NUM_OF_EQUATIONS_PER_ELEMENT];
    double tempY[NUM_OF_EQUATIONS_PER_ELEMENT];
    // each thread gets a number of equations
    if (threadIdx.x < NEQN){
        //calculate k0
        double k0 = h * F[threadIdx.x];

        //calculate k1
        tempY[threadIdx.x] = y[threadIdx.x] + k0 / 2;
        __syncthreads();

        dydt(t+h/2,tempY, g, tempF, NEQN, numODE);
        double k1 = h * tempF[threadIdx.x];

        //calculate k2
        tempY[threadIdx.x] = y[threadIdx.x] + k1 / 2;
        __syncthreads();

        dydt(t+h/2,tempY, g, tempF, NEQN, numODE);
        double k2 = h * tempF[threadIdx.x];

        //calculate k3
        tempY[threadIdx.x] = y[threadIdx.x] + k2;
        __syncthreads();

        dydt(t+h, tempY, g, tempF, NEQN, numODE);
        double k3 = h * tempF[threadIdx.x];

        //These K-values are for error checking
        tempY[threadIdx.x] = y[threadIdx.x] + (7*k0 + 10*k1 + k3)/27;
        __syncthreads();

        dydt(t + 2*h/3, tempY, g, tempF, NEQN, numODE);
        double k4 = h * tempF[threadIdx.x];

        tempY[threadIdx.x] = y[threadIdx.x] + (28*k0 - 125*k1 + 546*k2 + 54*k3 - 378*k4)/625;
        __syncthreads();

        dydt(t + 2*h/10, tempY, g, tempF, NEQN, numODE);
        double k5 = h * tempF[threadIdx.x];

        // combine using RK4 weighting
        yTemp[threadIdx.x]=y[threadIdx.x]+(k0 + 2*k1 +2*k2 + k3)/6;
        __syncthreads();
#ifdef RK4_ERROR
        //taking absolute difference between rk4 and riemann to determine eror
        double rk5Solution = y[threadIdx.x] + (14*k0 + 35*k3 + 162*k4 + 125*k5)/336;
        yErr[threadIdx.x] = fabs(yTemp[threadIdx.x] - rk5Solution); 
#else
        yErr[threadIdx.x] = 0; 
#endif

    }
    // make sure all the threads in the block have completed the step
    __syncthreads();

#endif

}

__device__ void
 innerstep ( double t, const double tEnd , double * g,
 double * y, int NEQN,int numODE, int method_flag) {

    int tid = threadIdx.x + ( blockDim.x * blockIdx.x);
    // maximum and minimum allowable step sizes
    const double hMax = fabs ( tEnd - t);
    const double hMin = 1.0e-20;

    // initial step size
    double h = 0.5 * fabs ( tEnd - t);

    // integrate until specified end time
    while (t < tEnd ) {
        printf("Step %.2f| %.2f \t",t,h);
        // for implicit solver, there would be a step size estimate here
        // maybe every 25 steps or something so overhead is invested so that
        // you don't reject steps as often. 
    
        // limit step size based on remaining time
        h = fmin ( tEnd - t, h);

#ifdef THREADTOELEMENTBINDING
        // y and error vectors temporary until error determined
        double yTemp [ NUM_OF_EQUATIONS_PER_ELEMENT ];
        double yErr [ NUM_OF_EQUATIONS_PER_ELEMENT ];
        
        // evaluate derivative
        double F[ NUM_OF_EQUATIONS_PER_ELEMENT ];

        // data should be stored s.t. memory access is coalesced
        // x0,x1,x2, ... | y0,y1,y2, ...
        dydt (t, (y + tid ), (g + tid ), F, NEQN, numODE);

        // take a trial step
        if (method_flag == 0) riemannStep ((y + tid), F, h, yTemp , yErr, NEQN, numODE);
        else if (method_flag ==1) rk4Step(t,
            (y + tid), F, h, (g + tid), yTemp , yErr, NEQN, numODE);
#else
        __shared__ double yTemp [ NUM_OF_EQUATIONS_PER_ELEMENT ];
        __shared__ double yErr [ NUM_OF_EQUATIONS_PER_ELEMENT ];
        
        // evaluate derivative
        double F[ NUM_OF_EQUATIONS_PER_ELEMENT ];
        // data should be stored s.t. memory acces is coalesced
        // x0,y0,z0 | x1,y1,z1 | ... 
        dydt (t, (y + (blockIdx.x * NEQN)), (g + (blockIdx.x * NEQN)), F, NEQN, numODE);

        // take a trial step
        if (method_flag == 0) riemannStep ((y + (blockIdx.x * NEQN)), F, h, yTemp , yErr, NEQN,numODE);
        else if (method_flag ==1) rk4Step(t,
            (y + (blockIdx.x * NEQN)), F, h, (g + (blockIdx.x * NEQN)), yTemp , yErr, NEQN,numODE);

#endif


#ifdef THREADTOELEMENTBINDING
        // calculate error
        double err = 0.0;
        int nanFlag = 0;

        for (int i = 0; i < (NEQN) ; ++i) {
            if ( isnan ( yErr [tid + i*numODE])) nanFlag = 1;
            // when we take the max, we necessarily limit each set of equations to proceed
            // at the rate of the *slowest* equation. If we relax this constraint (and maybe 
            // give each equation its own thread, we could see some serious optimization if 
            // shared memory was used correctly

            //Brian changed this
            //err = fmax(err, yErr[i]);
            ////printf("yErr = %f, y = %f, h = %f, f = %f\n", yErr[i], y[i], h, F[i] );
            //err = fmax(err, fabs(yErr[i]));
            err = fmax (err , yErr[tid + i*numODE]);

        }
        //printf("Error is calculated to be %f\n", err);

        //This might be needed, but I am killing it for now
        //err /= eps ;
        
        // check if error too large
        printf("error is: %.2f\n",err);
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
                //printf("step passed");
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
            for (int i = 0; i < NEQN ; ++i){
                y[tid + i*numODE] = yTemp [tid + i*numODE];
            }
        }

#else
        // automatically accept the step if we're using blocks, should do an error check 
        // and a synchronization in the future...
        y[tid]=yTemp[tid];
#endif
        __syncthreads();
    }
}

__global__ void
intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            double * gGlobal , double * yGlobal, int method_flag ) {

    // unique thread ID , based on local ID in block and block ID
    int tid = threadIdx.x + ( blockDim.x * blockIdx.x);

    // ensure thread within limit
    if (tid < numODE ) {
        innerstep(t, tEnd , gGlobal, yGlobal, NEQN, numODE, method_flag);
     }
}




