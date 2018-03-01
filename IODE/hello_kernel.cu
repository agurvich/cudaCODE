#include <stdio.h>
#include <math.h>
#include "hello.h"

# define UROUND (2.22e-16)
# define SAFETY 0.9
# define PGROW ( -0.2)
# define PSHRNK ( -0.25)
# define ERRCON (1.89e-4)
# define TINY (1.0e-30)
# define STEP_INCREASE 1.0
# define P1 0.9

const float eps = 1.0e-10;


#define EQUIL_LENGTH 2
#define SPRING_CONSTANT .04 // 1/25 -> w = 0.2 cycles per second

__device__ void dydt(double t, double * y, double * g, double * F, int NEQN){
    //printf("I have entered dydt\n");
    // r equation
    // F[0] = g[1]*SPRING_CONSTANT*(y[0]-EQUIL_LENGTH) + g[0];
    // //F[0] = .2*2*(y[0]-.04);

    // // phi equation
    // F[1] = g[1];

    F[0] = g[0];
    F[1] = (y[0]+g[1]*t)*y[1];
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
    double tempF[2];
    double tempY[2];
    //double tempTime = t + h;
    //printf("Entering rk4Step for time = %f and h = %f\n", t, h);
    for (int i=0; i < NEQN; i++){
        
        double k0 = h * F[i];
        //printf("Calculated k0 to be %f = %f (h) * %f (F[%d])\n", k0, h, F[i], i);
        __syncthreads();
        tempY[0] = y[0] + k0 / 2;
        tempY[1] = y[1] + k0 / 2;
        //f(t+ (h/2), y[i] + k0 / 2
        dydt(t+h/2,tempY, g, tempF, NEQN);
        __syncthreads();
        double k1 = h * tempF[i];
        

        tempY[0] = y[0] + k1 / 2;
        tempY[1] = y[1] + k1 / 2;
        //f(t + (h/2), y[i] + k1 / 2)
        dydt(t+h/2,tempY, g, tempF, NEQN);
        __syncthreads();
        double k2 = h * tempF[i];
        

        tempY[0] = y[0] + k2;
        tempY[1] = y[1] + k2;
        //f(t + h, y[i] + k2
        dydt(t+h, tempY, g, tempF, NEQN);
        __syncthreads();
        double k3 = h * tempF[i];
        
        //printf("I have calculated my k values to be %f, %f, %f, %f\n", k0,k1,k2,k3);
        __syncthreads();
        yTemp[i]=y[i]+(k0 + 2*k1 +2*k2 + k3)/6;

        //yErr[i]= yTemp[i]/(pow(((t + h)*(t + h)) / 4 + 1, 2)) - 1;
        yErr[i] = 0.0;
        __syncthreads();
        printf("I calculated the error to be:\n %f(yErr[%d]) = %f(yTemp[%d]) / (pow(%d * %d (tempTime) / 4 + 1, 2)) - 1\n", yErr[i], i, yTemp[i], i, t + h, t + h);
        __syncthreads();
    }
}

__device__ void
 rkckDriver ( double t, const double tEnd , double * g,
 double * y, int NEQN) {

 	// maximum and minimum allowable step sizes
 	const double hMax = fabs ( tEnd - t);
 	const double hMin = 1.0e-20;
	
 	// initial step size
 	double h = 0.5 * fabs ( tEnd - t);
    h = .05;
 	// integrate until specified end time
        while (t < tEnd ) {
                // for implicit solver, there would be a step size estimate here
                // maybe every 25 steps or something so overhead is invested so that
                // you don't reject steps as often. 
		
 		// limit step size based on remaining time
 		//h = fmin ( tEnd - t, h);
        __syncthreads();
        
        //h = .05;
        
        //printf("The min function determined that h is now %f\n", h);
        // __syncthreads();
        // printf("h is %f inside loop\n", h);
        
        //printf("tend (%f) - t(%f) = %f, h is now %f\n", tEnd, t, tEnd-t, h);
        // __syncthreads();
		
 		// y and error vectors temporary until error determined
 		double yTemp [ 2 ], yErr [ 2 ];
		
 		// evaluate derivative
 		double F[ 2 ];

                // function that takes the state and finds the derivative at the current time
                // for coupling just look within y!!! don't need to couple elements *yet*


 		dydt (t, y, g, F, NEQN);
///////////////////////////////// here is where you define your equations, one for each of NEQN
                //F[0] = 1;
                //F[1] = (y[0]+g[1]*t)*y[1];
///////////////////////////////// 

                /*
                if (threadIdx.x == 0) printf("t %.2f\t",t);
                if (threadIdx.x == 0) printf("y[0] %.2f\t",y[0]);
                if (threadIdx.x == 0) printf("F[0] %.2f (%.2f)\t",F[0],1.0);
                if (threadIdx.x == 0) printf("g[0] %.2f\n",g[0]);


                if (threadIdx.x == 0) printf("t %.2f\t",t);
                if (threadIdx.x == 0) printf("y[1] %.2f\t",y[1]);
                if (threadIdx.x == 0) printf("F[1] %.2f (%.2f)\t",F[1],y[0]+2*t-y[1]);
                if (threadIdx.x == 0) printf("g[1] %.2f\n",g[1]);
                */

 		// take a trial step
 		//riemannStep (y, F, h, yTemp , yErr, NEQN);
 		//printf("Jumping into rk4Step\n");


        //RK4 is not working when I just call it so I am improvising
        __syncthreads();
        rk4Step(t, y, F, h, g, yTemp , yErr, NEQN);
        __syncthreads();
        // printf("Declaring variables for rk4\n");
        // double tempF[2];
        // double tempY[2];
        // printf("Temptime delcare\n");
        // double tempTime = t + h;
        // printf("Beginning rk4 step for time %f\n", t);
        // for (int i=0; i < NEQN; i++){
            
        //     double k0 = h * F[i];
    
        //     tempY[0] = y[0] + k0 / 2;
        //     tempY[1] = y[1] + k0 / 2;
        //     //f(t+ (h/2), y[i] + k0 / 2
        //     dydt(t+h/2,tempY, g, tempF, NEQN);
        //     double k1 = h * tempF[i];
            
    
        //     tempY[0] = y[0] + k1 / 2;
        //     tempY[1] = y[1] + k1 / 2;
        //     //f(t + (h/2), y[i] + k1 / 2)
        //     dydt(t+h/2,tempY, g, tempF, NEQN);
        //     double k2 = h * tempF[i];
            
    
        //     tempY[0] = y[0] + k2;
        //     tempY[1] = y[1] + k2;
        //     //f(t + h, y[i] + k2
        //     dydt(t+h, tempY, g, tempF, NEQN);
        //     double k3 = h * tempF[i];
            
        //     yTemp[i]=y[i]+(k0 + 2*k1 +2*k2 + k3)/6;
    
        //     yErr[i]= yTemp[i]/(pow(tempTime*tempTime / 4 + 1, 2)) - 1;
            
        // }




 		//rkckStep (t, y, g, F, h, yTemp , yErr );
        //printf("Returned from rk4Step\n");
   //      for(int i = 0; i<NEQN; i++){
   //          printf("y[%d]: %f     yTemp[%d]: %f     Error is %f\n", i, y[i], i, yTemp[i], yErr[i]);		
 		// }
                // calculate error
 		double err = 0.0;
 		int nanFlag = 0;
 		for (int i = 0; i < NEQN ; ++i) {
 			if ( isnan ( yErr [i])) nanFlag = 1;
                        // when we take the max, we necessarily limit each set of equations to proceed
                        // at the rate of the *slowest* equation. If we relax this constraint (and maybe 
                        // give each equation its own thread, we could see some serious optimization if 
                        // shared memory was used correctly
 			err = fmax (err , fabs ( yErr [i] / ( fabs (y[i]) + fabs (h * F[i]) + TINY )));
 		}
        printf("Error is calculated to be %f\n", err);

        //This might be needed, but I am killing it for now
 		//err /= eps ;
		
 		// check if error too large
 		if (( err > 1.0) || isnan (err ) || ( nanFlag == 1)) {
            printf("Error is too large / wrong\n");
 		// step failed , error too large
 			if ( isnan (err ) || ( nanFlag == 1)) {
 				h *= P1;
 			} 
 			else {
 				h = fmax ( SAFETY * h * pow(err , PSHRNK ), P1 * h);
 			}
            printf("h is now %f\n", h);
 		} 
 		else {
 			// step accepted
            printf("Step accepted, t = %f will now be %f\n", t, t + h);
 			t += h;
 			if ( err > ERRCON ) {
 				h = SAFETY * h * pow(err , PGROW );
 			}
 			else {
 				h *= STEP_INCREASE;
 			}
 			// ensure step size is bounded
 			h = fmax (hMin , fmin (hMax , h));
 			for (int i = 0; i < NEQN ; ++i)
 				y[i] = yTemp [i];
 		}
 	}
 }

__global__ void
intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            double * gGlobal , double * yGlobal ) {

        
        //if (threadIdx.x == 0) printf("gGlobal[0] %.2f\n",gGlobal[0]);
        //if (threadIdx.x == 0) printf("gGlobal[1] %.2f\n",gGlobal[1]);
	// unique thread ID , based on local ID in block and block ID
	int tid = threadIdx.x + ( blockDim.x * blockIdx.x);

	// ensure thread within limit
        // (ABG): Each thread is given a system to work on
	if (tid < numODE ) {
                
                // for coupling elements, let's just do granular time steps here
                // and assume systems are constant within the integrator, will 
                // call integrator multiple times between t and tEnd
                // coupling matrix, pass stuff around as necessary
                // TODO

	 	// local array with initial values
		double yLocal [ 2 ];

	 	// constant parameter (s)
                // (ABG): doesn't need to be coalesced b.c. just 1 number 
		double gLocal[2]; 
                gLocal[0] = gGlobal[0];// [tid ];
                gLocal[1] = gGlobal[1];// [tid ];

	 	// load local array with initial values from global array
	 	for (int i = 0; i < NEQN ; ++i) {
	 		yLocal [i] = yGlobal [tid + numODE * i];
		}

	 	// call integrator for one time step
	 	rkckDriver(t, tEnd , gLocal, yLocal, NEQN);
	 	// update global array with integrated values
	 	for (int i = 0; i < NEQN ; ++i) {
	 		yGlobal [tid + numODE * i] = yLocal [i];
	 	}
	 }
}



