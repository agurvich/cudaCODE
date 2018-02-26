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


#define EQUIL_LENGTH = 2
#define SPRING_CONSTANT = .04 // 1/25 -> w = 0.2 cycles per second

__device__ void dydt(double t, double * y, double * g, double * F, int NEQN){
    // r equation
    //F[0] = g[1]*SPRING_CONSTANT*(y[0]-EQUIL_LENGTH) + g[0];
    F[0] = 0.2*2*(y[0]-.04);

    // phi equation
    F[1] = g[1];
}




 __device__ void riemannStep(double * y, double * F, double h, double * yTemp, double *  yErr,int NEQN){
    for (int i=0; i < NEQN; i++){
       yTemp[i]=y[i]+F[i]*h;
       // assume riemann is super awesome
       yErr[i]=0;
    }
 }

__device__ void rk4Step(double * y, double * F, double h, double * yTemp, double *  yErr,int NEQN){
    for (int i=0; i < NEQN; i++){

    /*
        k0 = F[i];
        k1 = <<<launch a new kernel?>>>dydt(t,y+k0*h/4)
        k2 = dydt(t,y+k1*h/2)
        k3 = dydt(t,y+k2*3/4*h)
        
        yTemp[i]=y[i]+(k0 + 2*k1 +2*k2 + k3)*h;

        yErr[i]=0;

    */
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
	
 	// integrate until specified end time
        while (t < tEnd ) {
                // for implicit solver, there would be a step size estimate here
                // maybe every 25 steps or something so overhead is invested so that
                // you don't reject steps as often. 
		
 		// limit step size based on remaining time
 		h = fmin ( tEnd - t, h);
		
 		// y and error vectors temporary until error determined
 		double yTemp [ 2 ], yErr [ 2 ];
		
 		// evaluate derivative
 		double F[ 2 ];

                // function that takes the state and finds the derivative at the current time
                // for coupling just look within y!!! don't need to couple elements *yet*


 		//dydt (t, y, g, F, NEQN);
///////////////////////////////// here is where you define your equations, one for each of NEQN
                F[0] = 1;
                F[1] = (y[0]+g[1]*t)*y[1];
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
 		riemannStep (y, F, h, yTemp , yErr, NEQN);
 		//rk4Step (t, y, g, F, h, yTemp , yErr );
 		//rkckStep (t, y, g, F, h, yTemp , yErr );
		
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
 		err /= eps ;
		
 		// check if error too large
 		if (( err > 1.0) || isnan (err ) || ( nanFlag == 1)) {
 		// step failed , error too large
 			if ( isnan (err ) || ( nanFlag == 1)) {
 				h *= P1;
 			} 
 			else {
 				h = fmax ( SAFETY * h * pow(err , PSHRNK ), P1 * h);
 			}
 		} 
 		else {
 			// step accepted
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



