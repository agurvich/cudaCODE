#include <stdio.h>
#include "hello.h"

# define UROUND (2.22e -16)
# define SAFETY 0.9
# define PGROW ( -0.2)
# define PSHRNK ( -0.25)
# define ERRCON (1.89e -4)
# define TINY (1.0e -30)
const double eps = 1.0e -10;

__global__ void
intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            const double * gGlobal , double * yGlobal ) {
	// unique thread ID , based on local ID in block and block ID
	int tid = threadIdx.x + ( blockDim.x * blockIdx.x);

	// ensure thread within limit
        // (ABG): Each thread is given a system to work on
	if (tid < numODE ) {

	 	// local array with initial values
		double yLocal [ NEQN ];

	 	// constant parameter (s)
                // (ABG): doesn't need to be coalesced b.c. just 1 number 
		double gLocal = gGlobal [tid ];

	 	// load local array with initial values from global array
	 	for (int i = 0; i < NEQN ; ++i) {
                // (ABG): super un-coalesced easily fixed?  
	 		yLocal [i] = yGlobal [tid + numODE * i];
		}

	 	// call integrator for one time step
	 	rkckDriver(t, tEnd , yLocal , gLocal );
	 	// update global array with integrated values
	 	for (int i = 0; i < NEQN ; ++i) {
                // (ABG): super un-coalesced easily fixed?  
	 		yGlobal [tid + numODE * i] = yLocal [i];
	 	}
	 }
}


 __device__ void
 rkckDriver ( double t, const double tEnd , const double g,
 double * y) {

 	// maximum and minimum allowable step sizes
 	const double hMax = fabs ( tEnd - t);
 	const double hMin = 1.0e -20;
	
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
 		double yTemp [ NEQN ], yErr [ NEQN ];
		
 		// evaluate derivative
 		double F[ NEQN ];
                //(* TODO implement this *)
                // function that takes the state and finds the derivative at the current time
 		dydt (t, y, g, F);
		
 		// take a trial step
        //(* TODO implement this *)
 		rkckStep (t, y, g, F, h, yTemp , yErr );
		
 		// calculate error
 		double err = 0.0;
 		int nanFlag = 0;
 		for (int i = 0; i < NEQN ; ++i) {
 			if ( isnan ( yErr [i])) nanFlag = 1;
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
 				h *= 5.0;
 			}
 			// ensure step size is bounded
 			h = fmax (hMin , fmin (hMax , h));
 			for (int i = 0; i < NEQN ; ++i)
 				y[i] = yTemp [i];
 		}
 	}
 }