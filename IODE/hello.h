
__global__ void intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            const double * gGlobal , double * yGlobal ) ;

//__device__ void riemannStep(double * y, double * F, double h, double * yTemp, double *  yErr,int NEQN);
//__device__ void dydt(double t, double * y, double * g, double * F, int NEQN);
/*__device__ void
 rkckDriver ( double t, const double tEnd , const double g,
 double * y, int NEQN);
*/
