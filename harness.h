
__device__ void dydt(double t, double * y, double * g, double * F, int NEQN);

__global__ void intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            double * gGlobal , double * yGlobal, int method_flag ) ;
void parseInputs(char* inputFile, double ** y, double * g, int* NEQN, int* numODE);

