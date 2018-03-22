#define THREADTOELEMENTBINDING

// there can be only one
#ifndef THREADTOELEMENTBINDING
#define BLOCKTOELEMENTBINDING
#endif



__global__ void intDriver ( const double t, const double tEnd , const int numODE , 
            const int NEQN,
            double * gGlobal , double * yGlobal, int method_flag ) ;
void parseInputs(char* inputFile, double ** y, double * g, int* NEQN, int* numODE);
void transposeInputs(double * y, double * g, int NEQN, int numODE);

