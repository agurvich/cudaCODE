// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
//#include <cutil.h>

// includes, kernels
#include "hello.h"

////////////////////////////////////////////////////////////////////////////////
// declarations, forward
double * allocateDeviceAndCopy(double * hostPointer,int memsize){
    double * devicePointer;
    cudaMalloc (( void **) & devicePointer , memsize);
    cudaMemcpy ( devicePointer , hostPointer ,memsize, cudaMemcpyHostToDevice );
    return devicePointer;
}

void generateSampleInput(int, int);
//extern "C"

#define INPUT_FMAX 2
#define INPUT_GMAX 10
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
void generateSampleInput(int numODE,int NEQN){
    FILE* stream = fopen("sample_input.txt", "w");
    //fprintf(stream,"%d\n",NEQN);
    double value;
    for (int i=0; i<numODE; i++){
        for (int j=0; j<NEQN;j++){
            if (j < NEQN/2){
                // between 0 and 1
                value = (float) rand()/(float)(RAND_MAX);
                value*=INPUT_FMAX;
            }
            else{ 
                // between 0 and 1
                value = (float) rand()/(float)(RAND_MAX);
                // between -0.5 and 0.5
                value-=0.5;
                value*=INPUT_GMAX;
            }
            // print this equation's initial condition to file
            fprintf(stream,"%f;",value);
        }
        // start a new element
        fprintf(stream,"\n");
    }
    return;
}

int main(int argc, char** argv) {
    // number of ode systems ("elements"), e.g. 10 million
    int numODE = 10;

    // number of equations, e.g. 157
    int NEQN = 2;
    generateSampleInput(numODE,NEQN);

    // the actual equations
    double * y;
    double * g;

    //This will initiallize stuff by reading the input file.
    parseInputs(argv[1], y, g, &NEQN, &numODE);


    // double y[2];
    // double g[2];
    // y[0] = 0;
    // y[1] = 0; // cycles per second, matches spring constant

    // g[0] = 1;
    // g[1] = 2;

    double tEnd = 14;//seconds

    double t0 = 0;
    double h = 0.75;// seconds

    // Format host matrix into 1-d array
    double * yHost ;
    yHost = ( double *) malloc ( numODE * NEQN * sizeof ( double ));
    //yHost[0] = y[0];
    //yHost[1] = y[1];

    double * gHost;
    gHost = (double *) malloc ( NEQN * sizeof(double));
    //gHost[0] = g[0];
    //gHost[1] = g[1];

    /*

    for (int i = 0; i < numODE ; ++i) {
        for (int j = 0; j < NEQN ; ++j) {
            yHost [i + numODE * j] = y[i][j];
        }
    }
    */

    // allocate memory on the device and copy over
    double * yDevice ;
    yDevice = allocateDeviceAndCopy(y,numODE * NEQN * sizeof ( double ));
    double * gDevice ;
    gDevice = allocateDeviceAndCopy(g,NEQN * sizeof ( double ));

    // setup grid dimensions
    int blockSize ;
    if ( numODE < 4194304) {
        blockSize = 64;
    } 
    else if ( numODE < 8388608) {
        blockSize = 128;
    } 
    else if ( numODE < 16777216) {
        blockSize = 256;
    }
    else {
        blockSize = 512;
    }

    //printf("%d threads/block\n",blockSize);
    dim3 dimBlock ( blockSize , 1);
    //printf("%d blocks\n",numODE/dimBlock.x);
    dim3 dimGrid ( numODE / dimBlock .x+1, 1);

    // set initial time
    double t = t0;
    double tNext = t + h;
    
    //printf("before intDriver %.2f %.2f\n",g[0],g[1]);
    while (t < tEnd ) {
        // transfer memory to GPU
        if (t!=t0){
            cudaMemcpy ( yDevice , yHost , numODE * NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
            cudaMemcpy ( gDevice , gHost , NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
        }
        
        intDriver <<<numODE , 1 >>> (t, tNext , numODE , NEQN, gDevice , yDevice );
        
         // transfer memory back to CPU
        cudaMemcpy (yHost , yDevice , numODE * NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );
        cudaMemcpy (gHost , gDevice , NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );

        // for each system
        for (int j=0; j<numODE; j++){
            printf("%.4f ",t);
            for (int i=0; i<NEQN;i++){
                printf("%.4f ",yHost[i+2*j]);
        }
        printf("\n");
    }


         
        t = tNext ;
        tNext += h;
    }
    //printf("after intDriver %.2f %.2f\n",g[0],g[1]); 
    
     cudaFree ( gDevice );
     cudaFree ( yDevice );
}

