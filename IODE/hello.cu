// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
//#include <cutil.h>

// includes, kernels
#include "hello.h"
#include "input.c"
////////////////////////////////////////////////////////////////////////////////
// declarations, forward
double * allocateDeviceAndCopy(double * hostPointer,int memsize){
    double * devicePointer;
    cudaMalloc (( void **) & devicePointer , memsize);
    cudaMemcpy ( devicePointer , hostPointer ,memsize, cudaMemcpyHostToDevice );
    return devicePointer;
}

void generateSampleInput(int, int, char*);
//extern "C"

#define INPUT_FMAX 2
#define INPUT_GMAX 10
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
void generateSampleInput(int numODE,int NEQN, char* inputFile){
    FILE* stream = fopen(inputFile, "w");
    //fprintf(stream,"%d\n",NEQN);
    double value;
    for (int i=0; i<numODE; i++){
        for (int j=0; j<2*NEQN;j++){
            if (j < NEQN/2){
                // between 0 and 1
                value = (float) rand()/(float)(RAND_MAX);
                // between -0.5 and 0.5
                value-=0.5;
                value*=INPUT_FMAX;
            }
            else{ 
                // between 0 and 1
                value = (float) rand()/(float)(RAND_MAX);
                value*=INPUT_GMAX;
            }
            // print this equation's initial condition to file
            fprintf(stream,"%f\t",value);
        }
        // start a new element
        fprintf(stream,"\n");
    }
    fclose(stream);
    return;
}

int main(int argc, char** argv) {
    char* inputFile;
    char* outputFile;
    if(argc == 1){
        printf("Producing random input and output files\n");
        inputFile = "sample_input.txt";
        outputFile = "sample_output.txt";
    }
    else if(argc == 2){
        inputFile = "sample_input.txt";
        outputFile = argv[1];
        printf("Producing random input and writing it to %s\n", argv[1]);
    }
    else if(argc == 3){
        inputFile = argv[1];
        outputFile = argv[2];
        printf("Using %s as input file and %s as output file\n", inputFile, outputFile);
    }
    // number of ode systems ("elements"), e.g. 10 million
    int numODE = 4096;

    // number of equations, e.g. 157
    int NEQN = 3;

    if(argc < 3){ 
        generateSampleInput(numODE,NEQN, inputFile);
        printf("Generated Input\n");
    }

    // the actual equations' initial conditions
    double * y;
    double * g;

    //This will initiallize stuff by reading the input file.
    //This will set NEQN and numODE to be correct
    getStats(inputFile, &NEQN, &numODE);

    //printf("Setting up my y and g arrays based on a NEQN of %d and a numODE of %d\n", NEQN, numODE);
    y = (double *)malloc(sizeof(double) * (NEQN) * (numODE));
    g = (double *)malloc(sizeof(double) * (NEQN) * (numODE));


    parseInputs(inputFile, y, g, &NEQN, &numODE);
    printf("Finished filling y and g matrices with their numbers\n");
    printf("Parsed Inputs\n");
    // printf("numODE = %d, NEQN = %d\n", numODE, NEQN);
    // printf("Printing my y and g arrays\n Y: ");
    for(int k = 0; k < 16; k++){
        printf("%f ", y[NEQN*k]);
    }

    printf("\n");

    for(int k = 0; k < 16; k++){
        printf("%f ", g[NEQN*k]);
    }

    // printf("\nG: ");
    // for(int k = 0; k < NEQN * numODE; k++){
    //     printf("%f ", g[k]);
    // }
    // printf("\n\n");

    // double y[2];
    // double g[2];
    // y[0] = 0;
    // y[1] = 0; // cycles per second, matches spring constant

    // g[0] = 1;
    // g[1] = 2;

    double tEnd = 10;//seconds

    double t0 = 0;
    double h = 0.05;// seconds

    // // Format host matrix into 1-d array
    // double * yHost ;
    // yHost = ( double *) malloc ( numODE * NEQN * sizeof ( double ));
    //yHost[0] = y[0];
    //yHost[1] = y[1];

    // double * gHost;
    // gHost = (double *) malloc ( NEQN * sizeof(double));
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
    printf("Copied over my y and g to CUDA\n");
    // setup grid dimensions
    int blockSize, gridSize ;
    if ( numODE <= 1024) {
        blockSize = numODE;
        gridSize = 1;
    } 
    else{
        blockSize = 1024;
        gridSize = numODE / blockSize + 1;
    }

    //printf("%d threads/block\n",blockSize);
    dim3 dimBlock ( blockSize , 1);
    //printf("%d blocks\n",numODE/dimBlock.x);
    dim3 dimGrid ( gridSize, 1);

    // set initial time
    double t = t0;
    double tNext = t + h;
    
    FILE* outputStream = fopen(outputFile,"w");
    //printf("before intDriver %.2f %.2f\n",g[0],g[1]);
    while (t < tEnd ) {
        // write to output file before integrating, otherwise get minor phase shift
        if(1){
            //printf("System\tTime\ty0(t)\ty1(t)\n______________________________\n");
            for (int j=0; j<numODE; j++){
                fprintf(outputStream,"%.4f\t",t);
                for (int i=0; i<NEQN;i++){
                    fprintf(outputStream,"%.4f\t",y[i+NEQN*j]);
                }
                fprintf(outputStream,"\n");
            }
        }
        // transfer memory to GPU
        if (t!=t0){
            cudaMemcpy ( yDevice , y , numODE * NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
            cudaMemcpy ( gDevice , g , NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
        }

        intDriver <<<dimBlock , dimGrid >>> (t, tNext , numODE , NEQN, gDevice , yDevice );
        
         // transfer memory back to CPU
        cudaMemcpy (y , yDevice , numODE * NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );
        cudaMemcpy (g , gDevice , NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );

      
        t = tNext ;
        tNext += h;
    }
    //printf("after intDriver %.2f %.2f\n",g[0],g[1]); 
    printf("Done!\n");
    
     cudaFree ( gDevice );
     cudaFree ( yDevice );
}

