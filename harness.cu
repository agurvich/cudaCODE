// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
//#include <cutil.h>

// includes, kernels
#include "ODE_CPU.c"
#include "harness.h"
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
    int method_flag=1;
    int out_flag = 1;
    double outer_dt=0.05;

    // number of ode systems ("elements"), e.g. 10 million
    int numODE = 2048;

    // number of equations, e.g. 157
    int NEQN = 3;

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
    else if(argc == 5){
        inputFile = argv[1];
        outputFile = argv[2];
        printf("Using %s as input file and %s as output file\n", inputFile, outputFile);
        method_flag = atoi(argv[3]);
        outer_dt = atof(argv[4]);
        printf("Using method %d with step size %.2f\n",method_flag,outer_dt);
    }

    else if(argc == 6){
        out_flag=0;
        printf("Making timing data, writing no output \n",method_flag,outer_dt);

        inputFile = argv[1];
        outputFile = argv[2];
        printf("Using %s as input file and %s as output file\n", inputFile, outputFile);

        method_flag = atoi(argv[3]);
        outer_dt = atof(argv[4]);
        printf("Using method %d with step size %.2f\n",method_flag,outer_dt);

        numODE = atoi(argv[5]);

    }

    // overwrite sampleInput
    if(argc < 3 || 1){ 
        generateSampleInput(numODE,NEQN, inputFile);
        printf("Generated Input\n");
    }

    // the actual equations' initial conditions
    double * y;
    double * g;

    //This will initiallize stuff by reading the input file.
    //This will set NEQN and numODE to be correct
    getStats(inputFile, &NEQN, &numODE);

    y = (double *)malloc(sizeof(double) * (NEQN) * (numODE));
    g = (double *)malloc(sizeof(double) * (NEQN) * (numODE));


    parseInputs(inputFile, y, g, &NEQN, &numODE);
    printf("Finished filling y and g matrices with their numbers\n");
    printf("Parsed Inputs\n");

    double tEnd = 10;//seconds

    double t0 = 0;
    double h = outer_dt;// seconds

    // allocate memory on the device and copy over
    double * yDevice ;
    yDevice = allocateDeviceAndCopy(y,numODE * NEQN * sizeof ( double ));
    double * gDevice ;
    gDevice = allocateDeviceAndCopy(g,numODE * NEQN * sizeof ( double ));
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

    dim3 dimBlock ( blockSize , 1);
    dim3 dimGrid ( gridSize, 1);

    // set initial time
    double t = t0;
    double tNext = t + h; 
    FILE* outputStream = fopen(outputFile,"w"); 
    clock_t init_time = clock();
    while (t < tEnd ) {
        // write to output file before integrating, otherwise get minor phase shift
        if(out_flag){
            for (int j=0; j<numODE; j++){
                fprintf(outputStream,"%.4f\t",t);
                for (int i=0; i<NEQN;i++){
                    fprintf(outputStream,"%.4f\t",y[i+NEQN*j]);
                }
                fprintf(outputStream,"\n");
            }
        }
        if (method_flag < 3){
            // transfer memory to GPU
            if (t!=t0){
                cudaMemcpy ( yDevice , y , numODE * NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
                cudaMemcpy ( gDevice , g , numODE * NEQN * sizeof ( double ), cudaMemcpyHostToDevice );
            }

            intDriver <<<dimBlock , dimGrid >>> (t, tNext , numODE , NEQN, gDevice , yDevice, method_flag);
            
             // transfer memory back to CPU
            cudaMemcpy (y , yDevice , numODE * NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );
            cudaMemcpy (g , gDevice , NEQN * sizeof ( double ), cudaMemcpyDeviceToHost );
        }
        else{
            CPU_intDriver(t, tNext,  numODE, NEQN,  g, y,method_flag);
        }

      
        t = tNext ;
        tNext += h;
    }
    fprintf(outputStream,"\n %f \n",1.0*(clock()-init_time)/CLOCKS_PER_SEC);
    printf("Done!\n");
    
    cudaFree ( gDevice );
    cudaFree ( yDevice );
}

