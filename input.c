#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*The format of the elements should be as follows: 

Y values | G Value  \n
.....				\n
\n

*/


void newElement(char* line, double * y, double * g, int* NEQN, int whichEle){
    char* tok;
    unsigned loc = 0;
    unsigned loc2 = 0;
	//printf("Handling line %d\n", whichEle);
	for (tok = strtok(line, "\t"); tok && *tok; tok = strtok(NULL, "\t\n")){
		
		if(loc < (*NEQN)){
			//printf("Dealing with %s and placing it in index %d * %d + %d = %d in Y\n", tok, whichEle, *NEQN, loc, (whichEle * (*NEQN)) + loc);
			y[(whichEle * (*NEQN)) + loc] = atof(tok);
			loc++;
		}
		else{
			//printf("Dealing with %s and placing it in index %d * %d + %d = %d in G\n", tok, whichEle, *NEQN, loc2, (whichEle * (*NEQN)) + loc2);
			g[(whichEle * (*NEQN)) + loc2] = atof(tok);
			loc2++;
		}
	}
}

void getStats(char* inputFile,int* NEQN, int* numODE){
	FILE* oStream = fopen(inputFile, "r");
	*NEQN = 0;
	*numODE = 0;
	char firstLineDone = 0;
	char c;
	unsigned count = 0;
	//printf("Opening file %s, numODE = %d, NEQN = %d \n", inputFile, *numODE, *NEQN);
	for (c = getc(oStream); c != EOF; c = getc(oStream)){
        if (c == '\n'){ // Increment count if this character is newline
        	//printf("Found a new line!\n");
            *numODE = *numODE + 1;
            firstLineDone = 1;
        }
        if (c == '\t' && !firstLineDone){
        	count++;
        }
    }
    fclose(oStream);
    //printf("%d Number of Elements\n", *numODE);
    (*NEQN) = count/2;
    //printf("There are %d equations per element, count = %d\n", *NEQN, count);
}

void parseInputs(char* inputFile, double * y, double * g, int* NEQN, int* numODE){
	
	char line[1024];
	int firstLine = 1;
	int whichEle = 0;

    FILE* stream = fopen(inputFile, "r");
	//This will seqentially go through each line of the input file and then create a new element per-line.
	while(fgets(line, 1024, stream)){
		char* tmp = strdup(line);
		
			//actually create the element
			//printf("Making new element\n");
			newElement(tmp, y, g, NEQN, whichEle);
			//Increment which element to work on
			whichEle++;
		
	}

	fclose(stream);
	return;
}
