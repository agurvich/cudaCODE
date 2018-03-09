#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*The format of the elements should be as follows: 

NEQN | Y values | G Value 

*/


void newElement(char* line, double * y, double * g, int* NEQN, int whichEle){
	const char* tok;
	int loc = 0;
	for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n")){
		if(loc < (*NEQN)){
			y[(whichEle * (*NEQN)) + loc] = atof(tok);
		}
		else{
			g[(whichEle * (*NEQN)) + loc] = atof(tok);
		}
		loc++;
	}
}

//The first line of the input file will state how many equations/element and number of elements there are
void parameterSetup(char* line, double * y, double * g, int* NEQN, int* numODE){
	const char* tok;
	int loc = 0;

	//Creates a 2D array of doubles
	double yTemp[(*numODE) * (*NEQN)];
	double gTemp[(*numODE) * (*NEQN)];
	

	//set given y to yTemp
	y = yTemp;

	//set g array
	g = gTemp;

	//done
	return;

}


void parseInputs(char* inputFile, double * y, double * g, int* NEQN, int* numODE){
	FILE* stream = fopen(inputFile, "r");
	FILE* oStream = fopen(inputFile, "r");
	char line[1024];
	int firstLine = 1;
	int whichEle = 0;
	uint32_t count = 0;
	char c;
	for (c = getc(oStream); c != EOF; c = getc(oStream)){
        if (c == '\n'){ // Increment count if this character is newline
            *numODE++;
        }
    }
    for (c = getc(oStream); c != '\n'; c = getc(oStream)){
        if (c == ';'){ // Increment count if this character is newline
            count++;
        }
    }
    *NEQN = count >> 1;

	//This will seqentially go through each line of the input file and then create a new element per-line.
	while(fgets(line, 1024, stream)){
		char* tmp = strdup(line);
		if(firstLine){
			parameterSetup(tmp, y, g, NEQN, numODE);
			firstLine = 0;
		}
		
			//actually create the element
			newElement(tmp, y, g, *NEQN, whichEle);
			//Increment which element to work on
			whichEle++;
		
	}
	return;
}
