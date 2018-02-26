#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*The format of the elements should be as follows: 

NEQN | Y values | G Value 

*/


void newElement(char* line, ** y, double * g, int* NEQN, int whichEle){
	const char* tok;
	int loc = 0;
	for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n")){
		if(loc < (*NEQN)){
			y[whichEle][loc] = atof(tok);
		}
		else{
			g[whichEle] = atof(tok);
		}
		loc++;
	}
}

//The first line of the input file will state how many equations/element and number of elements there are
void parameterSetup(char* line, double ** y, double * g, int* NEQN, int* numODE){
	const char* tok;
	int loc = 0;

	//Retrieve things needed for malloc.
	for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n")){
		//This will state how many elements in the input
		if(loc == 0){
			(*numODE) = atoi(tok);
		}
		//This will say the equations/element
		else if(loc == 1){
			(*NEQN) = atoi(tok);
		}
		loc++;
	}
	//Creates a 2D array of doubles
	double *yTemp[*numODE];
	for(int i = 0; i < *numODE; i++){
		yTemp[i] = (double *)calloc(*NEQN * sizeof(double));
	}

	//set given y to yTemp
	y = yTemp;

	//set g array
	g = (double *)calloc(NEQN * sizeof(double));

	//done
	return;

}


void parseInputs(char* inputFile, double ** y, double * g, int* NEQN, int* numODE){
	FILE* stream = fopen(inputFile, "r");
	char line[1024];
	int firstLine = 1;
	int whichEle = 0;
	//This will seqentially go through each line of the input file and then create a new element per-line.
	while(fgets(line, 1024, stream)){
		char* tmp = strdup(line);
		if(firstLine){
			parameterSetup(tmp, y, g, NEQN, numODE);
			firstLine = 0;
		}
		else{
			//actually create the element
			newElement(tmp, y, g, *NEQN, whichEle);
			//Increment which element to work on
			whichEle++;
		}
		
	}
	return;
}