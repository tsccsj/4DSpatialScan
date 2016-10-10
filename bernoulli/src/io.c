#include <stdio.h>
#include <stdlib.h>

int getNumPoints(FILE * file) {
	rewind(file);
	int count = 0;
	
	double x1, y1, x2, y2;
	int nCas, nCon;

	while(EOF != fscanf(file, "%lf,%lf,%lf,%lf,%d,%d\n", &x1, &y1, &x2, &y2, &nCas, &nCon)) {
		count ++;
	}
	
	return count;
}

void readFile(FILE * file, double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int & casCount, int & conCount) {
	rewind(file);

	int locID = 0;

	casCount = 0;
	conCount = 0;

	while(EOF != fscanf(file, "%lf,%lf,%lf,%lf,%d,%d\n", x1 + locID, y1 + locID, x2 + locID, y2 + locID, nCass + locID, nCons + locID)) {
		casCount += nCass[locID];
		conCount += nCons[locID];
		locID ++;
	}	

	return;
} 
