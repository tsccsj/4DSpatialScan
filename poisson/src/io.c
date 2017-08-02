/**
 * io.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/01/2017}
 */

#include <stdio.h>
#include <stdlib.h>

int getNumPoints(FILE * file) {
	rewind(file);
	int count = 0;
	
	double x1, y1, x2, y2;
	int nPop, nEvent;

	while(EOF != fscanf(file, "%lf,%lf,%lf,%lf,%d,%d\n", &x1, &y1, &x2, &y2, &nPop, &nEvent)) {
		count ++;
	}
	
	return count;
}

void readFile(FILE * file, double * x1, double * y1, double * x2, double * y2, int * nPop, int * nEvent, int & popCount, int & eventCount) {
	rewind(file);

	int locID = 0;

	popCount = 0;
	eventCount = 0;

	while(EOF != fscanf(file, "%lf,%lf,%lf,%lf,%d,%d\n", x1 + locID, y1 + locID, x2 + locID, y2 + locID, nPop + locID, nEvent + locID)) {
		popCount += nPop[locID];
		eventCount += nEvent[locID];
		locID ++;
	}	

	return;
} 
