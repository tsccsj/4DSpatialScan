#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <omp.h>
#include "scan.h"

using namespace std;

void simulateEvents(int * nPop, int * simEvents, int locCount, int eventCount) {
	
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::discrete_distribution<int> d (nPop, nPop + locCount);

	for(int i = 0; i < locCount; i++) {
		simEvents[i] = 0;
	}
	for(int i = 0; i < eventCount; i++) {
		simEvents[d(rng)] ++;
	}
}

int * monteCarlo(double * x1, double * y1, double * x2, double * y2, int * nPop, int * popInW, int locCount, int popCount, int eventCount, double wSize, int wCount, int highLow, double elimIntersectOD, double * clusterLL, int nClusters, int nSim) {

	int * nExtreme;

	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int * simEvents;
	if(NULL == (simEvents = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int * simEventInW;
	double * simll;

	if(NULL == (simEventInW = (int *) malloc (locCount * wCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simll = (double *) malloc (locCount * wCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	double simMaxLL;

	for(int i = 0; i < nSim; i++) {
	
		simulateEvents(nPop, simEvents, locCount, eventCount);	

		getECountOnly(x1, y1, x2, y2, simEvents, locCount, wSize, wCount, simEventInW, elimIntersectOD);
		
		loglikelihood(simll, popInW, simEventInW, locCount * wCount, popCount, eventCount, highLow);

		simMaxLL = -9999; //Possion's LL is larger than 0
		
		for(int k = 0; k < locCount * wCount; k++) {
			if(simll[k] > 0 && simll[k] > simMaxLL) {
				simMaxLL = simll[k];
			}
		}

		if(simMaxLL > 0) { 
			for(int j = 0; j < nClusters; j++) {
				if(simMaxLL > clusterLL[j]) {
					nExtreme[j] ++;
				}
			}
		}
		
	}


	free(simEventInW);
	free(simll);

	free(simEvents);

	return nExtreme;

}

int * monteCarloOld(double * x1, double * y1, double * x2, double * y2, int * nPop, int locCount, int popCount, int eventCount, int * clusterEvent, int * center, double * cRadius, bool * highCluster, int nClusters, int nSim) {
	
	int * nExtreme;

	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int * simEvents;
	if(NULL == (simEvents = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nSim; i++) {

		simulateEvents(nPop, simEvents, locCount, eventCount);	
#pragma omp parallel for

		for(int j = 0; j < nClusters; j++) { 
			double x1C = x1[center[j]];
			double y1C = y1[center[j]];
			double x2C = x2[center[j]];
			double y2C = y2[center[j]];
			double rad2 = cRadius[j] * cRadius[j];
			int simEventInc = 0;

			for(int k = 0; k < locCount; k++) { 
				if((x1[k] - x1C) * (x1[k] - x1C) + (y1[k] - y1C) * (y1[k] - y1C) + (x2[k] - x2C) * (x2[k] - x2C) + (y2[k] - y2C) * (y2[k] - y2C) <= rad2) {
					simEventInc += simEvents[k];
				}
			}

			if(highCluster[j] && (simEventInc >= clusterEvent[j])) {
				nExtreme[j] ++;
			}
			else if(!highCluster[j] && (simEventInc <= clusterEvent[j])) {
				nExtreme[j] ++;
			}
			
		}

	}

		

	free(simEvents);

	return nExtreme;
}
