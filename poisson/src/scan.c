#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void getPECount(double * x1, double * y1, double * x2, double * y2, int * nPop, int * nEvent, int locCount, double wSize, int wCount, int * popInW, int * eventInW, double elimIntersectOD) {
	double distance;
	int minWindow;

	for(int i = 0; i < locCount * wCount; i++) {
		popInW[i] = 0;
		eventInW[i] = 0;
	}

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]) + (x2[i] - x2[j]) * (x2[i] - x2[j]) + (y2[i] - y2[j]) * (y2[i] - y2[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				popInW[i * wCount + k] += nPop[j];
				eventInW[i * wCount + k] += nEvent[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((x1[i] - x2[i]) * (x1[i] - x2[i]) + (y1[i] - y2[i]) * (y1[i] - y2[i])) / elimIntersectOD;
			int maxWindow = (int)(ceil(ODDistance / wSize));
			if(maxWindow > 0)
				maxWindow --;
			for(int k = maxWindow; k < wCount; k++) {
				popInW[i * wCount + k] = -1;
				eventInW[i * wCount + k] = -1;
			}
		}
	}

	return;
}

void getECountOnly(double * x1, double * y1, double * x2, double * y2, int * nEvent, int locCount, double wSize, int wCount, int * eventInW, double elimIntersectOD) {
	double distance;
	int minWindow;

	for(int i = 0; i < locCount * wCount; i++) {
		eventInW[i] = 0;
	}

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]) + (x2[i] - x2[j]) * (x2[i] - x2[j]) + (y2[i] - y2[j]) * (y2[i] - y2[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				eventInW[i * wCount + k] += nEvent[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((x1[i] - x2[i]) * (x1[i] - x2[i]) + (y1[i] - y2[i]) * (y1[i] - y2[i])) / elimIntersectOD;
			int maxWindow = (int)(ceil(ODDistance / wSize));
			if(maxWindow > 0)
				maxWindow --;
			for(int k = maxWindow; k < wCount; k++) {
				eventInW[i * wCount + k] = -1;
			}
		}
	}

	return;
}


void loglikelihood(double * ll, int * popInW, int * eventInW, int totalWindow, int popCount, int eventCount, int highLow) {
	
	double llTemp;
	double event, pop, expEvent;
	bool highCluster = true;
	bool lowCluster = true;
	if(highLow == 1)
		lowCluster = false;
	else if(highLow == -1)
		highCluster = false;


#pragma omp parallel for private(event, pop, expEvent, llTemp)
	for(int i = 0; i < totalWindow; i++) {

		pop = popInW[i];
		event = eventInW[i];
		expEvent = pop * eventCount / popCount;		
		
		if(pop < 1 || pop > popCount - 2) {
			ll[i] = -9999;
		}
		else if(event >= expEvent) { //High cluster of events
			if(highCluster) {
				llTemp = event * log(event/expEvent);
				if(event < eventCount) {
					llTemp += (eventCount - event) * log((eventCount - event)/(eventCount - expEvent));
				}
				ll[i] = llTemp;
			}
			else
				ll[i] = -9999;
		}
		else { //Low cluster of cases
			if(lowCluster) {
				llTemp = (eventCount - event) * log((eventCount - event)/(eventCount - expEvent));
				if(event > 0) {
					llTemp += event * log(event/expEvent);
				}
				ll[i] = llTemp;
			}
			else
				ll[i] = -9999;
		}


	}			
}


void findTopNCluster(double * x1, double * y1, double * x2, double * y2, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int aCenter = -1;
	int aRadius = -1;

	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < wCount; j++) {
			if(ll[i * wCount + j] > -9990) {
				if(aCenter < 0) {
					aCenter = i;
					aRadius = j;
				}
				else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
					aCenter = i;
					aRadius = j;
				}
			}
		}
	}

	center[0] = aCenter;
	radius[0] = aRadius;	
	cLL[0] = ll[aCenter * wCount + aRadius];

	double lastX1, lastY1, lastX2, lastY2, lastRad;
	lastX1 = x1[aCenter];
	lastY1 = y1[aCenter];
	lastX2 = x2[aCenter];
	lastY2 = y2[aCenter];

	lastRad = (aRadius + 1) * wSize;

	double distance;
	int maxWindow;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < locCount; i++) {
			distance = sqrt((x1[i] - lastX1) * (x1[i] - lastX1) + (y1[i] - lastY1) * (y1[i] - lastY1) + (x2[i] - lastX2) * (x2[i] - lastX2) + (y2[i] - lastY2) * (y2[i] - lastY2)) - lastRad;
			maxWindow = ceil(distance / wSize) - 1;
			if(maxWindow < 0)
				maxWindow = 0;
			for(int j = maxWindow; j < wCount; j++) {
				ll[i * wCount + j] = 1;
			}			
		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < locCount; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] > -9990) {
					if(aCenter < 0) {
						aCenter = i;
						aRadius = j;
					}
					else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
						aCenter = i;
						aRadius = j;
					}
				}
			}
		}
		center[c] = aCenter;
		radius[c] = aRadius;
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * wCount + aRadius];
		}
		else {
			break;
		}

		lastX1 = x1[aCenter];
		lastY1 = y1[aCenter];
		lastX2 = x2[aCenter];
		lastY2 = y2[aCenter];

		lastRad = (aRadius + 1) * wSize;

	}

	return;	
}
