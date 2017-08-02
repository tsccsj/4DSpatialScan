/**
 * scan.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/01/2017}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void getCCCount(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, double elimIntersectOD) {
	double distance;
	int minWindow;

	for(int i = 0; i < locCount * wCount; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]) + (x2[i] - x2[j]) * (x2[i] - x2[j]) + (y2[i] - y2[j]) * (y2[i] - y2[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				conInW[i * wCount + k] += nCons[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((x1[i] - x2[i]) * (x1[i] - x2[i]) + (y1[i] - y2[i]) * (y1[i] - y2[i])) / elimIntersectOD;
			int maxWindow = (int)(ceil(ODDistance / wSize));
			if(maxWindow > 0)
				maxWindow --;
			for(int k = maxWindow; k < wCount; k++) {
				casInW[i * wCount + k] = -1;
			}
		}
	}

	return;
}

void loglikelihood(double * ll, int * casInW, int * conInW, int totalWindow, int casCount, int conCount, int highLow) {
	double cas, con, tot;
	double llTemp;
	int totCount = casCount + conCount;
	bool highCluster = true;
	bool lowCluster = true;
	if(highLow == 1)
		lowCluster = false;
	else if(highLow == -1)
		highCluster = false;

#pragma omp parallel for private(cas, con, tot, llTemp)
	for(int i = 0; i < totalWindow; i++) {
		cas = casInW[i];
		con = conInW[i];
		tot = cas + con;

		if(cas == -1) {
			ll[i] = 1;
		}
		else if(cas * conCount > con * casCount) { //High cluster of cases
			if(highCluster) {
				llTemp = cas * log(cas/tot);
				if(con > 0)
					llTemp += con * log(con/tot);
				if(casCount > cas)
					llTemp += (casCount - cas) * log((casCount - cas)/(totCount - tot));
				if(conCount > con)
					llTemp += (conCount - con) * log((conCount - con)/(totCount - tot));
				ll[i] = llTemp;
			}
			else
				ll[i] = 1;
		}
		else { //Low cluster of cases
			if(lowCluster) {
				llTemp = con * log(con/tot);
				if(cas > 0)
					llTemp += cas * log(cas/tot);
				if(casCount > cas)
					llTemp += (casCount - cas) * log((casCount - cas)/(totCount - tot));
				if(conCount > con)
					llTemp += (conCount - con) * log((conCount - con)/(totCount - tot));
				ll[i] = llTemp;
			}
			else
				ll[i] = 1;
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
			if(ll[i * wCount + j] < 0) {
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
				if(ll[i * wCount + j] < 0) {
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
