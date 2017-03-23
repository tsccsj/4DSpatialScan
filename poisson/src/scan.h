#ifndef SCANH
#define SCANH

void getPECount(double * x1, double * y1, double * x2, double * y2, int * nPop, int * nEvent, int locCount, double wSize, int wCount, int * popInW, int * eventInW, double elimIntersectOD);
void loglikelihood(double * ll, int * popInW, int * eventInW, int totalWindow, int popCount, int eventCount, int highLow);
void findTopNCluster(double * x1, double * y1, double * x2, double * y2, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters);

#endif
