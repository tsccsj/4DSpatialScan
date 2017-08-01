#ifndef MCH
#define MCH

int * monteCarlo(double * x1, double * y1, double * x2, double * y2, int * nPop, int * popInW, int locCount, int popCount, int eventCount, double wSize, int wCount, int highLow, double elimIntersectOD, double * clusterLL, int nClusters, int nSim);
int * monteCarloOld(double * x1, double * y1, double * x2, double * y2, int * nPop, int locCount, int popCount, int eventCount, int * clusterEvent, int * center, double * cRadius, bool * highCluster, int nClusters, int nSim);

#endif
