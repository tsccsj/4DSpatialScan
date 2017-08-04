#ifndef MCH
#define MCH

int * monteCarlo(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, double wSize, int wCount, double elimIntersectOD, int highLow, double * clusterLL, int nClusters, int nSim);
int * monteCarloOld(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim);

#endif
