#ifndef MCH
#define MCH

int * monteCarloOld(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim);
int * monteCarlo(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, double wSize, int wCount, double elimIntersectOD, bool highLow, double * clusterLL, bool * highCluster, int nClusters, int nSim);

#endif
