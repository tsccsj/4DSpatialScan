# 4DSpatialScan
An approach to find clusters (local excesses of events) in origin-destination movement datasets and more generally spatial interactions involving two locations. 
It is an extension of Martin Kulldorff's Spatial Scan Statistics(1997) to find the spatial clusters of spatial interactions (e.g. OD movement).
Two models are implemented: a *Poisson Model* and a *Bernoulli Model*.

## Publication
Gao, Y., Ting, L., Wang, S., Jeong, MH and Soltani, K., 2018. A multidimensional spatial scan statistics approach t movement pattern comparison. International Journal of Geographical Information Science, doi: 10.1080/13658816.2018.1426859.
https://www.tandfonline.com/doi/full/10.1080/13658816.2018.1426859

## Bernoulli model
A Bernoulli model handles spatial interactions that are in either one of two states (e.g. the migration of young people vs. the migration of senior people), which is often used to compare the spatial distributions of two types of spatial interactions, such as a case-control study.
### To run
4DScanBer InputFile WindowInc WindowCount NumberofClusters NumberofMonteCarlo HighOrLowIndicator ElimIntersectOD
 1. InputFile: a CSV without header that has 4 column: xorigin, yorigin, xdestination, ydestination, numberOfCases and numberOfControls for each location
 2. WindowInc: the increment of scan windows
 3. WindowCount: the number of scan windows put at each location, and thus the maximum scan window is (WindowInc * WindowCount)
 4. NumberofClusters: number of clusters to detect
 5. NumberofMonteCarlo: number of Monte Carlo simulation
 6. HighOrLowIndicator:
	* 1:  high value clusters only
	* -1: low value clusters only
	* 0:  both
 7. ElimIntersectOD: Whether clusters with intersecting origin and destination are allowed:
	* -1:  Allow intersect (no limit)
	* 1.0: Center of origin (destination) not in the area of destination (origin)
	* 2.0: No intersection allowed

## Poisson model
A Poisson model deals with the number of spatial interactions occurring in a time interval and a pair of OD spatial regions.
### To run
4DScanPoi InputFile WindowInc WindowCount NumberofClusters NumberofMonteCarlo HighOrLowIndicator ElimIntersectOD
 1. InputFile: a CSV without header that has 6 column: xorigin, yorigin, xdestination, ydestination, numberOfCases and intensity for each location
 2. WindowInc: the increment of scan windows
 3. WindowCount: the number of scan windows put at each location, and thus the maximum scan window is (WindowInc * WindowCount)
 4. NumberofClusters: number of clusters to detect
 5. NumberofMonteCarlo: number of Monte Carlo simulation
 6. HighOrLowIndicator:
	* 1: high value clusters only
	* -1: low value clusters only
	* 0: both
 7. ElimIntersectOD: Whether clusters with intersecting origin and destination are allowed:
	* -1:  Allow intersect (no limit)
	* 1.0: Center of origin (destination) not in the area of destination (origin)
	* 2.0: No intersection allowed

## Reference
Kulldorff, M., 1997. A spatial scan statistic. Communications in Statistics-Theory and methods, 26(6), pp.1481-1496.
 
