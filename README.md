# DDC-R
DDC Algorithm Coded in R

Data Density based Clustering

Hyde, R.; Angelov, P., "Data density based clustering," Computational Intelligence (UKCI), 2014 14th UK Workshop on , vol., no., pp.1,7, 8-10 Sept. 2014

doi: 10.1109/UKCI.2014.6930157

Downloadable from: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6930157&isnumber=6930143

**Please note, this was a quick and dirty port from the original Matlab code and is not intended to be complete or final.

Files:

DDC.Rproj - Project file

DDC.R - DDC algorithm implementation including the basic algorithm, together with plot functionality for 2D and 3D results. This file should be separated to different function files. The DDC algorithm is contained in lines 40-114. Note that some calculations are repeated to demonstrate clarity with the journal paper. The code could be simplified by separation of code to functions and reapeatedly calling them.

gaussian5000.csv - data file containg 5 random clusters of gaussian distributed data of 5000+ samples for cluster. This is the type of data distributions best suited to DDC.

ChainLink3DNoise.csv - standard test dataset of 3D 'chain-link' type data groups. Although DDC will created correct, pure clusters this data shows the limitations of DDC in working with non-standard data distrbutions. (To work with this type of data, see DDCAS, Data Density based Clustering for Arbitrary Shapes.)
