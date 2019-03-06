# Implementation of DDC by R Hyde 2017
# Currently does not include cluster axes fitting to data
rm(list=ls())
library(ggplot2)
library(fields)
library(scatterplot3d)

##### To change data sets and use suitable radii comment /  uncomment the following lines #####
# AllData = array(read.csv("Gaussian5000.csv", skip=1))
# InitR = c(0.1,0.05) # test multiple radii
# InitR = array(InitR, dim=c(1,length(InitR)))
AllData = array(read.csv("ChainLink3DNoise.csv", skip=1))
InitR = array(c(0.2)) # test with single radius

Verbose = 1 # not implemented
Merge = 0 # not implemented

idx = dim(AllData)
Data = data.frame(AllData[,1:idx[2]-1])

  # Initialize
  if (length(InitR) == 1) { # use same radii across each dimension if only one provided
    
    if (Verbose == 1) {
      print('Using equal radii')
      
    }
    
    InitR = array(InitR, dim=c(1,ncol(Data)))
  }  else if ( length(InitR)<ncol(Data) ) {
    stop('Radii not single or equal to data dimensions')
  # exit function here %%%%%
  }

Results = data.frame( Cluster = array(0, nrow(Data)), Centre=array(0, nrow(Data)), Radii=array(0, c(nrow(Data), length(InitR)) ) )


  NumClusters = 0; # initial number of clusters

  while (sum(Results["Cluster"] == 0)>0){
    print(sum(Results["Cluster"] == 0))
    NumClusters = NumClusters+1;
    # find initial centre
    Unassigned = which(Results["Cluster"] == 0)
    CentreIdx = attr( which.min( rowSums( scale(Data[Unassigned,], scale=FALSE)^2 ) ), "names") # find data closest to mean
    Results[CentreIdx,"Centre"] = NumClusters # assign it as temp centre of cluster
    Results[CentreIdx,"Cluster"] = NumClusters # assign it to cluster
    Results[CentreIdx, c(3:(2+length(InitR)) )] = InitR
    # include data
    Included = which ( rdist(Data[CentreIdx,],Data[Unassigned,]) < sqrt(rowSums(InitR^2)) )
    Included = Unassigned[Included]
    Results[Included,"Cluster"] = NumClusters
    
    ClusterDists = rdist(Data[CentreIdx,],Data[Included,])
    ClMean = mean(ClusterDists) # mean distances
    ClStdDev = sd(ClusterDists) #std dev distances
    Results[Included,"Cluster"] = 0 # remove cluster assign
    Included = Included[ClusterDists < (ClMean + 3*ClStdDev)] # re-assign keeping only those within 3 sigma
    if (length(Included)>1){
      Results[Included, "Cluster"] = NumClusters
    
      # move centre to final position (mean of temporary clustered data)
      CentreIdx = attr( which.min( rowSums( scale(Data[Included,], scale=FALSE)^2 ) ), "name") # find data closest to cluster mean
      # remove temporary cluster assignments and temp centre
      Results[which(Results["Centre"]==NumClusters),"Centre"]=0 # remove temp centre
      Results[which(Results["Cluster"]==NumClusters),"Centre"]=0 # remove temporary cluster assignment of the centre
      Results[which(Results["Cluster"]==NumClusters),"Cluster"]=0 # de-assign all data
      # re-assign to the new centre
      Results[CentreIdx,"Centre"] = NumClusters # assign it as temp centre of cluster
      Results[CentreIdx, c(3:(2+length(InitR)) )] = InitR # set cluster radii
      Results[CentreIdx,"Cluster"] = NumClusters # assign it to cluster
      # find new data
      Included = which ( rdist(Data[CentreIdx,],Data[Unassigned,]) < sqrt(rowSums(InitR^2)) ) # find any data filling into new cluster position
      Included = Unassigned[Included]
      Results[Included,"Cluster"] = NumClusters # assign final data
      ClusterDists = rdist(Data[CentreIdx,],Data[Included,]) # distances from final cluster centre to cluster data
      ClMean = mean(ClusterDists) # mean distances
      ClStdDev = sd(ClusterDists) #std dev distances
      Results[Included,"Cluster"] = 0 # remove cluster assign
      Included = Included[ClusterDists < (ClMean + 3*ClStdDev)] # re-assign keeping only those within 3 sigma
      if (length(Included)>1){
        Results[Included,"Cluster"] = NumClusters # assign final data
      } else {
        Results[CentreIdx,"Centre"] = NumClusters # assign it as temp centre of cluster
        Results[CentreIdx,"Cluster"] = NumClusters # assign it to cluster
        Results[CentreIdx, c(3:(2+length(InitR)) )] = InitR
      }
    } else {
      Results[CentreIdx,"Centre"] = NumClusters # assign it as temp centre of cluster
      Results[CentreIdx,"Cluster"] = NumClusters # assign it to cluster
      Results[CentreIdx, c(3:(2+ncol(InitR)) )] = InitR
    }
    
    # update radii to match data
    Results[CentreIdx, c(3:(2+length(InitR)) )] =
      apply( abs( sweep(Data[Included,], 2, as.numeric(Data[CentreIdx,]),"-") ), 2, max)
    # re-assign data within the new radii, under 3-sigma deviation
    Included = which ( rdist(Data[CentreIdx,],Data[Unassigned,]) < sqrt(rowSums(InitR^2)) ) # find any data filling into new cluster position
    Included = Unassigned[Included]
    Results[Included,"Cluster"] = NumClusters # assign final data
    ClusterDists = rdist(Data[CentreIdx,],Data[Included,]) # distances from final cluster centre to cluster data
    ClMean = mean(ClusterDists) # mean distances
    ClStdDev = sd(ClusterDists) #std dev distances
    Results[Included,"Cluster"] = 0 # remove cluster assign
    Included = Included[ClusterDists < (ClMean + 3*ClStdDev)] # re-assign keeping only those within 3 sigma
    if (length(Included)>1){
      Results[Included,"Cluster"] = NumClusters # assign final data
    } else {
      Results[CentreIdx,"Centre"] = NumClusters # assign it as temp centre of cluster
      Results[CentreIdx,"Cluster"] = NumClusters # assign it to cluster
      Results[CentreIdx, c(3:(2+length(InitR)) )] = InitR
    }
  
  }
  
  if (length(InitR)==2){
    ColPal = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    ggplot(data=Data, aes(x=X1, y=X2, color=Results["Cluster"])
         ) + geom_point(size=1) + scale_color_gradientn(colours = sample(ColPal, max(Results$Cluster)))
  } else if (length(InitR)==3){
    ColPal = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    PC = sample(ColPal, max(Results$Cluster)) # set plot colors
    scatterplot3d(Data[,1], Data[,2], Data[,3], color=PC[Results$Cluster])
    
  }
