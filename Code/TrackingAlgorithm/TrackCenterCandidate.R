TrackCenterCandidate <- function(compfolder,centertotrack,currentimeidx,resolution,timestep,SLPVarName){
  
  #compfolder: string with path to input folder
  #currentimeidx: timestamp index
  #resolution: nc grid resolution in km
  #timestep: time interval between consecutive timestamps in hours 
  
  
  
  slpnetcdf <- paste(compfolder,'/output/outputfile-slp.nc',sep='')
  
  library(ncdf4)  
  library(oce)
  
  ncinslp <- nc_open(slpnetcdf)
  slp <- ncvar_get(ncinslp,SLPVarName,start = c(1,1,currentimeidx),count=c(-1,-1,3))
  
  centers <- c()
  slps <- c()
  previnstspeeds <- c()
  
  for(timeidx in 1:3){
    slpt <- slp[,,timeidx]
    slpSM <- suppressMessages(oce::matrixSmooth(slpt))
    slpSMsorted <- sort(c(slpSM), decreasing = FALSE)
    if(timeidx==1 | timeidx==2){
      center <- centertotrack
      previnstspeed <- 0
    } else{
      centeredSLPmatrixidxs <- which((dd(slpSM,tail(centers,n=1))*resolution)<500,arr.ind=TRUE)
      centeredSLPmatrixvalues <- slpSM[centeredSLPmatrixidxs]
      sortedvaluesidxs <- sort(centeredSLPmatrixvalues, index.return=TRUE)$ix
      sortedlistofcenterstotry <- centeredSLPmatrixidxs[sortedvaluesidxs,,drop=FALSE]
      i <- 1
      center <- sortedlistofcenterstotry[i,,drop=FALSE]
      while(conditionalSLP(slpSM,tail(centers,n=2)[1,,drop=FALSE],tail(centers,n=1),center,tail(previnstspeeds,n=1),timestep,resolution) == FALSE){    
        i <- i+1
        center <- sortedlistofcenterstotry[i,,drop=FALSE]
      }
      previnstspeed <- (sqrt( ( (tail(centers,n=1))[1] - center[1] )^2 + ( (tail(centers,n=1))[2] - center[2] )^2) )*resolution/timestep
    }
    previnstspeeds <- c(previnstspeeds, previnstspeed)
    centers <- rbind(centers, center)
    slps <- c(slps, slpSM[center[1],center[2]])
  }  
 
  slps <- slps[c(1,3)]
  previnstspeeds <- previnstspeeds[c(1,3)]
  centers <- centers[c(1,3),]
  
  track <- list()
  track$slp <- slps
  track$centerspeed <- previnstspeeds
  track$centers <- centers
  return(track)
  
} 







