ParallelFindCenters <- function(compfolder, currentimeidx, TimestepDt, Resolution, SmoothingPasses, SLPThreshold, ProductQuantileLowerLimit, VorticityThreshold, CalculateZeroVortRadiusThreshold, CalculateZeroVortRadiusNPoints,  IfCheckZeroVortSymm, ZeroVortRadiusMaxAllowedAsymm, ZeroVortRadiusMinSymmDirs, ZeroVortRadiusUpperLimit, ZeroVortRadiusLowerLimit, SLPminsClustersMinIBdistance, MaxNumberOfDifferentClusters, MinPointsNumberInCluster, IfCheckHartParamsConditions, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, PressureVertLevelDimName, LonDimName, LonVarName, LatDimName, LatVarName, SLPVarName, U10VarName, V10VarName, ZVarName){

  
  isthesame <- function(a,b){(dd(a,b)*Resolution)<SLPminsClustersMinIBdistance}
  reducelist <- function(list){sapply(list,function(y){sapply(list,function(x){isthesame(y,x)},USE.NAMES = FALSE)})}
  
  #Get slp slice
  slpT <-  ncvar_get(nc_open(paste(compfolder,'/output/outputfile-slp.nc',sep='')),SLPVarName,start=c(1,1,currentimeidx),count=c(-1,-1,1))

  gradT <- suppressMessages(oce::grad(slpT, x=seq(1,dim(slpT)[1], by=1)*Resolution, y=seq(1,dim(slpT)[2],by=1)*Resolution))
  gradTx <- gradT$gx 
  gradTy <- gradT$gy
  gradTxx <- (suppressMessages(oce::grad(gradTx, x=seq(1,dim(slpT)[1], by=1)*Resolution, y=seq(1,dim(slpT)[2],by=1)*Resolution)))$gx
  gradTyy <- (suppressMessages(oce::grad(gradTy, x=seq(1,dim(slpT)[1], by=1)*Resolution, y=seq(1,dim(slpT)[2],by=1)*Resolution)))$gy
  laplacianT <- gradTxx + gradTyy
  
  #Get 10m wind components slices
  uv10x <- ncvar_get(nc_open(paste(compfolder,'/output/outputfile-uvmet10-U.nc',sep='')),U10VarName,start=c(1,1,currentimeidx),count=c(-1,-1,1))
  uv10y <- ncvar_get(nc_open(paste(compfolder,'/output/outputfile-uvmet10-V.nc',sep='')),V10VarName,start=c(1,1,currentimeidx),count=c(-1,-1,1))
  uv10xx <- nc_open(paste(compfolder,'/output/outputfile-uvmet10-U.nc',sep=''))
  
  if(uv10xx$dim[[LonDimName]]$vals[2]>uv10xx$dim[[LonDimName]]$vals[1]){revx <- FALSE} else{revx <- TRUE}
  if(uv10xx$dim[[LatDimName]]$vals[2]>uv10xx$dim[[LatDimName]]$vals[1]){revy <- FALSE} else{revy <- TRUE}
  xseq <- seq(1,dim(slpT)[1], by=1)
  if(revx==TRUE){xseq <- rev(xseq)}
  yseq <- seq(1,dim(slpT)[2],by=1)
  if(revy==TRUE){yseq <- rev(yseq)}
  uv10curl <- curl(uv10x, uv10y, x=xseq*Resolution, 
                   y=yseq*Resolution, geographical = FALSE, method = 1)$curl
  
  #Compute the product of SLP laplacian and the curl of the 10m wind speed components
  product <- laplacianT*uv10curl
  product <- suppressMessages(oce::matrixSmooth(product, passes = SmoothingPasses)) 
  #Remove values on the borders to prevent artifacts to appear in the percentile list (i.e., allow blending area)
  for(a in 1:nrow(product)){
    for(b in 1:ncol(product)){
      if(a<10 | a>((dim(slpT)[1])-10) | b<10 | b>((dim(slpT)[2])-10)){product[a,b] <- 0}
    }
  }
  
  product[slpT>SLPThreshold] <- 0 #First filter: set product to zero if SLP is above a given threshold
  
  #Select all the points above a certain quantile as cyclone center candidates 
  LAPsorted <- sort(c(product), decreasing = TRUE)
  productthreshold <- as.numeric(quantile(product,probs = c(ProductQuantileLowerLimit)))
  LAPsortedp <- LAPsorted[LAPsorted>productthreshold]
  
  if(length(LAPsortedp)>0){
    centercandidates <- sapply(LAPsortedp,function(x){which(product==x, arr.ind=TRUE)})
    cclist <- split(centercandidates, rep(1:ncol(centercandidates), each = nrow(centercandidates)))
  } else{cclist <- list()}
  
  #Second filter: select the points with vorticity above a certain threshold
  vorts <- lapply(cclist,function(x){uv10curl[x[1],x[2]]})
  cclist <- cclist[vorts>VorticityThreshold]
  
  #Check whether found center candidates fulfill vorticity closed circulation (vorticity symmetry)
  if(IfCheckZeroVortSymm==TRUE){
    if(length(cclist)>0){
      idxstokeep <- c()
      for(i in 1:length(cclist)){
        #track <- TrackCenterCandidate(compfolder,cclist[[i]],currentimeidx,Resolution,timestep)
        zerovortradius <- CalculateZeroVortRadius(Resolution, uv10curl, cclist[[i]], threshold=CalculateZeroVortRadiusThreshold, npoints=CalculateZeroVortRadiusNPoints, maxasymm = ZeroVortRadiusMaxAllowedAsymm, minsymmdirs = ZeroVortRadiusMinSymmDirs)
        print(paste("Timestep:",currentimeidx,". Center:",as.numeric(cclist[[i]][1]),as.numeric(cclist[[i]][2]),". ZeroVortRadius:",as.numeric(zerovortradius),"km","\n"))
	      if(zerovortradius<ZeroVortRadiusUpperLimit & zerovortradius > ZeroVortRadiusLowerLimit){idxstokeep <- c(idxstokeep, i)}
      }
      cclist <- cclist[idxstokeep]
    }
  }
  cclistcp <- cclist
  
  #Find different centers in a given radius
  if(length(cclist)>1){
    false <- c(0)
    keep <- c()
    while(length(false)>0 && length(keep)<MaxNumberOfDifferentClusters){
      tfvecrec <- reducelist(cclist)
      true <- min(as.numeric(names(which(tfvecrec[1,]==TRUE))))
      if((true[1] %in% keep)!=TRUE){keep <- c(keep,true[1])}
      false <- as.numeric(names(which(tfvecrec[1,]==FALSE)))
      if(length(false)>1){
        falselist <- split(centercandidates[,false], rep(false, each = nrow(centercandidates[,false])))
        cclist <- falselist
      }
      if(length(false)==1){
        if((false %in% keep)!=TRUE){keep <- c(keep,false)}
        false <- c()
      }
    }

  #Compute the final list of independent found centers (one per group)
    centercandidateslist <- split(centercandidates[,keep],rep(keep,each=2))
  } else if(length(cclist)==1){
    centercandidateslist <- cclist
  } else{
    centercandidateslist <- list()
  }
  
    #Now create groups, as many as number of different centers.
  #Find the center candidates closer than X km to each different center.
  
  #Track each center two timesteps forward, and check if it fulfills the conditions to be medicane
  #-the Hart conditions-. If so, track it; if not, remove it.

  centersdf <- c()
  if(length(centercandidateslist)>0){
    for(centeridx in 1:length(centercandidateslist)){
      candidategroup <- cclistcp[as.numeric(which(sapply(cclistcp,function(x){dd(x,centercandidateslist[[centeridx]])*Resolution})<SLPminsClustersMinIBdistance))]
      print(centeridx)
      #print(length(candidategroup))
      sortedbySLPcandidategroup <- candidategroup[names(sort(unlist(lapply(candidategroup, function(x){slpT[x[1],x[2]]}))))]
      candidategroup <- sortedbySLPcandidategroup
      
      if(length(candidategroup)>MinPointsNumberInCluster){
        if(IfCheckHartParamsConditions==TRUE){
          ll <- 0
          im <- list()
          im$logical <- FALSE
          while(im$logical==FALSE && ll<length(candidategroup)){
            ll <- ll + 1
            center <- candidategroup[[ll]]
            closerthanDDminclusteridx <- min(which((as.numeric(lapply(centercandidateslist, function(x){dd(x,center)*Resolution}))<SLPminsClustersMinIBdistance)==TRUE))
            if(closerthanDDminclusteridx==centeridx){
	            zerovortradius2 <- CalculateZeroVortRadius(Resolution, uv10curl, center, threshold=CalculateZeroVortRadiusThreshold, npoints=CalculateZeroVortRadiusNPoints, maxasymm = ZeroVortRadiusMaxAllowedAsymm, minsymmdirs = ZeroVortRadiusMinSymmDirs)
              im <- isMedicane(compfolder, center, currentimeidx, timestep = TimestepDt, resolution = Resolution, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, PressureVertLevelDimName, SLPVarName, ZVarName)
              print(ll)
              print(center)
              print(im)
              if(im$logical==TRUE){
                if(length(centersdf)==0){
                  centersdf <- rbind(centersdf, center)
                } else{
                  if( (min(dd(centersdf,center))*Resolution) > SLPminsClustersMinIBdistance ){
                    centersdf <- rbind(centersdf, center)
                  }
                }
              }
            }
          }
        } 
        if(IfCheckHartParamsConditions==FALSE){ 
           center <- candidategroup[[1]]
           if(length(centersdf)==0){
             centersdf <- rbind(centersdf, center)
           } else{
             if( (min(dd(centersdf,center))*Resolution) > SLPminsClustersMinIBdistance ){
               centersdf <- rbind(centersdf, center)
             }
           }
        }
      }

    }
  }
  
  if(length(centersdf)==0){
    toreturn <- c(NA,NA)
  } else if(length(centersdf)>=1){
    toreturn <- unique(centersdf)
  }

  return(toreturn)
  
}
