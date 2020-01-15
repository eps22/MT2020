ParallelMasterTracker <- function(compfolder, Resolution, timestep, ncores=NULL, ParamsList){ 
  
  SmoothingPasses <- ParamsList$SmoothingPasses 
  ProductQuantileLowerLimit  <- ParamsList$ProductQuantileLowerLimit
  CalculateRadiusThreshold <- ParamsList$CalculateRadiusThreshold
  CalculateRadiusNPoints <- ParamsList$CalculateRadiusNPoints
  IfCheckZGradSymm <- ParamsList$IfCheckZGradSymm 
  ZGradSymmPressureLevel <- ParamsList$ZGradSymmPressureLevel      
  ZGradSymmNumberOfDirections <- ParamsList$ZGradSymmNumberOfDirections       
  ZGradSymmMeasure <- ParamsList$ZGradSymmMeasure               
  ZGradSymmthreshold <- ParamsList$ZGradSymmthreshold          
  SLPminsClustersMinIBdistance <- ParamsList$SLPminsClustersMinIBdistance
  MaxNumberOfDifferentClusters <- ParamsList$MaxNumberOfDifferentClusters
  MinPointsNumberInCluster <- ParamsList$MinPointsNumberInCluster 
  IfCheckHartParamsConditions <- ParamsList$IfCheckHartParamsConditions       
  HartConditionsTocheck <- ParamsList$HartConditionsTocheck          
  Blowerpressurelevel <- ParamsList$Blowerpressurelevel                
  Bupperpressurelevel <- ParamsList$Bupperpressurelevel              
  Bmultiplemeasure <- ParamsList$Bmultiplemeasure  
  Bdirections <- ParamsList$Bdirections                            
  Bthreshold <- ParamsList$Bthreshold                        
  LTWlowerpressurelevel <- ParamsList$LTWlowerpressurelevel         
  LTWupperpressurelevel <- ParamsList$LTWupperpressurelevel               
  UTWlowerpressurelevel <- ParamsList$UTWlowerpressurelevel              
  UTWupperpressurelevel <- ParamsList$UTWupperpressurelevel              
  
  
  source('./dd.R')
  source('./ParallelFindCenters.R')
  source('./isMedicane.R')
  source('./TrackCenterCandidate.R')
  source('./CalculateRadius.R')
  source('./CheckZGradSymm.R')
  source('./conditionalSLP.R')
  source('./FulfillsHart.R')
  source('./ThermalSymmB.R')
  source('./ThermalWind.R')
  
  library(ncdf4)
  library(oce)
  
  name <- strsplit(tail(strsplit(compfolder, '/')[[1]],n=1),'to')[[1]][1]
  
  allslp <- nc_open(paste(compfolder,'/output/outputfile-slp.nc',sep=''))
  if(!is.null(ncores)){
    if(ncores==1){
      writeLines(c(""), paste(compfolder,'/findmedicanes.log',sep=''))
      medicanes <- list()
      for(currentimeidx in 1:(allslp$dim$Time$len-2)){
        sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
        cat(paste("Finding center on timestep ",currentimeidx," of ",(allslp$dim$Time$len-2),"\n"))
        sink()
        medicanes[[currentimeidx]] <- ParallelFindCenters(compfolder, currentimeidx, ProductQuantileLowerLimit, SmoothingPasses, timestep, Resolution, SLPminsClustersMinIBdistance, MaxNumberOfDifferentClusters, MinPointsNumberInCluster, CalculateRadiusThreshold, CalculateRadiusNPoints, IfCheckZGradSymm, ZGradSymmPressureLevel, ZGradSymmNumberOfDirections, ZGradSymmMeasure, ZGradSymmthreshold, IfCheckHartParamsConditions, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel)
      }
      saveRDS(medicanes, file=paste('../',name,'-track.RData',sep=''))
      sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
      cat(paste("Medicanes finding algorithm completed successfully.","\n"))
      sink()
    }
  } else{
    library(foreach)
    library(doParallel)
    
    saveeach <- 100
    exact <- floor((allslp$dim$Time$len-2)/saveeach)
    rest <- (allslp$dim$Time$len-2)%%saveeach
 
    writeLines(c(""), paste(compfolder,'/findmedicanes.log',sep=''))
    if(IfCheckZGradSymm==TRUE){
      writeLines(c(""), paste(compfolder,'/checkzvortradius.log',sep=''))
    }
    for(i in 1:(exact+1)){
    initial <- saveeach*(i-1)+1
    sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
    cat(paste("Starting tracking for batch",i," of ",(exact+1),"\n"))
    sink()
    if(is.null(ncores)){
      cores=detectCores()
      cl <- makeCluster(cores[1]-1) 
    } else{
      cl <- makeCluster(ncores) 
    }
    registerDoParallel(cl)
    
    pckgs <- c('dd','ParallelFindCenters','isMedicane','TrackCenterCandidate','CalculateRadius','CheckZGradSymm','conditionalSLP','FulfillsHart','ThermalSymmB','ThermalWind')
    if(i<(exact+1)){suminterv <- (saveeach-1)} else if(i==(exact+1)){suminterv <- (rest-1)}
    medicanes <- foreach(currentimeidx=initial:(initial+suminterv), .packages = c( 'ncdf4', 'oce'), .export = pckgs) %dopar% {
      sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
      cat(paste("Finding center on timestep ",currentimeidx," of ",(allslp$dim$Time$len-2),"\n"))
      sink()
      ParallelFindCenters(compfolder, currentimeidx, ProductQuantileLowerLimit, SmoothingPasses, timestep, Resolution, SLPminsClustersMinIBdistance, MaxNumberOfDifferentClusters, MinPointsNumberInCluster, CalculateRadiusThreshold, CalculateRadiusNPoints, IfCheckZGradSymm, ZGradSymmPressureLevel, ZGradSymmNumberOfDirections, ZGradSymmMeasure, ZGradSymmthreshold, IfCheckHartParamsConditions, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel)
    }
    
    sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
    cat(paste("Batch ",i," calculation completed successfully. Writing result to file ...","\n"))
    sink() 
    stopCluster(cl)
    
    if(i==1){currentlist <- medicanes} else{
      previouslist <- readRDS(paste('../',name,'-track.RData',sep=''))
      currentlist <- append(previouslist,medicanes)
    }
    saveRDS(currentlist, file=paste('../',name,'-track.RData',sep=''))
    }
    sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
    cat(paste("Medicanes tracking algorithm completed successfully.","\n"))
    sink()
  }
  
 }
