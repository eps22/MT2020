ParallelMasterTracker <- function(compfolder, ncores=NULL, ParamsList){ 
  
  InitTime <- ParamsList$InitTime
  FinalTime <- ParamsList$FinalTime
  Resolution <- ParamsList$Resolution     
  TimestepDt <- ParamsList$TimestepDt
  LonDimName <- ParamsList$LonDimName
  LonVarName <- ParamsList$LonVarName    
  LatDimName <- ParamsList$LatDimName
  LatVarName <- ParamsList$LatVarName                      
  TimeDimName <- ParamsList$TimeDimName                   
  PressureVertLevelDimName <- ParamsList$PressureVertLevelDimName         
  SLPVarName <- ParamsList$SLPVarName                        
  U10VarName <- ParamsList$U10VarName                      
  V10VarName <- ParamsList$V10VarName                         
  ZVarName <- ParamsList$ZVarName 
  SmoothingPasses <- ParamsList$SmoothingPasses 
  SLPThreshold <- ParamsList$SLPThreshold
  ProductQuantileLowerLimit  <- ParamsList$ProductQuantileLowerLimit
  VorticityThreshold <- ParamsList$VorticityThreshold
  CalculateZeroVortRadiusThreshold <- ParamsList$CalculateZeroVortRadiusThreshold
  CalculateZeroVortRadiusNPoints <- ParamsList$CalculateZeroVortRadiusNPoints
  IfCheckZeroVortSymm <- ParamsList$IfCheckZeroVortSymm 
  ZeroVortRadiusMaxAllowedAsymm <- ParamsList$ZeroVortRadiusMaxAllowedAsymm
  ZeroVortRadiusMinSymmDirs <- ParamsList$ZeroVortRadiusMinSymmDirs
  ZeroVortRadiusUpperLimit <- ParamsList$ZeroVortRadiusUpperLimit
  ZeroVortRadiusLowerLimit <- ParamsList$ZeroVortRadiusLowerLimit
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
  source('./CalculateZeroVortRadius.R')
  source('./conditionalSLP.R')
  source('./FulfillsHart.R')
  source('./ThermalSymmB.R')
  source('./ThermalWind.R')
  
  suppressMessages(suppressWarnings(library(ncdf4)))
  suppressMessages(suppressWarnings(library(oce)))
  
  name <- strsplit(tail(strsplit(compfolder, '/')[[1]],n=1),'to')[[1]][1]
  
  allslp <- nc_open(paste(compfolder,'/output/outputfile-slp.nc',sep=''))
  
  triminittime <- strsplit(InitTime,'')[[1]]
  if(sum(grepl("[[:digit:]]", triminittime)) == length(triminittime)){InitTime <- as.numeric(InitTime)
  } else if(InitTime=='initial'){InitTime <- 1
  } else if(InitTime!='initial' & sum(grepl("[[:digit:]]", triminittime)) != length(triminittime)){
    stop('Specified initial timestep is not a valid one. Stopping execution, please check FindMedicanes.namelist.')
    InitTime <- 1
  }
  
  trimfinaltime <- strsplit(FinalTime,'')[[1]]
  if(sum(grepl("[[:digit:]]", trimfinaltime)) == length(trimfinaltime)){FinalTime <- as.numeric(FinalTime)
  } else if(FinalTime=='final'){FinalTime <- allslp$dim[[TimeDimName]]$len-2
  } else if(FinalTime!='final' & sum(grepl("[[:digit:]]", trimfinaltime)) != length(trimfinaltime)){
    stop('Specified final timestep is not a valid one. Stopping execution, please check FindMedicanes.namelist.')
    FinalTime <- allslp$dim[[TimeDimName]]$len-2
  }
  
  if((InitTime <=0) | (InitTime > allslp$dim[[TimeDimName]]$len-2)){
    stop('Specified initial timestep is beyond limits. Stopping execution, please check FindMedicanes.namelist.')
    InitTime <- 1
  }
  
  if((FinalTime <=0) | (FinalTime > allslp$dim[[TimeDimName]]$len-2)){
    stop('Specified final timestep is beyond limits. Stopping execution, please check FindMedicanes.namelist.')
    FinalTime <- allslp$dim[[TimeDimName]]$len-2
  }
  
  if((FinalTime < InitTime)){
    stop('Speficied final timestep is lower than initial timestep. Stopping execution, please check FindMedicanes.namelist.')
    InitTime <- 1
    FinalTime <- allslp$dim[[TimeDimName]]$len-2
  }
  
  
  allz <- nc_open(paste(compfolder,'/output/outputfile-z.nc',sep=''))
  
  if(Blowerpressurelevel<=Bupperpressurelevel){stop('Error: selected lower pressure level is not below the upper one and may lead to inconsistent medicanes detections. Stopping execution, please check FindMedicanes.namelist.')}
  if((Blowerpressurelevel %in% allz$dim[[PressureVertLevelDimName]]$vals)==FALSE){
    stop('Error: selected lower pressure level for the calculation of thermal symmetry parameter B is not included in the geopotential height pressure levels. Stopping execution, please check FindMedicanes.namelist.')
  }     
  if((Bupperpressurelevel %in% allz$dim[[PressureVertLevelDimName]]$vals)==FALSE){
    stop('Error: selected upper pressure level for the calculation of thermal symmetry parameter B is not included in the geopotential height pressure levels. Stopping execution, please check FindMedicanes.namelist.')
  }     
  
  if(LTWlowerpressurelevel<=LTWupperpressurelevel){
    stop('Selected lower pressure level for the calculation of lower tropospheric thermal wind calculation is above the upper one and may lead to inconsistent medicanes detections. Stopping execution, please check FindMedicanes.namelist.')
  }
  if(UTWlowerpressurelevel<=UTWupperpressurelevel){
    stop('Selected lower pressure level for the calculation of upper tropospheric thermal wind calculation is above the upper one and may lead to inconsistent medicanes detections. Stopping execution, please check FindMedicanes.namelist.')
  }
  if((LTWlowerpressurelevel %in% allz$dim[[PressureVertLevelDimName]]$vals)==FALSE){
    stop('Error: selected lower pressure level for the calculation of lower tropospheric thermal wind is not included in the geopotential height pressure levels. Stopping execution, please check FindMedicanes.namelist.')
  }    
  if((LTWupperpressurelevel %in% allz$dim[[PressureVertLevelDimName]]$vals)==FALSE){
    stop('Error: selected upper pressure level for the calculation of lower tropospheric thermal wind is not included in the geopotential height pressure levels. Stopping execution, please check FindMedicanes.namelist.')
  } 
  if((UTWlowerpressurelevel %in% allz$dim[[PressureVertLevelDimName]]$vals)==FALSE){
    stop('Error: selected lower pressure level for the calculation of upper tropospheric thermal wind is not included in the geopotential height pressure levels. Stopping execution, please check FindMedicanes.namelist.')
  }    
  if((UTWupperpressurelevel %in% allz$dim[[PressureVertLevelDimName]]$vals)==FALSE){
    stop('Error: selected upper pressure level for the calculation of upper tropospheric thermal wind is not included in the geopotential height pressure levels. Stopping execution, please check FindMedicanes.namelist.')
  } 
  

  if(!is.null(ncores) && ncores==0){
    stop('Error: 0 cores is not a valid number of cores to run the algorithm. Stopping execution.')
  } else if(!is.null(ncores) && ncores==1){
    cat("\n","Running medicanes tracking algorithm. Progress can be monitored in outfolder/findmedicanes.log.","\n")
    writeLines(c(""), paste(compfolder,'/findmedicanes.log',sep=''))
    
    medicanes <- list()
    for(currentimeidx in InitTime:FinalTime){
      sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
      cat(paste("Finding center on timestep ",currentimeidx," of ",(allslp$dim[[TimeDimName]]$len-2),"\n"))
      sink()
      medicanes[[currentimeidx]] <- ParallelFindCenters(compfolder, currentimeidx, TimestepDt, Resolution, SmoothingPasses, SLPThreshold, ProductQuantileLowerLimit, VorticityThreshold, CalculateZeroVortRadiusThreshold, CalculateZeroVortRadiusNPoints, IfCheckZeroVortSymm, ZeroVortRadiusMaxAllowedAsymm, ZeroVortRadiusMinSymmDirs, ZeroVortRadiusUpperLimit, ZeroVortRadiusLowerLimit, SLPminsClustersMinIBdistance, MaxNumberOfDifferentClusters, MinPointsNumberInCluster, IfCheckHartParamsConditions, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, PressureVertLevelDimName, LonDimName, LonVarName, LatDimName, LatVarName, SLPVarName, U10VarName, V10VarName, ZVarName)
    }
    saveRDS(medicanes, file=paste('../',name,'-track.RData',sep=''))
    sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
    cat(paste("Medicanes tracking algorithm completed successfully.","\n"))
    sink()
    cat("\n","Algorithm execution completed successfully","\n")
  } else{
    suppressMessages(suppressWarnings(library(foreach)))
    suppressMessages(suppressWarnings(library(doParallel)))
   
    saveeach <- 100
    exact <- floor((FinalTime-InitTime+1)/saveeach)
    rest <- (FinalTime-InitTime+1)%%saveeach
    if(rest==0){exact <- exact - 1}
 
    cat("\n","Running medicanes tracking algorithm. Progress can be monitored in outfolder/findmedicanes.log.","\n")

    writeLines(c(""), paste(compfolder,'/findmedicanes.log',sep=''))
    
    for(i in 1:(exact+1)){
      
    initial <- saveeach*(i-1)+InitTime
    sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
    cat(paste("Starting tracking for batch",i," of ",(exact+1),"\n"))
    sink()
    
    cores=detectCores()
    
    if(is.null(ncores)){
      cl <- makeCluster(cores[1]-1) 
    } else{
      if(ncores>cores[1]){
        ncores <- cores[1]-1
        cl <- makeCluster(ncores) 
      } else if(ncores==cores[1]){
        warning('Warning: At least one core should be reserved in order to prevent system overloading.')
        cl <- makeCluster(ncores) 
      } else{
        cl <- makeCluster(ncores) 
      }
    }
    registerDoParallel(cl)
    
    pckgs <- c('dd','ParallelFindCenters','isMedicane','TrackCenterCandidate','CalculateZeroVortRadius','conditionalSLP','FulfillsHart','ThermalSymmB','ThermalWind')
    if(exact==0){suminterv <- (FinalTime-InitTime)
    } else if(exact>0 & i<(exact+1)){suminterv <- (saveeach-1)
    } else if(exact>0 & i==(exact+1)){suminterv <- (rest-1)}
    medicanes <- foreach(currentimeidx=initial:(initial+suminterv), .packages = c( 'ncdf4', 'oce'), .export = pckgs) %dopar% {
      sink(paste(compfolder,'/findmedicanes.log',sep=''), append=TRUE)
      cat(paste("Finding center on timestep ",currentimeidx," of ",FinalTime,"\n"))
      sink()
      ParallelFindCenters(compfolder, currentimeidx, TimestepDt, Resolution, SmoothingPasses, SLPThreshold, ProductQuantileLowerLimit, VorticityThreshold, CalculateZeroVortRadiusThreshold, CalculateZeroVortRadiusNPoints, IfCheckZeroVortSymm, ZeroVortRadiusMaxAllowedAsymm, ZeroVortRadiusMinSymmDirs, ZeroVortRadiusUpperLimit, ZeroVortRadiusLowerLimit, SLPminsClustersMinIBdistance, MaxNumberOfDifferentClusters, MinPointsNumberInCluster, IfCheckHartParamsConditions, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, PressureVertLevelDimName, LonDimName, LonVarName, LatDimName, LatVarName, SLPVarName, U10VarName, V10VarName, ZVarName)
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
    cat("\n","Algorithm execution completed successfully","\n")
  }
  
 }
