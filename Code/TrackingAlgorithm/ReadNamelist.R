ReadNamelist <- function(){
  
  fileslist <- list.files('TrackingAlgorithm')
  if(("FindMedicanes.namelist" %in% fileslist)==FALSE){stop('Namelist not found. Please check it is inside the TrackingAlgorithm folder.')}
  
  abc <- readLines(con = 'TrackingAlgorithm/FindMedicanes.namelist')
  abc <- abc[unlist(sapply(trimws(abc),function(x){nchar(x)>1},USE.NAMES = FALSE))]
  
  paramslist <- list()
  paramslist$InitTime <- trimws(strsplit(abc[1],'=')[[1]][2])             
  paramslist$FinalTime <- trimws(strsplit(abc[2],'=')[[1]][2])
  paramslist$Resolution <- as.numeric(strsplit(abc[3],'=')[[1]][2])           
  paramslist$TimestepDt <- as.numeric(strsplit(abc[4],'=')[[1]][2])
  paramslist$LonDimName <- trimws(strsplit(abc[5],'=')[[1]][2])           
  paramslist$LonVarName <- trimws(strsplit(abc[6],'=')[[1]][2])   
  paramslist$LatDimName <- trimws(strsplit(abc[7],'=')[[1]][2])   
  paramslist$LatVarName <- trimws(strsplit(abc[8],'=')[[1]][2])                     
  paramslist$TimeDimName <- trimws(strsplit(abc[9],'=')[[1]][2])                       
  paramslist$PressureVertLevelDimName <- trimws(strsplit(abc[10],'=')[[1]][2])           
  paramslist$SLPVarName <- trimws(strsplit(abc[11],'=')[[1]][2])                        
  paramslist$U10VarName <- trimws(strsplit(abc[12],'=')[[1]][2])                        
  paramslist$V10VarName <- trimws(strsplit(abc[13],'=')[[1]][2])                         
  paramslist$ZVarName <- trimws(strsplit(abc[14],'=')[[1]][2])                           
  paramslist$SmoothingPasses <- as.numeric(strsplit(abc[15],'=')[[1]][2])
  paramslist$SLPThreshold <- as.numeric(strsplit(abc[16],'=')[[1]][2])
  paramslist$ProductQuantileLowerLimit <- as.numeric(strsplit(abc[17],'=')[[1]][2])
  paramslist$VorticityThreshold <- as.numeric(strsplit(abc[18],'=')[[1]][2])
  paramslist$CalculateZeroVortRadiusThreshold <- trimws(strsplit(abc[19],'=')[[1]][2])
  paramslist$CalculateZeroVortRadiusNPoints <- as.numeric(strsplit(abc[20],'=')[[1]][2])
  paramslist$IfCheckZeroVortSymm <- as.logical(trimws(strsplit(abc[21],'=')[[1]][2]))
  paramslist$ZeroVortRadiusMaxAllowedAsymm <- as.numeric(strsplit(abc[22],'=')[[1]][2])
  paramslist$ZeroVortRadiusMinSymmDirs <- as.numeric(strsplit(abc[23],'=')[[1]][2])
  paramslist$ZeroVortRadiusUpperLimit <- as.numeric(strsplit(abc[24],'=')[[1]][2])
  paramslist$ZeroVortRadiusLowerLimit <- as.numeric(strsplit(abc[25],'=')[[1]][2])
  paramslist$SLPminsClustersMinIBdistance <- as.numeric(strsplit(abc[26],'=')[[1]][2])
  paramslist$MaxNumberOfDifferentClusters <- as.numeric(strsplit(abc[27],'=')[[1]][2])
  paramslist$MinPointsNumberInCluster <- as.numeric(strsplit(abc[28],'=')[[1]][2])
  paramslist$IfCheckHartParamsConditions <- as.logical(trimws(strsplit(abc[29],'=')[[1]][2]))
  paramslist$HartConditionsTocheck <- as.numeric(strsplit(trimws(strsplit(abc[30],'=')[[1]][2]),',')[[1]])
  paramslist$Blowerpressurelevel <- as.numeric(strsplit(abc[31],'=')[[1]][2])
  paramslist$Bupperpressurelevel <- as.numeric(strsplit(abc[32],'=')[[1]][2])
  paramslist$Bmultiplemeasure <- trimws(strsplit(abc[33],'=')[[1]][2])
  paramslist$Bdirections <- as.numeric(strsplit(abc[34],'=')[[1]][2])
  paramslist$Bthreshold <- as.numeric(strsplit(abc[35],'=')[[1]][2])
  paramslist$LTWlowerpressurelevel <- as.numeric(strsplit(abc[36],'=')[[1]][2])
  paramslist$LTWupperpressurelevel <- as.numeric(strsplit(abc[37],'=')[[1]][2])
  paramslist$UTWlowerpressurelevel <- as.numeric(strsplit(abc[38],'=')[[1]][2])
  paramslist$UTWupperpressurelevel <- as.numeric(strsplit(abc[39],'=')[[1]][2])
  
  return(paramslist)

}

