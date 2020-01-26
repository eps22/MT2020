ReadNamelist <- function(){
  
  fileslist <- list.files('TrackingAlgorithm')
  if(("FindMedicanes.namelist" %in% fileslist)==FALSE){stop('Namelist not found. Please check it is inside the TrackingAlgorithm folder.')}
  
  abc <- readLines(con = 'TrackingAlgorithm/FindMedicanes.namelist')
  
  paramslist <- list()
  paramslist$SmoothingPasses <- as.numeric(strsplit(abc[2],'=')[[1]][2])
  paramslist$SLPThreshold <- as.numeric(strsplit(abc[3],'=')[[1]][2])
  paramslist$ProductQuantileLowerLimit <- as.numeric(strsplit(abc[4],'=')[[1]][2])
  paramslist$VorticityThreshold <- as.numeric(strsplit(abc[5],'=')[[1]][2])
  paramslist$CalculateZeroVortRadiusThreshold <- trimws(strsplit(abc[6],'=')[[1]][2])
  paramslist$CalculateZeroVortRadiusNPoints <- as.numeric(strsplit(abc[7],'=')[[1]][2])
  paramslist$IfCheckZeroVortSymm <- as.logical(trimws(strsplit(abc[8],'=')[[1]][2]))
  paramslist$ZeroVortRadiusMaxAllowedAsymm <- as.numeric(strsplit(abc[9],'=')[[1]][2])
  paramslist$ZeroVortRadiusMinSymmDirs <- as.numeric(strsplit(abc[10],'=')[[1]][2])
  paramslist$ZeroVortRadiusUpperLimit <- as.numeric(strsplit(abc[11],'=')[[1]][2])
  paramslist$ZeroVortRadiusLowerLimit <- as.numeric(strsplit(abc[12],'=')[[1]][2])
  paramslist$SLPminsClustersMinIBdistance <- as.numeric(strsplit(abc[13],'=')[[1]][2])
  paramslist$MaxNumberOfDifferentClusters <- as.numeric(strsplit(abc[14],'=')[[1]][2])
  paramslist$MinPointsNumberInCluster <- as.numeric(strsplit(abc[15],'=')[[1]][2])
  paramslist$IfCheckHartParamsConditions <- as.logical(trimws(strsplit(abc[16],'=')[[1]][2]))
  paramslist$HartConditionsTocheck <- as.numeric(strsplit(trimws(strsplit(abc[17],'=')[[1]][2]),',')[[1]])
  paramslist$Blowerpressurelevel <- as.numeric(strsplit(abc[18],'=')[[1]][2])
  paramslist$Bupperpressurelevel <- as.numeric(strsplit(abc[19],'=')[[1]][2])
  paramslist$Bmultiplemeasure <- trimws(strsplit(abc[20],'=')[[1]][2])
  paramslist$Bdirections <- as.numeric(strsplit(abc[21],'=')[[1]][2])
  paramslist$Bthreshold <- as.numeric(strsplit(abc[22],'=')[[1]][2])
  paramslist$LTWlowerpressurelevel <- as.numeric(strsplit(abc[23],'=')[[1]][2])
  paramslist$LTWupperpressurelevel <- as.numeric(strsplit(abc[24],'=')[[1]][2])
  paramslist$UTWlowerpressurelevel <- as.numeric(strsplit(abc[25],'=')[[1]][2])
  paramslist$UTWupperpressurelevel <- as.numeric(strsplit(abc[26],'=')[[1]][2])
  
  return(paramslist)

}

