ReadNamelist <- function(){
  
  fileslist <- list.files('TrackingAlgorithm')
  if(("FindMedicanes.namelist" %in% fileslist)==FALSE){stop('Namelist not found. Please check it is inside the TrackingAlgorithm folder.')}
  
  abc <- readLines(con = 'TrackingAlgorithm/FindMedicanes.namelist')
  
  paramslist <- list()
  paramslist$SmoothingPasses <- as.numeric(strsplit(abc[2],'=')[[1]][2])
  paramslist$ProductQuantileLowerLimit <- as.numeric(strsplit(abc[3],'=')[[1]][2])
  paramslist$CalculateRadiusThreshold <- trimws(strsplit(abc[4],'=')[[1]][2])
  paramslist$CalculateRadiusNPoints <- as.numeric(strsplit(abc[5],'=')[[1]][2])
  paramslist$IfCheckZGradSymm <- as.logical(trimws(strsplit(abc[6],'=')[[1]][2]))
  paramslist$ZGradSymmPressureLevel <- as.numeric(strsplit(abc[7],'=')[[1]][2])
  paramslist$ZGradSymmNumberOfDirections <- as.numeric(strsplit(abc[8],'=')[[1]][2])
  paramslist$ZGradSymmMeasure <- trimws(strsplit(abc[9],'=')[[1]][2])
  paramslist$ZGradSymmthreshold <- as.numeric(strsplit(abc[10],'=')[[1]][2])   
  paramslist$SLPminsClustersMinIBdistance <- as.numeric(strsplit(abc[11],'=')[[1]][2])
  paramslist$MaxNumberOfDifferentClusters <- as.numeric(strsplit(abc[12],'=')[[1]][2])
  paramslist$MinPointsNumberInCluster <- as.numeric(strsplit(abc[13],'=')[[1]][2])
  paramslist$IfCheckHartParamsConditions <- as.logical(trimws(strsplit(abc[14],'=')[[1]][2]))
  paramslist$HartConditionsTocheck <- as.numeric(strsplit(trimws(strsplit(abc[15],'=')[[1]][2]),',')[[1]])
  paramslist$Blowerpressurelevel <- as.numeric(strsplit(abc[16],'=')[[1]][2])
  paramslist$Bupperpressurelevel <- as.numeric(strsplit(abc[17],'=')[[1]][2])
  paramslist$Bmultiplemeasure <- trimws(strsplit(abc[18],'=')[[1]][2])
  paramslist$Bdirections <- as.numeric(strsplit(abc[19],'=')[[1]][2])
  paramslist$Bthreshold <- as.numeric(strsplit(abc[20],'=')[[1]][2])
  paramslist$LTWlowerpressurelevel <- as.numeric(strsplit(abc[21],'=')[[1]][2])
  paramslist$LTWupperpressurelevel <- as.numeric(strsplit(abc[22],'=')[[1]][2])
  paramslist$UTWlowerpressurelevel <- as.numeric(strsplit(abc[23],'=')[[1]][2])
  paramslist$UTWupperpressurelevel <- as.numeric(strsplit(abc[24],'=')[[1]][2])
  
  return(paramslist)

}

