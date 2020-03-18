FulfillsHart <- function(compfolder, currentimeidx, consecutivecenters, resolution, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, PressureVertLevelDimName, ZVarName){
  
  conditions <- c()
  values <- c()
  
  if((1 %in% HartConditionsTocheck)==TRUE){
    Bparam <- ThermalSymmB(compfolder, currentimeidx, consecutivecenters, resolution, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, zerovortradius2, PressureVertLevelDimName, ZVarName)
    if(!is.na(Bparam) && (Bparam<Bthreshold)){BSymmetric <- TRUE} else {BSymmetric <- FALSE} 
    conditions <- c(conditions, BSymmetric)
    values <- c(values, Bparam)
  } 
  
  if((1 %in% HartConditionsTocheck)==TRUE && BSymmetric==TRUE){
    
    if((2 %in% HartConditionsTocheck)==TRUE && (3 %in% HartConditionsTocheck)==TRUE){
      TW <- ThermalWind(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW=TRUE, upperTW=TRUE, PressureVertLevelDimName, ZVarName)
      conditions <- c(conditions, TW$LTW>0, TW$UTW>0)
      values <- c(values, TW$LTW, TW$UTW)
    } 
    
    if((2 %in% HartConditionsTocheck)==TRUE && (3 %in% HartConditionsTocheck)==FALSE){
      TW <- ThermalWind(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW=TRUE, upperTW=FALSE, PressureVertLevelDimName, ZVarName)
      conditions <- c(conditions, TW$LTW>0)  
      values <- c(values, TW$LTW)
    } 
    
    if((2 %in% HartConditionsTocheck)==FALSE && (3 %in% HartConditionsTocheck)==TRUE){
      TW <- ThermalWind(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW=FALSE, upperTW=TRUE, PressureVertLevelDimName, ZVarName)
      conditions <- c(conditions, TW$LTW>0)  
      values <- c(values, TW$UTW)
    } 
    
    if((4 %in% HartConditionsTocheck)==TRUE){
      if((2 %in% HartConditionsTocheck)==TRUE && (3 %in% HartConditionsTocheck)==TRUE){ conditions <- c(conditions, TW$LTW>TW$UTW) } else{
        stop('Requesting the 4th Hart condition while allowing negative thermal wind values is physically inconsistent for tropical-like cyclones.')
      }
    } 
    
    result <- list()
    if(sum(conditions)==length(HartConditionsTocheck)){result$logical <- TRUE} else{result$logical <- FALSE}
    
    triad0 <- rep(NA, 3)
    triad0[c(1,2,3) %in% HartConditionsTocheck] <- values
    result$triad <- triad0
  } else if((1 %in% HartConditionsTocheck)==TRUE && BSymmetric==FALSE){
    result <- list()
    result$logical <- FALSE
    
    triad0 <- rep(NA, 3)
    result$triad <- triad0
  } else if((1 %in% HartConditionsTocheck)==FALSE){
    
    if((2 %in% HartConditionsTocheck)==TRUE && (3 %in% HartConditionsTocheck)==TRUE){
      TW <- ThermalWind(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW=TRUE, upperTW=TRUE, PressureVertLevelDimName, ZVarName)
      conditions <- c(conditions, TW$LTW>0, TW$UTW>0)
      values <- c(values, TW$LTW, TW$UTW)
    } 
    
    if((2 %in% HartConditionsTocheck)==TRUE && (3 %in% HartConditionsTocheck)==FALSE){
      TW <- ThermalWind(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW=TRUE, upperTW=FALSE, PressureVertLevelDimName, ZVarName)
      conditions <- c(conditions, TW$LTW>0)  
      values <- c(values, TW$LTW)
    } 
    
    if((2 %in% HartConditionsTocheck)==FALSE && (3 %in% HartConditionsTocheck)==TRUE){
      TW <- ThermalWind(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW=FALSE, upperTW=TRUE, PressureVertLevelDimName, ZVarName)
      conditions <- c(conditions, TW$LTW>0)  
      values <- c(values, TW$UTW)
    } 
    
    if((4 %in% HartConditionsTocheck)==TRUE){
      if((2 %in% HartConditionsTocheck)==TRUE && (3 %in% HartConditionsTocheck)==TRUE){ conditions <- c(conditions, TW$LTW>TW$UTW) } else{
        stop('Requesting the 4th Hart condition while allowing negative thermal wind values is physically inconsistent for tropical-like cyclones.')
      }
    } 
    
    result <- list()
    if(sum(conditions)==length(HartConditionsTocheck)){result$logical <- TRUE} else{result$logical <- FALSE}
    
    triad0 <- rep(NA, 3)
    triad0[c(1,2,3) %in% HartConditionsTocheck] <- values
    result$triad <- triad0
  }
  
  return(result)
}
