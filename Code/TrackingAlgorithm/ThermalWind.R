ThermalWind <- function(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW, upperTW, PressureVertLevelDimName, ZVarName){
  
  library(ncdf4)
  
  if(lowerTW==TRUE){
  
    zfieldnc <- paste(compfolder,'/output/outputfile-z.nc',sep='')
    
    zheight <- nc_open(zfieldnc)
    zz <- ncvar_get(zheight,ZVarName,start = c(1,1,1,currentimeidx), count=c(-1,-1,-1,1))
   
    #Get requested plevels indices
    plevs <- zheight$dim[[PressureVertLevelDimName]]$vals
  
    rangelowidx1 <- as.numeric(which(plevs==LTWlowerpressurelevel,arr.ind = TRUE))
    rangelowidx2 <- as.numeric(which(plevs==LTWupperpressurelevel,arr.ind = TRUE))
    
    deltaZlow <- c()
    for(plevidx in rangelowidx1:rangelowidx2){
      zPLEVEL <- zz[,,plevidx]
      
      currentcenter <- consecutivecenters[2,]
        
      INradiuspoints <- which((dd(zPLEVEL,currentcenter)*resolution)<zerovortradius2,arr.ind=TRUE)
      zVradiuskm <- c()
      for(i in 1:nrow(INradiuspoints)){
        zVradiuskm <- c(zVradiuskm,zPLEVEL[INradiuspoints[i,1],INradiuspoints[i,2]])
      }
      zVradiuskmMIN <- min(zVradiuskm)
      zVradiuskmMINposition <- INradiuspoints[which(zVradiuskm==min(zVradiuskmMIN),arr.ind = TRUE),]
      if(length(zVradiuskmMINposition)>2){zVradiuskmMINposition <- zVradiuskmMINposition[1,,drop=FALSE]}
      zVradiuskmMAX <- max(zVradiuskm)
      zVradiuskmMAXposition <- INradiuspoints[which(zVradiuskm==max(zVradiuskmMAX),arr.ind = TRUE),]
      if(length(zVradiuskmMAXposition)>2){zVradiuskmMAXposition <- zVradiuskmMAXposition[1,,drop=FALSE]}
      maxmindistance <- abs(dd(zVradiuskmMAXposition, zVradiuskmMINposition))#*resolution
      plevdeltaZ <- (zVradiuskmMAX - zVradiuskmMIN)/maxmindistance
      
      deltaZlow <- c(deltaZlow,plevdeltaZ)
    }
    
    LOWERVT <- as.numeric(lm(as.numeric(deltaZlow) ~ log(plevs[rangelowidx1:rangelowidx2]))$coefficients[2])
    
  } else{LOWERVT <- NA}
 
  if(upperTW==TRUE){
    zfieldnc <- paste(compfolder,'/output/outputfile-z.nc',sep='')
    
    zheight <- nc_open(zfieldnc)
    zz <- ncvar_get(zheight,ZVarName,start = c(1,1,1,currentimeidx), count=c(-1,-1,-1,1))
    
    #Get requested plevels indices
    plevs <- zheight$dim[[PressureVertLevelDimName]]$vals
    
    rangehighidx1 <- as.numeric(which(plevs==UTWlowerpressurelevel,arr.ind = TRUE))
    rangehighidx2 <- as.numeric(which(plevs==UTWupperpressurelevel,arr.ind = TRUE))
    
    deltaZhigh <- c()
    for(plevidx in rangehighidx1:rangehighidx2){
      zPLEVEL <- zz[,,plevidx]
      
      currentcenter <- consecutivecenters[2,]
      
      INradiuspoints <- which((dd(zPLEVEL,currentcenter)*resolution)<zerovortradius2,arr.ind=TRUE)
      zVradiuskm <- c()
      for(i in 1:nrow(INradiuspoints)){
        zVradiuskm <- c(zVradiuskm,zPLEVEL[INradiuspoints[i,1],INradiuspoints[i,2]])
      }
      zVradiuskmMIN <- min(zVradiuskm)
      zVradiuskmMINposition <- INradiuspoints[which(zVradiuskm==min(zVradiuskmMIN),arr.ind = TRUE),]
      if(length(zVradiuskmMINposition)>2){zVradiuskmMINposition <- zVradiuskmMINposition[1,,drop=FALSE]}
      zVradiuskmMAX <- max(zVradiuskm)
      zVradiuskmMAXposition <- INradiuspoints[which(zVradiuskm==max(zVradiuskmMAX),arr.ind = TRUE),]
      if(length(zVradiuskmMAXposition)>2){zVradiuskmMAXposition <- zVradiuskmMAXposition[1,,drop=FALSE]}
      maxmindistance <- abs(dd(zVradiuskmMAXposition, zVradiuskmMINposition))#*resolution
      plevdeltaZ <- (zVradiuskmMAX - zVradiuskmMIN)/maxmindistance
      
      deltaZhigh <- c(deltaZhigh,plevdeltaZ)
    }
    
    UPPERVT <- as.numeric(lm(as.numeric(deltaZhigh) ~ log(plevs[rangehighidx1:rangehighidx2]))$coefficients[2])
    
  } else{UPPERVT <- NA}
  
  TW <- list()
  TW$LTW <- LOWERVT
  TW$UTW <- UPPERVT
  return(TW)
  
} 
