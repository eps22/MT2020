ThermalWind <- function(compfolder, currentimeidx, consecutivecenters, resolution, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, lowerTW, upperTW){
  
  library(ncdf4)
  
  if(lowerTW==TRUE | upperTW==TRUE){
    
    rangelow <- c(LTWlowerpressurelevel, LTWupperpressurelevel) #hPa
    rangehigh <- c(UTWlowerpressurelevel, UTWupperpressurelevel) #hPa
      
    if(as.numeric(rangelow[1])<as.numeric(rangelow[2])){cat('Selected lower pressure level for the calculation of lower tropospheric thermal wind calculation is above the upper one and may lead to inconsistent medicanes detections')}
    if(as.numeric(rangehigh[1])<as.numeric(rangehigh[2])){cat('Selected lower pressure level for the calculation of upper tropospheric thermal wind calculation is above the upper one and may lead to inconsistent medicanes detections')}
      
  
    zfieldnc <- paste(compfolder,'/output/outputfile-z.nc',sep='')
    
    zheight <- nc_open(zfieldnc)
    zz <- ncvar_get(zheight,'height',start = c(1,1,1,currentimeidx), count=c(-1,-1,-1,1))
   
    #Get requested plevels indices
    plevs <- zheight$dim$interp_level$vals
  
    
    rangelowidx1 <- as.numeric(which(plevs==rangelow[1],arr.ind = TRUE))
    rangelowidx2 <- as.numeric(which(plevs==rangelow[2],arr.ind = TRUE))
    rangehighidx1 <- as.numeric(which(plevs==rangehigh[1],arr.ind = TRUE))
    rangehighidx2 <- as.numeric(which(plevs==rangehigh[2],arr.ind = TRUE))
    
    plevs <- plevs
    
    if(length(rangelowidx1)==0){stop('Error: selected lower pressure level for the calculation of lower tropospheric thermal wind is not a valid pressure level.')}
    if(length(rangelowidx2)==0){stop('Error: selected higher pressure level for the calculation of lower tropospheric thermal wind is not a valid pressure level.')}
    if(length(rangehighidx1)==0){stop('Error: selected lower pressure level for the calculation of upper tropospheric thermal wind is not a valid pressure level.')}
    if(length(rangehighidx2)==0){stop('Error: selected higher pressure level for the calculation of upper tropospheric thermal wind is not a valid pressure level.')}
    
    deltaZ <- c()
    maxmindistances <- c()
    for(plevidx in 1:length(plevs)){
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
      
      deltaZ <- c(deltaZ,plevdeltaZ)
    }
  
  }
  
  if(lowerTW==TRUE){
    LOWERVT <- as.numeric(lm(as.numeric(rev(deltaZ[(rangelowidx1+1):(rangelowidx2+1)])) ~ log(rev(plevs[rangelowidx1:rangelowidx2])))$coefficients[2])
  } else{LOWERVT <- NA}
  
  if(upperTW==TRUE){
    UPPERVT <- as.numeric(lm(as.numeric(rev(deltaZ[(rangehighidx1+1):(rangehighidx2+1)])) ~ log(rev(plevs[rangehighidx1:rangehighidx2])))$coefficients[2])
  } else{UPPERVT <- NA}
  
  TW <- list()
  TW$LTW <- LOWERVT
  TW$UTW <- UPPERVT
  return(TW)
  
} 
