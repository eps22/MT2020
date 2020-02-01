expandTrackData <- function(folder, resolution, trackingdf, expandN){ #returns extended DF of dots
  
  dottss0 <- trackingdf[,c('timestep','x','y')]

  #Read SLP file
  slp <- ncvar_get(nc_open(paste(getwd(),'/',folder,'/output/outputfile-slp.nc',sep='')),'slp')
  
  initial <- min(dottss0$timestep)
  if((initial-expandN)<1){initial2 <- 1} else{initial2 <- initial-expandN}
  final <- max(dottss0$timestep)
  if(dim(slp)[3]<(final+expandN)){final2 <- dim(slp)[3]} else{final2 <- final+expandN}
  
  dottss1 <- trackingdf[,c('timestep','min.inZVradius.slp.x.position','min.inZVradius.slp.y.position')]
  colnames(dottss1) <- colnames(dottss0)
  
  for(i in c(initial2:(initial-1),(final+1):final2)){
    minidx <- which.min(abs(trackingdf$timestep-i))
    if(length(minidx)>1){minidx <- minidx[1]}
    radius <- trackingdf[minidx,'outer.radius.km']
    center <- trackingdf[minidx,c('x','y')]
    INradiuspoints <- which((dd(slp[,,1],as.numeric(center))*resolution)<radius,arr.ind=TRUE)
    slpT <- slp[,,i]
    minINradiusslpposition <- INradiuspoints[which(slpT[INradiuspoints]==min(slpT[INradiuspoints]), arr.ind = TRUE),]
    dottss1 <- rbind(dottss1,c(i,minINradiusslpposition[1],minINradiusslpposition[2]))
  }
  
  #Add 400000 to all SLP dots timesteps so that they don't get connected with medicane center dots
  dottss1$timestep <- dottss1$timestep#+400000
  
  #dottss2 <- rbind(dottss0,dottss1)
  
  return(dottss1)
  
}
