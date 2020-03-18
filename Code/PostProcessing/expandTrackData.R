expandTrackData <- function(folder, resolution, trackingdf, expandN){ #returns extended DF of dots
  
  abc <- readLines(con = './TrackingAlgorithm/FindMedicanes.namelist')
  abc <- abc[unlist(sapply(trimws(abc),function(x){nchar(x)>1},USE.NAMES = FALSE))]
  SLPVarName <- trimws(strsplit(abc[11],'=')[[1]][2])

  dottss0 <- trackingdf[,c('timestep','x','y')]

  #Read SLP file
  slp <- ncvar_get(nc_open(paste(getwd(),'/',folder,'/output/outputfile-slp.nc',sep='')),SLPVarName)
  
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
  
  return(dottss1)
  
}
