getTrackData <- function(folder, complete=NULL, save=NULL){
 
  abc <- readLines(con = './TrackingAlgorithm/FindMedicanes.namelist')
  abc <- abc[unlist(sapply(trimws(abc),function(x){nchar(x)>1},USE.NAMES = FALSE))]
  
  Resolution <- as.numeric(strsplit(abc[3],'=')[[1]][2])           
  TimestepDt <- as.numeric(strsplit(abc[4],'=')[[1]][2])
  LonDimName <- trimws(strsplit(abc[5],'=')[[1]][2])
  LonVarName <- trimws(strsplit(abc[6],'=')[[1]][2])
  LatDimName <- trimws(strsplit(abc[7],'=')[[1]][2])
  LatVarName <- trimws(strsplit(abc[8],'=')[[1]][2])
  TimeDimName <- trimws(strsplit(abc[9],'=')[[1]][2])
  PressureVertLevelDimName <- trimws(strsplit(abc[10],'=')[[1]][2])
  SLPVarName <- trimws(strsplit(abc[11],'=')[[1]][2])
  U10VarName <- trimws(strsplit(abc[12],'=')[[1]][2])
  V10VarName <- trimws(strsplit(abc[13],'=')[[1]][2])
  ZVarName <- trimws(strsplit(abc[14],'=')[[1]][2])
  CalculateZeroVortRadiusThreshold <- trimws(strsplit(abc[19],'=')[[1]][2])
  CalculateZeroVortRadiusNPoints <- as.numeric(strsplit(abc[20],'=')[[1]][2])
  ZeroVortRadiusMaxAllowedAsymm <- as.numeric(strsplit(abc[22],'=')[[1]][2])
  ZeroVortRadiusMinSymmDirs <- as.numeric(strsplit(abc[23],'=')[[1]][2])
  HartConditionsTocheck <- as.numeric(strsplit(trimws(strsplit(abc[30],'=')[[1]][2]),',')[[1]])
  Blowerpressurelevel <- as.numeric(strsplit(abc[31],'=')[[1]][2])
  Bupperpressurelevel <- as.numeric(strsplit(abc[32],'=')[[1]][2])
  Bmultiplemeasure <- trimws(strsplit(abc[33],'=')[[1]][2])
  Bdirections <- as.numeric(strsplit(abc[34],'=')[[1]][2])
  Bthreshold <- as.numeric(strsplit(abc[35],'=')[[1]][2])
  LTWlowerpressurelevel <- as.numeric(strsplit(abc[36],'=')[[1]][2])
  LTWupperpressurelevel <- as.numeric(strsplit(abc[37],'=')[[1]][2])
  UTWlowerpressurelevel <- as.numeric(strsplit(abc[38],'=')[[1]][2])
  UTWupperpressurelevel <- as.numeric(strsplit(abc[39],'=')[[1]][2])

  	
  trackfile <- list.files(folder)[grepl('-track.RData',list.files(folder))]
  trackinglist <- readRDS(paste(folder,'/',trackfile,sep=''))
  
  idx <- 0
  trackingdf <- do.call('rbind', lapply(trackinglist, function(x){
    idx <<- idx + 1
    if(length(x)==2){
      ab <- t(as.data.frame(c(x[1],x[2])))
    } else{
      ab <- apply(x, 2, function(y){c(y[1],y[2])})
    }
    trep <- as.data.frame(rep(idx,length(ab)/2))
    return(cbind(trep,ab))
  }))
  
  
  row.names(trackingdf) <- NULL
  colnames(trackingdf) <- c('timestep','x','y')
  trackingdf <- trackingdf[!is.na(trackingdf$x),]
 
  library(ncdf4) 
  ncinslp <- nc_open(paste(folder,'/output/outputfile-slp.nc',sep=''))
  lonsc <- ncvar_get(ncinslp, LonVarName)
  latsc <- ncvar_get(ncinslp, LatVarName)
  if(!is.null(dim(lonsc)) & length(dim(lonsc))>1){lon <- lonsc; lat <- latsc} else{
    lon <- matrix(rep(lonsc, times=length(latsc)),nrow=length(lonsc),ncol=length(latsc))
    lat <- t(matrix(rep(latsc, times=length(lonsc)),nrow=length(latsc),ncol=length(lonsc)))
  }
  slp <- ncvar_get(ncinslp, SLPVarName)

  numofdates <- ncinslp$dim[[TimeDimName]]$vals + 1
  initdate <- strsplit((strsplit(ncinslp$dim[[TimeDimName]]$units, ' ')[[1]])[3], '-')[[1]]
  inithour <- strsplit((strsplit(ncinslp$dim[[TimeDimName]]$units, ' ')[[1]])[4],':')[[1]]
  init <- c(initdate, inithour)
  alldates <- seq.POSIXt(do.call(ISOdate,as.list(init)), by="hour", length.out = max(numofdates))
  alldates <- alldates[numofdates]
  seldates <- alldates[trackingdf$timestep]
  seldates <- data.frame(seldates)
  
  lons <- c()
  lats <- c()
  slpmins <- c()
  for(i in 1:nrow(trackingdf)){
    timestep <- trackingdf[i,'timestep']
    point <- c(trackingdf[i,'x',drop=FALSE], trackingdf[i,'y',drop=FALSE])
    lons <- c(lons, round(lon[point$x,point$y],3))
    lats <- c(lats, round(lat[point$x,point$y],3))
    slpmins <- c(slpmins, round(slp[point$x,point$y,timestep],2))
  }
  lons <- data.frame(lons)
  lats <- data.frame(lats)
  slpmins <- data.frame(slpmins)
  
  trackingdf <- cbind.data.frame(seldates,trackingdf$timestep, trackingdf$x, lons, trackingdf$y, lats, slpmins)
  colnames(trackingdf) <- c('Date','timestep','x','Lon','y','Lat','center.slp.hPa')
  
  abc <- list()
  abc$trackdf <- trackingdf
  abc$lon <- lon
  abc$lat <- lat
 
  if(is.null(save) | save==FALSE){ 
    if(is.null(complete) | complete==FALSE){
      return(abc)
    }
  }

  if(!is.null(complete) && complete==TRUE){
    source('./TrackingAlgorithm/isMedicane.R')
    source('./TrackingAlgorithm/TrackCenterCandidate.R')
    source('./TrackingAlgorithm/FulfillsHart.R')
    source('./TrackingAlgorithm/conditionalSLP.R')
    source('./TrackingAlgorithm/CalculateZeroVortRadius.R')
    source('./TrackingAlgorithm/ThermalSymmB.R')
    source('./TrackingAlgorithm/ThermalWind.R')
    source('./TrackingAlgorithm/dd.R')
    
    library(oce)

    ncinuv10U <- nc_open(paste(folder,'/output/outputfile-uvmet10-U.nc',sep=''))
    ncinuv10V <- nc_open(paste(folder,'/output/outputfile-uvmet10-V.nc',sep=''))
    uv10U <- ncvar_get(ncinuv10U, U10VarName)
    uv10V <- ncvar_get(ncinuv10V, V10VarName)
    
    if(ncinuv10U$dim[[LonDimName]]$vals[2]>ncinuv10U$dim[[LonDimName]]$vals[1]){revx <- FALSE} else{revx <- TRUE}
    if(ncinuv10U$dim[[LatDimName]]$vals[2]>ncinuv10U$dim[[LatDimName]]$vals[1]){revy <- FALSE} else{revy <- TRUE}

    bs <- c()
    vtls <- c()
    vtus <- c()
    centerwinds <- c()
    maxINradiuswinds <- c()
    minINradiuswinds <- c()
    minINradiusslps <- c()
    minINradiusslpxpositions <- c()
    minINradiusslplons <- c()
    minINradiusslpypositions <- c()
    minINradiusslplats <- c()
    minslpcenterds <- c()
    categories <- c()
    maxwindradius <- c()
    zerovortradiuseq <- c()

    for(i in 1:nrow(trackingdf)){

      timestamp <- trackingdf[i,'timestep']
      point <- c(trackingdf[i,'x'], trackingdf[i,'y'])
      
      uv10Ut <- uv10U[,,timestamp]
      uv10Vt <- uv10V[,,timestamp]
      uv10MODt <- sqrt(uv10Ut*uv10Ut + uv10Vt*uv10Vt)
      
      centerwinds <- c(centerwinds, round(uv10MODt[point[1],point[2]],2))

      slpt <- slp[,,timestamp]
      xseq <- seq(1,dim(slpt)[1], by=1)
      if(revx==TRUE){xseq <- rev(xseq)}
      yseq <- seq(1,dim(slpt)[2],by=1)
      if(revy==TRUE){yseq <- rev(yseq)}
      uv10curl <- curl(uv10Ut, uv10Vt, x=xseq*Resolution,
                   y=yseq*Resolution, geographical = FALSE, method = 1)$curl
      
      zerovortradius <- CalculateZeroVortRadius(Resolution, uv10curl, point, threshold=CalculateZeroVortRadiusThreshold, npoints=CalculateZeroVortRadiusNPoints, maxasymm = ZeroVortRadiusMaxAllowedAsymm, minsymmdirs=ZeroVortRadiusMinSymmDirs)
      zerovortradiuseq <- c(zerovortradiuseq, round(zerovortradius,1))
      INradiuspoints <- which((dd(slp[,,1],point)*Resolution)<zerovortradius,arr.ind=TRUE)

      im <- isMedicane(compfolder=folder, center=point, currentimeidx=timestamp, timestep=TimestepDt, resolution=Resolution, HartConditionsTocheck=HartConditionsTocheck, Blowerpressurelevel=Blowerpressurelevel, Bupperpressurelevel=Bupperpressurelevel, Bmultiplemeasure=Bmultiplemeasure, Bdirections=Bdirections, Bthreshold=Bthreshold, LTWlowerpressurelevel=LTWlowerpressurelevel, LTWupperpressurelevel=LTWupperpressurelevel, UTWlowerpressurelevel=UTWlowerpressurelevel, UTWupperpressurelevel=UTWupperpressurelevel, zerovortradius, PressureVertLevelDimName, SLPVarName, ZVarName)
      
      bs <- c(bs, round(as.numeric(im$triad[1]),1))
      vtls <- c(vtls, round(as.numeric(im$triad[2]),1))
      vtus <- c(vtus, round(as.numeric(im$triad[3]),1))

      uv10MODINradiuspoints <- uv10MODt[INradiuspoints]
      maxINradiuswinds <- c(maxINradiuswinds, round(max(uv10MODINradiuspoints),2)) 
      minINradiuswinds <- c(minINradiuswinds, round(min(uv10MODINradiuspoints),2)) 
      
      maxwindpoint <- INradiuspoints[which(uv10MODINradiuspoints==max(uv10MODINradiuspoints), arr.ind=TRUE),]
      maxwindradius <- c(maxwindradius, round(dd(point,maxwindpoint)*Resolution,1)) 
     
      SLPINradiuspoints <- slpt[INradiuspoints]
     
      slpposition <- which(SLPINradiuspoints==min(SLPINradiuspoints), arr.ind = TRUE)
      if(length(slpposition)>1){slpposition <- slpposition[1]}
      minINradiusslpposition <- INradiuspoints[slpposition,]
      #print(minINradiusslpposition)
      minINradiusslpxpositions <- c(minINradiusslpxpositions, minINradiusslpposition[1])
      minINradiusslpypositions <- c(minINradiusslpypositions, minINradiusslpposition[2])
      minINradiusslplons <- c(minINradiusslplons, round(lon[minINradiusslpposition[1],minINradiusslpposition[2]],3))
      minINradiusslplats <- c(minINradiusslplats, round(lat[minINradiusslpposition[1],minINradiusslpposition[2]],3))
      minINradiusslps <- c(minINradiusslps, round(min(SLPINradiuspoints),2))
      minslpcenterds <- c(minslpcenterds, round(dd(INradiuspoints[slpposition,],point)*Resolution,1))
      
      if(max(uv10MODINradiuspoints)<=62){category <- 'TD'} #<62 kmh: Tropical depression
      if(max(uv10MODINradiuspoints)>62 && max(uv10MODINradiuspoints)<=118){category <- 'TS'} #Tropical storm
      if(max(uv10MODINradiuspoints)>118 && max(uv10MODINradiuspoints)<=153){category <- 'H C1'} #Hurricane cat 1
      if(max(uv10MODINradiuspoints)>153 && max(uv10MODINradiuspoints)<=177){category <- 'H C2'} #Hurricane cat 2
      if(max(uv10MODINradiuspoints)>177 && max(uv10MODINradiuspoints)<=208){category <- 'H C3'} #Hurricane cat 3
      if(max(uv10MODINradiuspoints)>208 && max(uv10MODINradiuspoints)<=251){category <- 'H C4'} #Hurricane cat 4
      if(max(uv10MODINradiuspoints)>251){category <- 'H C5'} #Hurricane cat 5
      categories <- c(categories, category)

    }

    maxwindradius <- data.frame(maxwindradius)
    zerovortradiuseq <- data.frame(zerovortradiuseq)
    bs <- data.frame(bs)
    vtls <- data.frame(vtls)
    vtus <- data.frame(vtus)
    centerwinds <- data.frame(centerwinds)
    maxINradiuswinds <- data.frame(maxINradiuswinds)
    minINradiuswinds <- data.frame(minINradiuswinds)
    minINradiusslps <- data.frame(minINradiusslps)
    minINradiusslpxpositions <- data.frame(minINradiusslpxpositions)
    minINradiusslplons <- data.frame(minINradiusslplons)
    minINradiusslpypositions <- data.frame(minINradiusslpypositions)
    minINradiusslplats <- data.frame(minINradiusslplats)
    minslpcenterds <- data.frame(minslpcenterds)
    categories <- data.frame(categories)
    
    prevnames <- colnames(trackingdf)
    trackingdf <- cbind.data.frame(trackingdf, maxwindradius, zerovortradiuseq, bs, vtls, vtus, maxINradiuswinds, centerwinds, minINradiuswinds, minINradiusslpxpositions, minINradiusslplons, minINradiusslpypositions, minINradiusslplats, minINradiusslps, minslpcenterds, categories) 
    colnames(trackingdf) <- c(prevnames, 'inner.radius.km', 'outer.radius.km' ,'B.m', 'LTW','UTW', 'max.inZVradius.wind.kmh', 'center.wind.kmh', 'min.inZVradius.wind.kmh', 'min.inZVradius.slp.x.position','min.inZVradius.slp.Lon','min.inZVradius.slp.y.position','min.inZVradius.slp.Lat','min.inZVradius.slp.hPa','minSLP.center.distance.km', 'Saffir.Simpson.scale.category')
  }

  if(!is.null(save) && save==TRUE){
    foldername <- tail(strsplit(folder,'/')[[1]],1)
    write.csv(x=trackingdf,file=paste(folder,'/',foldername,'-trackingdf.csv',sep=''))
  } else{
    return(trackingdf)
  }
  
}
