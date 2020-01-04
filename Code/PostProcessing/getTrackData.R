getTrackData <- function(folder, complete=NULL, save=NULL, timestepdt=NULL, resolution=NULL){
  
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
  ncinslp <- nc_open(paste(folder,'/','output/outputfile-slp.nc',sep=''))
  lon <- ncvar_get(ncinslp, 'XLONG')
  lat <- ncvar_get(ncinslp, 'XLAT')
  slp <- ncvar_get(ncinslp,'slp')

  numofdates <- ncinslp$dim$Time$vals + 1
  initdate <- strsplit((strsplit(ncinslp$dim$Time$units, ' ')[[1]])[3], '-')[[1]]
  inithour <- strsplit((strsplit(ncinslp$dim$Time$units, ' ')[[1]])[4],':')[[1]]
  init <- c(initdate, inithour)
  alldates <- seq.POSIXt(do.call(ISOdate,as.list(init)), by="hour", length.out = max(numofdates))
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
 
  if(is.null(save)){ 
    return(abc)
  }
  if(save==FALSE){ 
    return(abc)
  }
 
  if(!is.null(complete) && complete==TRUE){
    if(is.null(resolution)){resolution <- 9}
    if(is.null(timestepdt)){timestepdt <- 1}
    source('./TrackingAlgorithm/isMedicane.R')
    source('./TrackingAlgorithm/TrackCenterCandidate.R')
    source('./TrackingAlgorithm/FulfillsHart.R')
    source('./TrackingAlgorithm/conditionalSLP.R')
    source('./TrackingAlgorithm/ThermalSymmB.R')
    source('./TrackingAlgorithm/ThermalWind.R')
    source('./TrackingAlgorithm/dd.R')
    
    ncinuv10MOD <- nc_open(paste(folder,'/','output/outputfile-uvmet10-MOD.nc',sep=''))
    ncinuv10U <- nc_open(paste(folder,'/','output/outputfile-uvmet10-U.nc',sep=''))
    ncinuv10V <- nc_open(paste(folder,'/','output/outputfile-uvmet10-V.nc',sep=''))
    uv10MOD <- ncvar_get(ncinuv10MOD, 'uvmet10')
    uv10U <- ncvar_get(ncinuv10U, 'uvmet10')
    uv10V <- ncvar_get(ncinuv10V, 'uvmet10')
    
    bs <- c()
    vtls <- c()
    vtus <- c()
    centerwinds <- c()
    max100kmwinds <- c()
    min100kmwinds <- c()
    min100kmslps <- c()
    minslpcenterds <- c()
    categories <- c()
    maxwindradius <- c()
    zerovortradius <- c()
    for(i in 1:nrow(trackingdf)){
      timestamp <- trackingdf[i,'timestep']
      point <- c(trackingdf[i,'x'], trackingdf[i,'y'])
      im <- isMedicane(compfolder=folder, center=point, currentimeidx=timestamp, timestep=timestepdt, resolution=resolution, HartConditionsTocheck=c(1,2,3,4), Blowerpressurelevel=900, Bupperpressurelevel=600, Bradius=100, Bmultiplemeasure='max', Bdirections=4, Bthreshold=10, LTWlowerpressurelevel=900, LTWupperpressurelevel=600, UTWlowerpressurelevel=600, UTWupperpressurelevel=300, TWradius=100)
      
      bs <- c(bs, round(as.numeric(im$triad[1]),1))
      vtls <- c(vtls, round(as.numeric(im$triad[2]),1))
      vtus <- c(vtus, round(as.numeric(im$triad[3]),1))
      
      INradiuspoints <- which((dd(slp[,,1],point)*resolution)<100,arr.ind=TRUE)
      uv10MODt <- uv10MOD[,,timestamp]
      uv10Ut <- uv10U[,,timestamp]
      uv10Vt <- uv10V[,,timestamp]
      slpt <- slp[,,timestamp]
      uv10x <- uv10Ut
      uv10y <- uv10Vt
      uv10curl <- curl(uv10x, uv10y, x=seq(1,dim(slpt)[1], by=1)*resolution, 
                       y=seq(1,dim(slpt)[2],by=1)*resolution, geographical = FALSE, method = 1)$curl

      centerwinds <- c(centerwinds, round(uv10MODt[point[1],point[2]],2))
      
      uv10MODINradiuspoints <- uv10MODt[INradiuspoints]
      max100kmwinds <- c(max100kmwinds, round(max(uv10MODINradiuspoints),2)) 
      min100kmwinds <- c(min100kmwinds, round(min(uv10MODINradiuspoints),2)) 
      
      maxwindpoint <- INradiuspoints[which(uv10MODINradiuspoints==max(uv10MODINradiuspoints), arr.ind=TRUE),]
      maxwindradius <- c(maxwindradius, round(dd(point,maxwindpoint)*resolution,1)) 
     
      center <- point

      if((dim(uv10Ut)[1]-center[1])<31){dright <- NA} else{valuesright <- uv10curl[center[1]:(center[1]+30),center[2]]; dright <- (min(which(valuesright<0.5))-0.5)*resolution}
      if(center[1]<31){dleft <- NA} else{valuesleft <- uv10curl[(center[1]-30):center[1],center[2]]; dleft <- (max(which(valuesleft<0.5))+0.5)*resolution}
      if((dim(uv10Ut)[2]-center[2])<31){dup <- NA} else{valuesup <- uv10curl[center[1],center[2]:(center[2]+30)]; dup <- (min(which(valuesup<0.5))-0.5)*resolution}
      if(center[2]<31){ddown <- NA} else{valuesdown <- uv10curl[center[1],(center[2]-30):center[2]]; ddown <- (max(which(valuesdown<0.5))+0.5)*resolution}
      
      vectorofds <- c(dright, dleft, dup, ddown)
      if(sum(!is.na(vectorofds))>1){ zerovortd <- mean(vectorofds, na.rm = TRUE)} else{zerovortd <- NA}
      zerovortradius <- c(zerovortradius, round(zerovortd,1))
      
      SLPINradiuspoints <- slpt[INradiuspoints]
      
      min100kmslps <- c(min100kmslps, round(min(SLPINradiuspoints),2))
      minslpcenterds <- c(minslpcenterds, round(dd(INradiuspoints[which(SLPINradiuspoints==min(SLPINradiuspoints), arr.ind = TRUE),],point)*resolution,1))
      
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
    zerovortradius <- data.frame(zerovortradius)
    bs <- data.frame(bs)
    vtls <- data.frame(vtls)
    vtus <- data.frame(vtus)
    centerwinds <- data.frame(centerwinds)
    max100kmwinds <- data.frame(max100kmwinds)
    min100kmwinds <- data.frame(min100kmwinds)
    min100kmslps <- data.frame(min100kmslps)
    minslpcenterds <- data.frame(minslpcenterds)
    categories <- data.frame(categories)
    
    prevnames <- colnames(trackingdf)
    trackingdf <- cbind.data.frame(trackingdf, maxwindradius, zerovortradius, bs, vtls, vtus, max100kmwinds, centerwinds, min100kmwinds, min100kmslps, minslpcenterds, categories) 
    colnames(trackingdf) <- c(prevnames, 'inner.radius.km', 'outer.radius.km' ,'B.m', 'LTW','UTW', 'max.100km.wind.kmh', 'center.wind.kmh', 'min.100km.wind.kmh', 'min.100km.slp.hPa','minSLP.center.distance.km', 'Saffir.Simpson.scale.category')
  }

  if(!is.null(save) && save==TRUE){
    write.csv(x=trackingdf,file=paste(folder,'/trackingdf.csv',sep=''))
  }
  
}
