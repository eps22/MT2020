getTrackData <- function(folder, complete=NULL, save=NULL, timestepdt=NULL, resolution=NULL, CalculateZeroVortRadiusThreshold='zero', CalculateZeroVortRadiusNPoints=30, ZeroVortRadiusMaxAllowedAsymm=300){
  
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
 
  if(is.null(save) | save==FALSE){ 
    if(is.null(complete) | complete==FALSE){
      return(abc)
    }
  }

  if(!is.null(complete) && complete==TRUE){
    if(is.null(resolution)){resolution <- 9}
    if(is.null(timestepdt)){timestepdt <- 1}
    source('./TrackingAlgorithm/isMedicane.R')
    source('./TrackingAlgorithm/TrackCenterCandidate.R')
    source('./TrackingAlgorithm/FulfillsHart.R')
    source('./TrackingAlgorithm/conditionalSLP.R')
    source('./TrackingAlgorithm/CalculateZeroVortRadius.R')
    source('./TrackingAlgorithm/ThermalSymmB.R')
    source('./TrackingAlgorithm/ThermalWind.R')
    source('./TrackingAlgorithm/dd.R')
    
    library(oce)

    ncinuv10MOD <- nc_open(paste(folder,'/output/outputfile-uvmet10-MOD.nc',sep=''))
    ncinuv10U <- nc_open(paste(folder,'/output/outputfile-uvmet10-U.nc',sep=''))
    ncinuv10V <- nc_open(paste(folder,'/output/outputfile-uvmet10-V.nc',sep=''))
    uv10MOD <- ncvar_get(ncinuv10MOD, 'uvmet10')
    uv10U <- ncvar_get(ncinuv10U, 'uvmet10')
    uv10V <- ncvar_get(ncinuv10V, 'uvmet10')
    
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

      uv10MODt <- uv10MOD[,,timestamp]
      centerwinds <- c(centerwinds, round(uv10MODt[point[1],point[2]],2))

      uv10Ut <- uv10U[,,timestamp]
      uv10Vt <- uv10V[,,timestamp]
      slpt <- slp[,,timestamp]
      uv10x <- uv10Ut
      uv10y <- uv10Vt
      uv10curl <- curl(uv10x, uv10y, x=seq(1,dim(slpt)[1], by=1)*resolution,
                       y=seq(1,dim(slpt)[2],by=1)*resolution, geographical = FALSE, method = 1)$curl

      zerovortradius <- CalculateZeroVortRadius(resolution, uv10curl, point, threshold=CalculateZeroVortRadiusThreshold, npoints=CalculateZeroVortRadiusNPoints, maxasymm = ZeroVortRadiusMaxAllowedAsymm, minsymmdirs=6)
      zerovortradiuseq <- c(zerovortradiuseq, round(zerovortradius,1))
      INradiuspoints <- which((dd(slp[,,1],point)*resolution)<zerovortradius,arr.ind=TRUE)

      im <- isMedicane(compfolder=folder, center=point, currentimeidx=timestamp, timestep=timestepdt, resolution=resolution, HartConditionsTocheck=c(1,2,3,4), Blowerpressurelevel=900, Bupperpressurelevel=600, Bmultiplemeasure='max', Bdirections=4, Bthreshold=10, LTWlowerpressurelevel=900, LTWupperpressurelevel=600, UTWlowerpressurelevel=600, UTWupperpressurelevel=300, zerovortradius)
      
      bs <- c(bs, round(as.numeric(im$triad[1]),1))
      vtls <- c(vtls, round(as.numeric(im$triad[2]),1))
      vtus <- c(vtus, round(as.numeric(im$triad[3]),1))

      uv10MODINradiuspoints <- uv10MODt[INradiuspoints]
      maxINradiuswinds <- c(maxINradiuswinds, round(max(uv10MODINradiuspoints),2)) 
      minINradiuswinds <- c(minINradiuswinds, round(min(uv10MODINradiuspoints),2)) 
      
      maxwindpoint <- INradiuspoints[which(uv10MODINradiuspoints==max(uv10MODINradiuspoints), arr.ind=TRUE),]
      maxwindradius <- c(maxwindradius, round(dd(point,maxwindpoint)*resolution,1)) 
     
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
      minslpcenterds <- c(minslpcenterds, round(dd(INradiuspoints[slpposition,],point)*resolution,1))
      
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
