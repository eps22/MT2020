PlotTrack <- function(folder, adjust, complete, connect, SLPexpandN=NULL, dthreshold=NULL, dtthreshold=NULL){
 	
  library(ncdf4)
  source('./PostProcessing/getTrackData.R')
  source('./PostProcessing/MatrixPlot.R')
  if(is.null(dthreshold)){dthreshold <- 100}
  if(is.null(dtthreshold)){dtthreshold <- 1}
 
  abc <- readLines(con = './TrackingAlgorithm/FindMedicanes.namelist')
  abc <- abc[unlist(sapply(trimws(abc),function(x){nchar(x)>1},USE.NAMES = FALSE))]

  Resolution <- as.numeric(strsplit(abc[3],'=')[[1]][2])

  tracklist <- getTrackData(paste(getwd(),'/',folder,sep=''), complete = FALSE, save = FALSE)
  trackingdf <- tracklist$trackdf
  lon <- tracklist$lon
  lat <- tracklist$lat
  
  zeromatrix <- matrix(0,nrow = nrow(lon), ncol = ncol(lat))
  
  dots <- trackingdf[,c('timestep','x','y')]
  ds <- 1
  dcol <- 'blue'
  dthreshold <- dthreshold
  dtthreshold <- dtthreshold

  if(complete==TRUE){
    source('./PostProcessing/expandTrackData.R')
    trackingdfcomp <- getTrackData(paste(getwd(),'/',folder,sep=''), complete = TRUE, save = FALSE)
    dots2 <- expandTrackData(folder,Resolution,trackingdfcomp,expandN=SLPexpandN)
    ds2 <- 3
    dcol2 <- 'red'
    dthreshold2 <- dthreshold
    dtthreshold2 <- dtthreshold
  } else{
    dots2 <- NULL
    ds2 <- 3
    dcol2 <- 'red'
    dthreshold2 <- dthreshold
    dtthreshold2 <- dtthreshold

  }
  
   if(adjust==TRUE){minlon <- min(trackingdf$Lon) - 1} else{minlon <- min(lon)}
   if(adjust==TRUE){maxlon <- max(trackingdf$Lon) + 1} else{maxlon <- max(lon)}
   if(adjust==TRUE){minlat <- min(trackingdf$Lat) - 1} else{minlat <- min(lat)}
   if(adjust==TRUE){maxlat <- max(trackingdf$Lat) + 1} else{maxlat <- max(lat)}

   print(folder)
   infofile <- paste(getwd(),'/',folder,'/output/outputfile-slp.nc',sep='')
   print(infofile)
   fileinfo <- system(paste('ncdump -h ',infofile, '| more',sep=''), wait=FALSE, intern=TRUE)
   projectioninfo <- fileinfo[grep('projection',fileinfo)]
   if(length(projectioninfo)>0){
      if(grep('LambertConformal',projectioninfo)==1){
         reflon <- as.character(round(as.numeric(strsplit(strsplit(projectioninfo,'stand_lon=')[[1]][2],',')[[1]][1]),3))
         reflat <- as.character(round(as.numeric(strsplit(strsplit(projectioninfo,'moad_cen_lat=')[[1]][2],',')[[1]][1]),3))
         crs=paste('+proj=laea +lat_0=',reflat,' +lon_0=',reflon,' +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',sep='')
      }
   } else{
      crs <- "+proj=longlat +datum=WGS84 +no_defs"
   }

  foldername <- tail(strsplit(folder,'/')[[1]],1)
  
  cairo_pdf(filename = paste(folder,'/',foldername,'-trackplot.pdf',sep=''), width = 10, height = 8)
  MatrixPlot(zeromatrix, lon, lat, crs=crs, xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat,
                            dots=dots, dots_symbol = ds, dots_cex = 1,
                            dots_color = dcol, dots_connect = connect, 
                            dots_connect_color='blue',resolution = Resolution,
                            dthreshold = dthreshold, dtthreshold = dtthreshold, legend = FALSE,
                            dthreshold2 = dthreshold, dtthreshold2 = dtthreshold,
                            dots2=dots2, dots2_symbol = ds2, dots2_cex = 1,
                            dots2_color = dcol2, dots2_connect = connect, dots2_connect_color = 'red')
  dev.off()
    
}
