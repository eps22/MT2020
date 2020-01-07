PlotTrack <- function(folder, connect, adjust, resolution=NULL, dthreshold=NULL, dtthreshold=NULL){
  
  library(ncdf4)
  source('./PostProcessing/getTrackData.R')
  source('./PostProcessing/MatrixPlotwDotsConnect.R')
  if(is.null(resolution)){resolution <- 9}
  if(is.null(dthreshold)){dthreshold <- 100}
  if(is.null(dtthreshold)){dtthreshold <- 1}
  
  tracklist <- getTrackData(paste(getwd(),'/',folder,sep=''), complete = FALSE, save = FALSE)
  trackingdf <- tracklist$trackdf
  lon <- tracklist$lon
  lat <- tracklist$lat
  
  zeromatrix <- matrix(0,nrow = nrow(lon), ncol = ncol(lat))
  
  dots <- trackingdf[,c('timestep','x','y')]

  if(adjust==TRUE){minlon <- min(lon[min(dots$x),]) - 1} else{minlon <- min(lon)}
  if(adjust==TRUE){maxlon <- max(lon[max(dots$x),]) + 1} else{maxlon <- max(lon)}
  if(adjust==TRUE){minlat <- min(lat[,min(dots$y)]) - 1} else{minlat <- min(lat)}
  if(adjust==TRUE){maxlat <- max(lat[,max(dots$y)]) + 1} else{maxlat <- max(lat)}
  
  infofile <- paste(getwd(),'/',folder,'/output/outputfile-slp.nc',sep='')
  fileinfo <- system(paste('ncdump -h ',infofile, '| more',sep=''), wait=FALSE, intern=TRUE)
  projectioninfo <- fileinfo[grep('projection',fileinfo)]
  reflon <- as.character(round(as.numeric(strsplit(strsplit(projectioninfo,'stand_lon=')[[1]][2],',')[[1]][1]),3))
  reflat <- as.character(round(as.numeric(strsplit(strsplit(projectioninfo,'moad_cen_lat=')[[1]][2],',')[[1]][1]),3))
  
  crs=paste('+proj=laea +lat_0=',reflat,' +lon_0=',reflon,' +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',sep='')
   
  print(crs)

  cairo_pdf(filename = paste(folder,'/trackplot.pdf',sep=''), width = 10, height = 8)
  MatrixPlotwDotsConnect(zeromatrix, lon, lat, crs=crs, xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat,
                        dots=dots, dots_symbol = 1, dots_cex=1,
                        dots_color = 'blue', dots_connect = connect, resolution = resolution,
                        dthreshold = dthreshold, dtthreshold = dtthreshold, legend = FALSE)
  dev.off()
  
}
