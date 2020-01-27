
completeinput <- commandArgs(trailingOnly = FALSE)

args <- completeinput[2:(length(completeinput)-1)]
print(args)

filestring <- args[grep('file=',args)]
if(length(filestring)==0){
  stop('No file was specified. Please check the input arguments string.')
} else{
  file=strsplit(strsplit(filestring,'file=')[[1]][2],' ')[[1]][1]
}
print(file)

fieldstring <- args[grep('field=',args)]  
if(length(fieldstring)==0){
  stop('No field was specified. Please check the input arguments string.')
} else{
  field=strsplit(strsplit(fieldstring,'field=')[[1]][2],' ')[[1]][1]
}
print(field)

vlevelstring <- args[grep('vlevel=',args)]
if(length(vlevelstring)==0){
  vlevel <- NA
} else{
  vlevel <- as.numeric(strsplit(strsplit(vlevelstring,'vlevel=')[[1]][2],' ')[[1]][1])
}
print(vlevel)

timestepstring <- args[grep('timestep=',args)]
if(length(timestepstring)==0){
  timestep <- NA
} else{
  timestep <- as.numeric(strsplit(strsplit(timestepstring,'timestep=')[[1]][2],' ')[[1]][1])
}
print(timestep)

xminstring <- args[grep('xmin=',args)]
if(length(xminstring)==0){
  xmin <- NA
} else{
  xmin <- as.numeric(strsplit(strsplit(xminstring,'xmin=')[[1]][2],' ')[[1]][1])
}
print(xmin)

xmaxstring <- args[grep('xmax=',args)]
if(length(xmaxstring)==0){
  xmax <- NA
} else{
  xmax <- as.numeric(strsplit(strsplit(xmaxstring,'xmax=')[[1]][2],' ')[[1]][1])
}
print(xmax)

yminstring <- args[grep('ymin=',args)]
if(length(yminstring)==0){
  ymin <- NA
} else{
  ymin <- as.numeric(strsplit(strsplit(yminstring,'ymin=')[[1]][2],' ')[[1]][1])
}
print(ymin)

ymaxstring <- args[grep('ymax=',args)]
if(length(ymaxstring)==0){
  ymax <- NA
} else{
  ymax <- as.numeric(strsplit(strsplit(ymaxstring,'ymax=')[[1]][2],' ')[[1]][1])
}
print(ymax)

xticksintstring <- args[grep('xticksint=',args)]
if(length(xticksintstring)==0){
  xticksint <- NULL
} else{
  xticksint <- as.numeric(strsplit(strsplit(xticksintstring,'xticksint=')[[1]][2],' ')[[1]][1])
}
print(xticksint)


yticksintstring <- args[grep('yticksint=',args)]
if(length(yticksintstring)==0){
  yticksint <- NULL
} else{
  yticksint <- as.numeric(strsplit(strsplit(yticksintstring,'yticksint=')[[1]][2],' ')[[1]][1])
}
print(yticksint)

colorpalettestring <- args[grep('colorpalette=',args)]
if(length(colorpalettestring)==0){
  colorpalette <- NA
} else{
  colorpalette <- strsplit(strsplit(colorpalettestring,'colorpalette=')[[1]][2],' ')[[1]][1]
  library(RColorBrewer)
  palette <- rev(colorRampPalette(brewer.pal(11, colorpalette))(100))
}
print(colorpalette)

barlimstring <- args[grep('barlims=',args)]
if(length(barlimstring)==0){
  barlims <- NULL
} else{
  barlims <- eval(parse(text=strsplit(strsplit(barlimstring,'barlims=')[[1]][2],' ')[[1]][1]))
}
print(barlims)

barbreakstring <- args[grep('barbreaks=',args)]
if(length(barbreakstring)==0){
  barbreaks <- NULL
} else{
  barbreaks <- eval(parse(text=strsplit(strsplit(barbreakstring,'barbreaks=')[[1]][2],' ')[[1]][1]))
  barbreaks[1] <- barbreaks[1]+0.1
  barbreaks[length(barbreaks)] <- barbreaks[length(barbreaks)]-0.1
}
print(barbreaks)

barlabstring <- args[grep('barlabs=',args)]
if(length(barlabstring)==0){
  barlabs <- NULL
} else{
  barlabs <- as.character(eval(parse(text=strsplit(strsplit(barlabstring,'barlabs=')[[1]][2],' ')[[1]][1])))
}
print(barlabs)

outnamestring <- args[grep('outname=',args)]
if(length(outnamestring)==0){
  outname <- './plotfield-outmap.pdf'
} else{
  outname <- strsplit(strsplit(outnamestring,'outname=')[[1]][2],' ')[[1]][1]
}
print(outname)

library(ncdf4)

inputnc <- nc_open(file)
print(length(inputnc$dim))
if(is.na(timestep)){
  stop('No timestep was specified. Please check the input arguments string.')
} else{
  if(length(inputnc$dim)==4 & is.na(vlevel)){
    stop('No vertical level was specified for the 4D field you are trying to plot. Please check the input arguments string.')
  } 
  if(length(inputnc$dim)==4 & !is.na(vlevel)){
    inputfield <- ncvar_get(inputnc,field,start=c(1,1,vlevel,timestep),count=c(-1,-1,1,1))
  }
  if(length(inputnc$dim)==3){
    inputfield <- ncvar_get(inputnc,field,start=c(1,1,timestep),count=c(-1,-1,1))
  }
}
print(dim(inputfield))

lon <- ncvar_get(inputnc,'XLONG')
lat <- ncvar_get(inputnc,'XLAT')

fileinfo <- system(paste('ncdump -h ',file, '| more',sep=''), wait=FALSE, intern=TRUE)
projectioninfo <- fileinfo[grep('projection',fileinfo)]
reflon <- as.character(round(as.numeric(strsplit(strsplit(projectioninfo,'stand_lon=')[[1]][2],',')[[1]][1]),3))
reflat <- as.character(round(as.numeric(strsplit(strsplit(projectioninfo,'moad_cen_lat=')[[1]][2],',')[[1]][1]),3))

crs=paste('+proj=laea +lat_0=',reflat,' +lon_0=',reflon,' +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',sep='')

source('MatrixPlot.R')

if(sum(is.na(c(xmin,xmax,ymin,ymax)))>=1){
  warning('One or more missing values in the specified geographical limits. Default limits will be used.')
  if(!is.na(colorpalette)){
    cairo_pdf(filename = outname, width = 10, height = 8) 
    MatrixPlot(inputfield,lon,lat,crs,x_ticks_interval=xticksint,y_ticks_interval=yticksint,bar_color_scale=palette,bar_lims=barlims,bar_breaks=barbreaks,bar_labels=barlabs)
    dev.off()
  } else{
    cairo_pdf(filename = outname, width = 10, height = 8)
    MatrixPlot(inputfield,lon,lat,crs,x_ticks_interval=xticksint,y_ticks_interval=yticksint,bar_lims=bar_lims,bar_breaks=barbreaks,bar_labels=barlabs)
    dev.off()
  }
} else{
  if(!is.na(colorpalette)){
    cairo_pdf(filename = outname, width = 10, height = 8)
    MatrixPlot(inputfield,lon,lat,crs,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x_ticks_interval=xticksint,y_ticks_interval=yticksint,bar_color_scale=palette,bar_lims=barlims,bar_breaks=barbreaks,bar_labels=barlabs)
    dev.off()
  } else{
    cairo_pdf(filename = outname, width = 10, height = 8)
    MatrixPlot(inputfield,lon,lat,crs,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x_ticks_interval=xticksint,y_ticks_interval=yticksint,bar_breaks=barbreaks,bar_lims=barlims,bar_labels=barlabs)
    dev.off()
  }
}
