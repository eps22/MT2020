MatrixPlot <- function(field, lon, lat, crs, resolution = NULL, dthreshold = NULL, dtthreshold = NULL,
                          dthreshold2 = NULL, dtthreshold2 = NULL, dthreshold3 = NULL, dtthreshold3 = NULL,
                          xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, contours = NULL, contours_field = NULL , 
                          contours_spacing = NULL, contours_breaks = NULL, contours_color = NULL , contours_labels_size = NULL,
                          contours_labels_color = NULL, contours_lty = NULL, contours_lwd = NULL, 
                          bar_lims = NULL, bar_breaks = NULL, bar_labels = NULL, 
                          x_ticks_interval = NULL, y_ticks_interval = NULL, legend = TRUE,
                          bar_color_scale = NULL, field_units = NULL, main_title = NULL, completewithSLPdots = FALSE,
                          dots = NULL, dots_symbol = NULL, dots_cex = NULL, dots_color = NULL, dots_connect = NULL, dots_connect_color = NULL,
                          dots2 = NULL, dots2_symbol = NULL, dots2_cex = NULL, dots2_color = NULL, dots2_connect = NULL, dots2_connect_color = NULL,
                          dots3 = NULL, dots3_symbol = NULL, dots3_cex = NULL, dots3_color = NULL, dots3_connect = NULL, dots3_connect_color = NULL,
                          arrows_field_x = NULL, arrows_field_y = NULL, arrows_alpha = NULL, arrows_pixels_subsampling = NULL, arrows_length_factor = NULL){
  
  suppressMessages(library(rnaturalearth))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(sf))
  suppressMessages(library(akima))
  suppressMessages(library(reshape2))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggspatial))
  suppressMessages(library(maptools))
  suppressMessages(library(rgeos))
  suppressMessages(library(rgdal))
  suppressMessages(library(sp)) 
  suppressMessages(library(raster)) 
  suppressMessages(library(latticeExtra))
  suppressMessages(library(metR))
  suppressMessages(library(scales))
  
  dd <- function(point1,point2){
    if(is.null(nrow(point1)) && length(point1)==2){
      dist <- sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2)
      return(as.numeric(dist))
    }
    if(nrow(point1)==1){
      dist <- sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2)
      return(as.numeric(dist))
    }
    if(nrow(point1)>1 && ncol(point1) == 2){
      dist <- c()
      for(ll in 1:nrow(point1)){
        pointt <- point1[ll,]
        dist <- c(dist,sqrt((pointt[1]-point2[1])^2 + (pointt[2]-point2[2])^2))
      }
      return(as.numeric(dist))
    }
    if(nrow(point1)>1 && ncol(point1) > 1){
      dist <- matrix(0,nrow = nrow(point1),ncol=ncol(point1))
      for(ll in 1:nrow(point1)){
        for(bb in 1:ncol(point1)){
          pointt <- c(ll,bb)
          dist[ll,bb] <- sqrt((pointt[1]-point2[1])^2 + (pointt[2]-point2[2])^2)
        }
      }
      return(dist)
    }
  }
  
  isSequent <- function(point1, point2, resolution, dthreshold, dtthreshold){
    time1 <- point1[1]
    time2 <- point2[1]
    if(time2<=(time1+dtthreshold)){timesequent <- TRUE} else{timesequent <- FALSE}
    point1position <- c(point1[2],point1[3])
    point2position <- c(point2[2],point2[3])
    if((dd(point1position, point2position)*resolution)<dthreshold){spatialsequent <- TRUE} else{spatialsequent <- FALSE}
    if(timesequent && spatialsequent){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
  
  if(is.null(xmin)){xmin <- min(lon)}
  if(is.null(xmax)){xmax <- max(lon)}
  if(is.null(ymin)){ymin <- min(lat)}
  if(is.null(ymax)){ymax <- max(lat)}
  if(is.null(bar_color_scale)){bar_color_scale <- brewer.pal(11, "RdBu")}
  if(is.null(resolution)){resolution <- 9}
  if(is.null(dthreshold)){dthreshold <- 100}
  if(is.null(dtthreshold)){dtthreshold <- 1}
  if(is.null(field_units)){field_units <- ''}
  if(is.null(x_ticks_interval)){x_ticks_interval <- round((xmax-xmin)/5,0)}
  if(is.null(y_ticks_interval)){y_ticks_interval <- round((ymax-ymin)/5,0)}
  if(is.null(main_title)){main_title <- ''}
  if(is.null(dots_symbol)){dots_symbol <- 3}
  if(is.null(dots_cex)){dots_cex <- 4}
  if(is.null(dots_color)){dots_color <- 'black'}
  if(is.null(dots_connect_color)){dots_connect_color <- 'black'}
  if(is.null(dots2_symbol)){dots2_symbol <- 3}
  if(is.null(dots2_cex)){dots2_cex <- 4}
  if(is.null(dots2_color)){dots2_color <- 'black'}
  if(is.null(dots2_connect_color)){dots2_connect_color <- 'black'}
  if(is.null(dots3_symbol)){dots3_symbol <- 3}
  if(is.null(dots3_cex)){dots3_cex <- 4}
  if(is.null(dots3_color)){dots3_color <- 'black'}
  if(is.null(dots3_connect_color)){dots3_connect_color <- 'black'}
  if(is.null(arrows_alpha)){arrows_alpha <- 0.15}
  if(is.null(arrows_pixels_subsampling)){arrows_pixels_subsampling <- 2}
  if(is.null(arrows_length_factor)){arrows_length_factor <- 1.7}
  if(is.null(contours)){contours <- FALSE}
  if(contours==TRUE & is.null(contours_breaks)){contours_breaks <- round(seq(min(contours_field),max(contours_field),length.out = 10),0)}
  if(contours==TRUE & is.null(contours_spacing)){contours_spacing <- round((max(contours_field)-min(contours_field))/10,0)}
  if(contours==TRUE & is.null(contours_lty)){if(contours==TRUE){contours_lty <- 3}else{contours_lty <- 1}}
  if(contours==TRUE & is.null(contours_lwd)){if(contours==TRUE){contours_lwd <- 0.3}else{contours_lwd <- 0}}
  if(contours==TRUE & is.null(contours_labels_size)){if(contours==TRUE){contours_labels_size <- 3}else{contours_labels_size <- 0}}
  if(contours==TRUE & is.null(contours_color)){contours_color <- 'black'}
  if(contours==TRUE & is.null(contours_labels_color)){contours_labels_color <- 'white'}
  if(contours==TRUE & is.null(contours_field)){contours_field <- field}
  if(!is.null(bar_breaks) & is.null(bar_labels)){
    bar_labels <- as.character(bar_breaks)
  } else if(is.null(bar_breaks) & !is.null(bar_labels)){
    bar_breaks <- as.numeric(bar_labels)
  } else if(is.null(bar_breaks) & is.null(bar_labels)){
    bar_breaks <- waiver()
    bar_labels <- waiver()
  } 
  if(is.null(bar_lims)){bar_lims <- c(min(field),max(field))}
  
  #Create coordinates matrices if data projection is regular
  lonsc <- lon
  latsc <- lat
  if(is.null(dim(lon)) | length(dim(lon))==1){lon <- matrix(rep(lonsc, times=length(latsc)),nrow=length(lonsc),ncol=length(latsc))}
  if(is.null(dim(lat)) | length(dim(lat))==1){lat <- t(matrix(rep(latsc, times=length(lonsc)),nrow=length(latsc),ncol=length(lonsc)))}
  rm(lonsc,latsc)

  formatlonlabels <- function(vv){
    #vv should be a numeric vector of breaks coordinates
    lonformattedlabs <- sapply(vv, function(x){if(x>0){paste(as.character(x),'ºE',sep='')} else if(x<0){paste(as.character(x),'ºW',sep='')} else if(x==0){paste(as.character(x),'º',sep='')}})
    return(lonformattedlabs)  
  }
  
  formatlatlabels <- function(vv2){
    #vv2 should be a numeric vector of breaks coordinates
    latformattedlabs <- sapply(vv2, function(x){if(x>0){paste(as.character(x),'ºN',sep='')} else if(x<0){paste(as.character(x),'ºS',sep='')} else if(x==0){paste(as.character(x),'º',sep='')}})
    return(latformattedlabs)  
  }
  
  #Extract world map and crop it to make it lighter
  worldmap <- ne_coastline(scale = "medium", returnclass = "sf")
  world_cropped <- suppressMessages(suppressWarnings(
    st_crop(worldmap, xmin = xmin-30, xmax = xmax+30, ymin = ymin-4, ymax = ymax+4)))
  
  crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"
  crsLAEA <- crs
  
  #Adapt matrix data to projection
  if(crs==crsLONGLAT){
    
    df2 <- data.frame(Lon = c(lon), Lat = c(lat), FLD = c(field))
    
  } else{
    
    #Prepare data to plot with the correct projection
    df2 <- data.frame(Lon = c(lon), Lat = c(lat), FLD = c(field))
    coords <-  SpatialPoints(cbind(df2$Lon, df2$Lat), proj4string=CRS(crsLONGLAT))
    coords2 <- spTransform(coords, crsLAEA)
    df2 <- data.frame(Lon=coords2@coords[,'coords.x1'], Lat=coords2@coords[,'coords.x2'], FLD=df2$FLD)
    spdf <- SpatialPointsDataFrame(data.frame(x=df2$Lon, y=df2$Lat) , data=data.frame(z=df2$FLD))
    e <- extent(spdf)
    
    #Determine ratio between x and y dimensions
    ratio <- (e@xmax - e@xmin)/(e@ymax - e@ymin)
    dim1 <- 0.9*(dim(field)[2])
    dim2 <- floor(dim1*ratio)
    r <- raster(nrows=dim1, ncols=dim2 , ext=extent(spdf)) 
    rf <- rasterize(spdf, r, field="z", fun=mean)
    rdf <- data.frame(rasterToPoints(rf))    
    
    #Calculate the bounding box
    if(xmin!=min(lon) | xmax!=max(lon) | ymin!=min(lat) | ymax!=max(lat)){ 
      dfbounds <- data.frame(Lon = c(lon), Lat = c(lat), FLD = c(field))
      dfbounds <- dfbounds[dfbounds$Lon>=xmin & dfbounds$Lon<=xmax & dfbounds$Lat>=ymin & dfbounds$Lat<=ymax, ]
      coordsbounds <-  SpatialPoints(cbind(dfbounds$Lon, dfbounds$Lat), proj4string=CRS(crsLONGLAT))
      coords2bounds <- spTransform(coordsbounds, crsLAEA)
      dfbounds <- data.frame(Lon=coords2bounds@coords[,'coords.x1'], Lat=coords2bounds@coords[,'coords.x2'], FLD=dfbounds$FLD)
      spdfbounds <- SpatialPointsDataFrame(data.frame(x=dfbounds$Lon, y=dfbounds$Lat) , data=data.frame(z=dfbounds$FLD))
      ebounds <- extent(spdfbounds)
      
      # Determine ratio between x and y dimensions
      ratiobounds <- (ebounds@xmax - ebounds@xmin)/(ebounds@ymax - ebounds@ymin)
      rbounds <- raster(nrows=dim1, ncols=dim2 , ext=extent(spdfbounds)) 
      rfbounds <- rasterize(spdfbounds, rbounds, field="z", fun=mean)
      rdfbounds <- data.frame(rasterToPoints(rfbounds))   
      
    } else{
      
      rdfbounds <- rdf
      
    }
    
  }
  
  #Adapt contours data to projection
  if( contours==TRUE ){
    if( sum(field==contours_field)!=(dim(field)[1]*dim(field)[2]) ){ #check if contours_field is equal to field
      if(crs==crsLONGLAT){
        df2_contours <- data.frame(Lon = c(lon), Lat = c(lat), FLD = c(contours_field))
      } else{
        df_contours <- data.frame(Lon = c(lon), Lat = c(lat), FLD = c(contours_field))
        
        coords_contours <-  SpatialPoints(cbind(df_contours$Lon, df_contours$Lat), proj4string=CRS(crsLONGLAT))
        coords2_contours <- spTransform(coords_contours, crsLAEA)
        
        df_contours <- data.frame(Lon=coords2_contours@coords[,'coords.x1'], Lat=coords2_contours@coords[,'coords.x2'], FLD=df_contours$FLD)
        
        spdf_contours <- SpatialPointsDataFrame(data.frame(x=df_contours$Lon, y=df_contours$Lat) , data=data.frame(z=df_contours$FLD))
        e_contours <- extent(spdf_contours)
        
        # Determine ratio between x and y dimensions
        ratio_contours <- (e_contours@xmax - e_contours@xmin)/(e_contours@ymax - e_contours@ymin)
        r_contours <- raster(nrows=120, ncols=floor(120*ratio_contours) , ext=extent(spdf_contours))
        rf_contours <- rasterize(spdf_contours, r_contours, field="z", fun=mean)
        rdf_contours <- data.frame(rasterToPoints(rf_contours))    
      }
    } else{
      if(crs==crsLONGLAT){
        df2_contours <- df2
      } else{
        rdf_contours <- rdf
      }
      
    }
  }
  
  #Adapt arrows data to projection
  if( !is.null(arrows_field_x) & !is.null(arrows_field_y)){ #check if vector field is requested to be plotted with arrows
    
    df_arrows_x <- data.frame(Lon = c(lon), xidx = rep(1:nrow(lon),ncol(lon)), yidx=rep(1:ncol(lon),each=nrow(lon)), FLDx = c(arrows_field_x))
    df_arrows_y <- data.frame(xidx = rep(1:nrow(lon),ncol(lon)), Lat = c(lat), yidx=rep(1:ncol(lon),each=nrow(lon)), FLDy = c(arrows_field_y))
    df4 <- merge(df_arrows_x, df_arrows_y, by=c('xidx','yidx'))
    #Remove the unnecessary columns with a previous row sumsampling
    df4v2 <- dplyr::filter(df4,(xidx %% arrows_pixels_subsampling) ==0,(yidx %% arrows_pixels_subsampling) == 0)
    df4v2 <- df4v2[!names(df4v2) %in% c('xidx','yidx')]
    dgr.sub <- dplyr::mutate(df4v2,xend=Lon+(df4v2$FLDx/sd(df4v2$FLDx)^arrows_length_factor),
                             yend=Lat+(df4v2$FLDy/sd(df4v2$FLDy)^arrows_length_factor))
    dgr.sub <- dgr.sub[!names(dgr.sub) %in% c('FLDx','FLDy')]
    
    if(crs!=crsLONGLAT){
      
      dotsdf.beg <- data.frame(beglons=dgr.sub$Lon, beglats=dgr.sub$Lat)
      dotscoords.beg <-  SpatialPoints(cbind(dotsdf.beg$beglons, dotsdf.beg$beglats), proj4string=CRS(crsLONGLAT))
      dotscoords2.beg <- spTransform(dotscoords.beg, crs)
      dotsdf.beg <- data.frame(Lon=dotscoords2.beg@coords[,'coords.x1'], Lat=dotscoords2.beg@coords[,'coords.x2'])
      
      dotsdf.end <- data.frame(endlons=dgr.sub$xend, endlats=dgr.sub$yend)
      dotscoords.end <-  SpatialPoints(cbind(dotsdf.end$endlons, dotsdf.end$endlats), proj4string=CRS(crsLONGLAT))
      dotscoords2.end <- spTransform(dotscoords.end, crs)
      dotsdf.end <- data.frame(xend=dotscoords2.end@coords[,'coords.x1'], yend=dotscoords2.end@coords[,'coords.x2'])
      
      dgr.sub <- cbind(dotsdf.beg, dotsdf.end)
      
    }
  }
  
  
  #Plot  
  
  if(crs==crsLONGLAT){
    plot <- suppressWarnings(
      ggplot(data = world_cropped) +
        geom_raster(data = df2, aes(x = Lon, y = Lat, fill = FLD), show.legend=legend) + 
        geom_sf(fill='transparent', col='black') +
        coord_sf(crs = crsLONGLAT, xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE))
  } else{
    plot <- suppressWarnings(
      ggplot(data = world_cropped) +
        geom_raster(data = rdf, aes(x = x, y = y, fill = layer), show.legend=legend) +
        geom_sf(fill='transparent', col='black') +
        coord_sf(crs = crsLAEA ,xlim = c(min(rdfbounds$x), max(rdfbounds$x)), ylim = c(min(rdfbounds$y), max(rdfbounds$y)), expand = TRUE))
  }
  
  if(contours==TRUE){
    if(crs==crsLONGLAT){
      plot <- suppressWarnings(
        plot + stat_contour(data = df2_contours, aes(x = Lon, y = Lat, z = round(FLD,2)),
                            breaks=contours_breaks, binwidth = contours_spacing, #binwidth only used if breaks are not set
                            colour = contours_color, size = contours_lwd, show.legend = TRUE, 
                            linetype=contours_lty) +
          geom_text_contour(data = df2_contours, aes(x = Lon, y = Lat, z = FLD), size=contours_labels_size,
                            breaks=contours_breaks, check_overlap = TRUE, colour = contours_labels_color))
    } else{
    plot <- suppressWarnings(
      plot + stat_contour(data = rdf_contours, aes(x = x, y = y, z = round(layer,2)),
                          breaks=contours_breaks, binwidth = contours_spacing, #binwidth only used if breaks are not set
                          colour = contours_color, size = contours_lwd, show.legend = TRUE, 
                          linetype=contours_lty) +
        geom_text_contour(data = rdf_contours, aes(x = x, y = y, z = layer), size=contours_labels_size,
                          breaks=contours_breaks, check_overlap = TRUE, colour = contours_labels_color))
    }
  }
  
  if(!is.null(dots)){
    dslons <- c()
    dslats <- c()
    for(i in 1:nrow(dots)){
      dslons <- c(dslons, lon[dots[i,2],dots[i,3]])
      dslats <- c(dslats, lat[dots[i,2],dots[i,3]])
    }
    dotsdfSC <- data.frame(dslons, dslats)
    dotsdf <- data.frame(dslons, dslats)
    dotscoords <-  SpatialPoints(cbind(dotsdf$dslons, dotsdf$dslats), proj4string=CRS(crsLONGLAT))
    dotscoords2 <- spTransform(dotscoords, crsLAEA)
    
    dotsdf <- data.frame(Lon=dotscoords2@coords[,'coords.x1'], Lat=dotscoords2@coords[,'coords.x2'])
    
    plot <- suppressWarnings(
      plot + geom_point(data = dotsdf, aes(Lon, Lat), col=dots_color, pch=dots_symbol, cex=dots_cex))
    
    if(!is.null(dots_connect) && dots_connect==TRUE){
      if(nrow(dots)>1){  
        sequentdf <- c()
        for(l in 1:(nrow(dots)-1)){
          dtthreshold0 <- 0
          abcc <- FALSE
          while((dtthreshold0+1) <= dtthreshold && abcc==FALSE){
            dtthreshold0 <- dtthreshold0 + 1
            m <- l
            while(m<nrow(dots) && abcc==FALSE){
              m <- m+1
              abcc <- isSequent(as.numeric(dots[l,]),as.numeric(dots[m,]), 
                                resolution = resolution, dthreshold = dthreshold, dtthreshold = dtthreshold0)
              if(abcc==TRUE){sequentdf <- rbind(sequentdf, c(l,m))}
            }
          }
        }	    	    
        
        dotsdfsequent <- c() 
        for(l in 1:nrow(sequentdf)){dotsdfsequent <- rbind(dotsdfsequent,c(as.numeric(dotsdfSC[sequentdf[l,1],]),as.numeric(dotsdfSC[sequentdf[l,2],])))}
        dotsdfsequent <- data.frame(dotsdfsequent)
        colnames(dotsdfsequent) <- c('x','y','xend','yend')
        
        seqdf1coords <-  SpatialPoints(cbind(dotsdfsequent$x, dotsdfsequent$y), proj4string=CRS(crsLONGLAT))
        seqdf1coords2 <- spTransform(seqdf1coords, crsLAEA)
        seqdf1 <- data.frame(Lon=seqdf1coords2@coords[,'coords.x1'], Lat=seqdf1coords2@coords[,'coords.x2'])
        
        seqdf2coords <-  SpatialPoints(cbind(dotsdfsequent$xend, dotsdfsequent$yend), proj4string=CRS(crsLONGLAT))
        seqdf2coords2 <- spTransform(seqdf2coords, crsLAEA)
        seqdf2 <- data.frame(Lon=seqdf2coords2@coords[,'coords.x1'], Lat=seqdf2coords2@coords[,'coords.x2'])
        
        dotsdfsequent <- cbind(seqdf1,seqdf2)
        colnames(dotsdfsequent) <- c('x','y','xend','yend')
        
        plot <- suppressWarnings(
          plot + geom_segment(data = dotsdfsequent, aes(x=x, y=y, xend=xend, yend=yend), col=dots_connect_color))
        
      }
    }
  }
  
  if(!is.null(dots2)){
    dslons2 <- c()
    dslats2 <- c()
    for(i in 1:nrow(dots2)){
      dslons2 <- c(dslons2, lon[dots2[i,2],dots2[i,3]])
      dslats2 <- c(dslats2, lat[dots2[i,2],dots2[i,3]])
    }
    dotsdfSC2 <- data.frame(dslons2, dslats2)
    dotsdf2 <- data.frame(dslons2, dslats2)
    dotscoords22 <-  SpatialPoints(cbind(dotsdf2$dslons2, dotsdf2$dslats2), proj4string=CRS(crsLONGLAT))
    dotscoords22 <- spTransform(dotscoords22, crsLAEA)
    
    dotsdf2 <- data.frame(Lon=dotscoords22@coords[,'coords.x1'], Lat=dotscoords22@coords[,'coords.x2'])
    
    plot <- suppressWarnings(
      plot + geom_point(data = dotsdf2, aes(Lon, Lat), col=dots2_color, pch=dots2_symbol, cex=dots2_cex))
    
    if(!is.null(dots2_connect) && dots2_connect==TRUE){
      if(nrow(dots2)>1){
        sequentdf2 <- c()
        for(l2 in 1:(nrow(dots2)-1)){
          dtthreshold02 <- 0
          abcc2 <- FALSE
          while((dtthreshold02+1) <= dtthreshold2 && abcc2==FALSE){
            dtthreshold02 <- dtthreshold02 + 1
            m2 <- l2
            while(m2<nrow(dots2) && abcc2==FALSE){
              m2 <- m2+1
              abcc2 <- isSequent(as.numeric(dots2[l2,]),as.numeric(dots2[m2,]), 
                                 resolution = resolution, dthreshold = dthreshold2, dtthreshold = dtthreshold02)
              if(abcc2==TRUE){sequentdf2 <- rbind(sequentdf2, c(l2,m2))}
            }
          }
        }	    	    
        
        dotsdfsequent2 <- c() 
        for(l2 in 1:nrow(sequentdf2)){dotsdfsequent2 <- rbind(dotsdfsequent2,c(as.numeric(dotsdfSC2[sequentdf2[l2,1],]),as.numeric(dotsdfSC2[sequentdf2[l2,2],])))}
        dotsdfsequent2 <- data.frame(dotsdfsequent2)
        colnames(dotsdfsequent2) <- c('x','y','xend','yend')
        
        seqdf1coords2 <-  SpatialPoints(cbind(dotsdfsequent2$x, dotsdfsequent2$y), proj4string=CRS(crsLONGLAT))
        seqdf1coords22 <- spTransform(seqdf1coords2, crsLAEA)
        seqdf12 <- data.frame(Lon=seqdf1coords22@coords[,'coords.x1'], Lat=seqdf1coords22@coords[,'coords.x2'])
        
        seqdf2coords2 <-  SpatialPoints(cbind(dotsdfsequent2$xend, dotsdfsequent2$yend), proj4string=CRS(crsLONGLAT))
        seqdf2coords22 <- spTransform(seqdf2coords2, crsLAEA)
        seqdf22 <- data.frame(Lon=seqdf2coords22@coords[,'coords.x1'], Lat=seqdf2coords22@coords[,'coords.x2'])
        
        dotsdfsequent2 <- cbind(seqdf12,seqdf22)
        colnames(dotsdfsequent2) <- c('x','y','xend','yend')
        
        plot <- suppressWarnings(
          plot + geom_segment(data = dotsdfsequent2, aes(x=x, y=y, xend=xend, yend=yend), col=dots2_connect_color))
      }
    }
  }
  
  if(!is.null(dots3)){
    dslons3 <- c()
    dslats3 <- c()
    for(i in 1:nrow(dots3)){
      dslons3 <- c(dslons3, lon[dots3[i,2],dots3[i,3]])
      dslats3 <- c(dslats3, lat[dots3[i,2],dots3[i,3]])
    }
    dotsdfSC3 <- data.frame(dslons3, dslats3)
    dotsdf3 <- data.frame(dslons3, dslats3)
    dotscoords33 <-  SpatialPoints(cbind(dotsdf3$dslons3, dotsdf3$dslats3), proj4string=CRS(crsLONGLAT))
    dotscoords33 <- spTransform(dotscoords33, crsLAEA)
    
    dotsdf3 <- data.frame(Lon=dotscoords33@coords[,'coords.x1'], Lat=dotscoords33@coords[,'coords.x2'])
    
    plot <- suppressWarnings(
      plot + geom_point(data = dotsdf3, aes(Lon, Lat), col=dots3_color, pch=dots3_symbol, cex=dots3_cex))
    
    if(!is.null(dots3_connect) && dots3_connect==TRUE){
      if(nrow(dots3)>1){
        sequentdf3 <- c()
        for(l3 in 1:(nrow(dots3)-1)){
          dtthreshold03 <- 0
          abcc3 <- FALSE
          while((dtthreshold03+1) <= dtthreshold3 && abcc3==FALSE){
            dtthreshold03 <- dtthreshold03 + 1
            m3 <- l3
            while(m3<nrow(dots3) && abcc3==FALSE){
              m3 <- m3+1
              abcc3 <- isSequent(as.numeric(dots3[l3,]),as.numeric(dots3[m3,]), 
                                 resolution = resolution, dthreshold = dthreshold3, dtthreshold = dtthreshold03)
              if(abcc3==TRUE){sequentdf3 <- rbind(sequentdf3, c(l3,m3))}
            }
          }
        }	    	    
        
        dotsdfsequent3 <- c() 
        for(l3 in 1:nrow(sequentdf3)){dotsdfsequent3 <- rbind(dotsdfsequent3,c(as.numeric(dotsdfSC3[sequentdf3[l3,1],]),as.numeric(dotsdfSC3[sequentdf3[l3,2],])))}
        dotsdfsequent3 <- data.frame(dotsdfsequent3)
        colnames(dotsdfsequent3) <- c('x','y','xend','yend')
        
        seqdf1coords3 <-  SpatialPoints(cbind(dotsdfsequent3$x, dotsdfsequent3$y), proj4string=CRS(crsLONGLAT))
        seqdf1coords33 <- spTransform(seqdf1coords3, crsLAEA)
        seqdf13 <- data.frame(Lon=seqdf1coords33@coords[,'coords.x1'], Lat=seqdf1coords33@coords[,'coords.x2'])
        
        seqdf2coords3 <-  SpatialPoints(cbind(dotsdfsequent3$xend, dotsdfsequent3$yend), proj4string=CRS(crsLONGLAT))
        seqdf2coords33 <- spTransform(seqdf2coords3, crsLAEA)
        seqdf33 <- data.frame(Lon=seqdf2coords33@coords[,'coords.x1'], Lat=seqdf2coords33@coords[,'coords.x2'])
        
        dotsdfsequent3 <- cbind(seqdf13,seqdf33)
        colnames(dotsdfsequent3) <- c('x','y','xend','yend')
        
        plot <- suppressWarnings(
          plot + geom_segment(data = dotsdfsequent3, aes(x=x, y=y, xend=xend, yend=yend), col=dots3_connect_color))
      }
    }
  }
  
  if(!is.null(arrows_field_x) & !is.null(arrows_field_y)){
    plot <- suppressWarnings(
      plot + 
        geom_segment(data=dgr.sub,aes(x=Lon,y=Lat,xend=xend,yend=yend),arrow = arrow(length = unit(0.01, "npc")),
                     col="black",alpha=arrows_alpha))
  }
  
  plot <- suppressWarnings( 
    plot +
      guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", 
                                   title = field_units, title.position = "top", 
                                   frame.linewidth = 1.5, frame.linetype = 1, 
                                   ticks = TRUE, barheight = 15, ticks.linewidth = 1,
                                   barwidth = 1, title.vjust = 3)) +   
      scale_fill_gradientn(colours = bar_color_scale,
                           name = field_units,na.value="transparent", 
                           limits = c(bar_lims[1],bar_lims[2]), 
                           oob = scales::squish,
                           breaks = bar_breaks,
                           labels = bar_labels,
                           space = 2.5, guide = "colourbar") +
      scale_x_continuous(breaks=seq(round(min(lon),0),round(max(lon),0),by=x_ticks_interval),
                         labels = formatlonlabels) +
      scale_y_continuous(breaks=seq(round(min(lat),0),round(max(lat),0),by=y_ticks_interval), 
                         labels = formatlatlabels) +
      annotation_scale(location = "br", width_hint = 0.1, bar_cols = c("black", "white"), 
                       line_width = 1, height = unit(0.1, "cm"), text_pad = unit(0.15, "cm"), 
                       tick_height = 0.5) +
      theme_bw() + 
      theme(axis.title = element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.text.x.bottom = element_text(colour = " black", size = 23, face = "bold", vjust = 0.5),
            axis.text.y.left = element_text(colour = " black", size = 23, angle = 0, face = "bold"),
            axis.text.y = element_text(vjust=0.5),
            legend.title = element_text(colour = " black", size = 11, face = "bold"),
            legend.text = element_text(size = 10, colour = "black", face = "bold"),
            legend.margin = margin(6, 6, 6, 6),
            legend.text.align = 1,
            legend.box.spacing = unit(3, "lines"),
            legend.key.size = unit(1, "lines"),
            line = element_line(colour = "black"),
            plot.title = element_text(face = "bold", size = (15),
                                      hjust = 0.5, vjust = 8)) +
      
      labs(title = main_title))
  
  suppressMessages(print(plot))
  
}
