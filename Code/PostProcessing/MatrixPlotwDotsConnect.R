MatrixPlotwDotsConnect <- function(field, lon, lat, crs, resolution = NULL, dthreshold = NULL, dtthreshold = NULL,
                                   xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, contours = NULL,
                                   contours_field = NULL , contours_spacing = NULL, contours_breaks = NULL, 
                                   contours_color = NULL , contours_labels_size = NULL, contours_labels_color = NULL,
                                   contours_lty = NULL, contours_lwd = NULL, bar_lims = NULL, bar_breaks = NULL, bar_labels = NULL, 
                                   x_ticks_interval = NULL, y_ticks_interval = NULL, legend = TRUE,
                                   bar_color_scale = NULL, field_units = NULL, main_title = NULL, 
                                   dots = NULL, dots_symbol = NULL, dots_cex = NULL, dots_color = NULL, dots_connect = NULL, dots_connect_color = NULL){
  
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
  
  source('./PostProcessing/isSequent.R')
  
  if(is.null(xmin)){xmin <- min(lon)}
  if(is.null(xmax)){xmax <- max(lon)}
  if(is.null(ymin)){ymin <- min(lat)}
  if(is.null(ymax)){ymax <- max(lat)}
  if(is.null(bar_color_scale)){bar_color_scale <- brewer.pal(11, "RdBu")}
  if(is.null(resolution)){resolution <- 9}
  if(is.null(dthreshold)){dthreshold <- 100}
  if(is.null(dtthreshold)){dtthreshold <- 1}
  if(is.null(field_units)){field_units <- ''}
  if(is.null(x_ticks_interval)){x_ticks_interval <- 4}
  if(is.null(y_ticks_interval)){y_ticks_interval <- 4}
  if(is.null(main_title)){main_title <- ''}
  if(is.null(dots_symbol)){dots_symbol <- 3}
  if(is.null(dots_cex)){dots_cex <- 4}
  if(is.null(dots_color)){dots_color <- 'black'}
  if(is.null(dots_connect_color)){dots_connect_color <- 'black'}
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
  
  lonsc <- lon
  latsc <- lat
  #Create coordinates matrices if data projection is regular
  if(is.null(dim(lon))){lon <- matrix(rep(lonsc, times=length(latsc)),nrow=length(lonsc),ncol=length(latsc))}
  if(is.null(dim(lat))){lat <- t(matrix(rep(latsc, times=length(latsc)),nrow=length(latsc),ncol=length(lonsc)))}
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
  
  #bounding box polygon in long/lat projection, i.e. axis-aligned
  bb <- st_sfc(
    st_polygon(list(cbind(
      c(xmin, xmax, xmax, xmin, xmin), # x-coordinates (longitudes) of points A,B,C,D
      c(ymin, ymin, ymax, ymax, ymin)   # y-coordinates (latitudes) of points A,B,C,D
    ))),
    crs = crsLONGLAT)
  
  # now in LAEA projection
  laeabb <- st_transform(bb, crs = crsLAEA)
  
  # the extent of the bounding box in the new projection
  b <- st_bbox(laeabb)
  
  #Prepare data to plot with the correct projection
  df2 <- data.frame(Lon = c(lon), Lat = c(lat), FLD = c(field))
  
  coords <-  SpatialPoints(cbind(df2$Lon, df2$Lat), proj4string=CRS(crsLONGLAT))
  coords2 <- spTransform(coords, crsLAEA)
  
  df2 <- data.frame(Lon=coords2@coords[,'coords.x1'], Lat=coords2@coords[,'coords.x2'], FLD=df2$FLD)
  
  spdf <- SpatialPointsDataFrame(data.frame(x=df2$Lon, y=df2$Lat) , data=data.frame(z=df2$FLD))
  e <- extent(spdf)
  
  # Determine ratio between x and y dimensions
  ratio <- (e@xmax - e@xmin)/(e@ymax - e@ymin)
  r <- raster(nrows=180, ncols=floor(180*ratio) , ext=extent(spdf)) #MODIFICAR PARA QUE SEA AUTOMATICO
  rf <- rasterize(spdf, r, field="z", fun=mean)
  
  rdf <- data.frame(rasterToPoints(rf))    
  
  #Adapt contours data to projection
  if( contours==TRUE ){
    if( sum(field==contours_field)!=(dim(field)[1]*dim(field)[2]) ){ #check if contours_field is equal to field
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
    } else{
      rdf_contours <- rdf
    }
  }
  
  plot <- suppressWarnings(
    ggplot(data = world_cropped) +
      geom_raster(data = rdf, aes(x = x, y = y, fill = layer), show.legend=legend) +
      geom_sf(fill='transparent', col='black') +
      coord_sf(crs = crsLAEA ,xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]), expand = FALSE))
  
  
  if(contours==TRUE){
    plot <- suppressWarnings(
      plot + stat_contour(data = rdf_contours, aes(x = x, y = y, z = round(layer,2)),
                          breaks=contours_breaks, binwidth = contours_spacing, #binwidth only used if breaks are not set
                          colour = contours_color, size = contours_lwd, show.legend = TRUE, 
                          linetype=contours_lty) +
        geom_text_contour(data = rdf_contours, aes(x = x, y = y, z = layer), size=contours_labels_size,
                          breaks=contours_breaks, check_overlap = TRUE, colour = contours_labels_color))
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

