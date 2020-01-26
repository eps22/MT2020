CheckZGradSymm <- function(compfolder, consecutivecenters, currentimeidx, resolution, timestep, ZGradSymmPressureLevel, ZGradSymmNumberOfDirections, ZGradSymmMeasure, ZGradSymmthreshold, zerovortradius){
  
  if(ZGradSymmNumberOfDirections>8){ZGradSymmNumberOfDirections <- 8}
  
  zfieldnc <- paste(compfolder,'/output/outputfile-z.nc',sep='')
  
  zheight <- nc_open(zfieldnc)
  zz <- ncvar_get(zheight,'height',start = c(1,1,1,currentimeidx),count=c(-1,-1,-1,1))
  
  #Get requested plevel index
  plevs <- zheight$dim$interp_level$vals
  hpaPLEVELidx <- which(plevs==ZGradSymmPressureLevel, arr.ind = TRUE)
  
  if(length(hpaPLEVELidx)==0){stop('Error: selected pressure level for checking Z symmetry is not a valid pressure level.')}

  zPLEVEL <- zz[,,hpaPLEVELidx] 
  currentcenter <- consecutivecenters[1,]
  nextcenter <- consecutivecenters[2,]
  movementdirectionx <- as.numeric(nextcenter[1] - currentcenter[1])
  movementdirectiony <- as.numeric(nextcenter[2] - currentcenter[2])
  movementdirectionMOD <- sqrt(movementdirectionx^2 + movementdirectiony^2)
  
  if(movementdirectionMOD==0){
    isZPLEVELSymm <- FALSE
  } else{
    
    INradiuspoints <- which((dd(zPLEVEL,currentcenter)*resolution)<zerovortradius,arr.ind=TRUE)

    angles <- c()
    for(i in 1:nrow(INradiuspoints)){
      pointx <- as.numeric(INradiuspoints[i,1])
      pointy <- as.numeric(INradiuspoints[i,2])
      pointMOD <- sqrt(pointx^2 + pointy^2)
      
      pointtocentervecx <- pointx-currentcenter[1]
      pointtocentervecy <- pointy-currentcenter[2]
      pointtocentervecMOD <- sqrt(pointtocentervecx^2 + pointtocentervecy^2)
      if(pointtocentervecMOD != 0){
        dotproduct <- pointtocentervecx*movementdirectionx + pointtocentervecy*movementdirectiony
        cosangle <- dotproduct/(pointtocentervecMOD*movementdirectionMOD)
        zcross <- movementdirectionx*pointtocentervecy - pointtocentervecx*movementdirectiony
        if(zcross < 0){angle <- 2*pi - acos(cosangle)} 
        if(zcross >= 0){angle <- acos(cosangle)}
      } else{angle <- NA}
      angles <- c(angles, angle)
    }
    anglesdeg <- angles*360/(2*pi)
    INradiuspoints <- INradiuspoints[!is.na(anglesdeg),]
    anglesdeg <- as.numeric(anglesdeg[!is.na(anglesdeg)])

    in0points <- INradiuspoints[which(abs(anglesdeg-0)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in0points)==0){max0d <- NA; whichmax0d <- c(NA, NA)} else{
    max0d <- max(dd(in0points,currentcenter))*resolution
    whichmax0d <- in0points[which.max(dd(in0points,currentcenter)),]
    }
    in45points <- INradiuspoints[which(abs(anglesdeg-45)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in45points)==0){max45d <- NA; whichmax45d <- c(NA, NA)} else{
    max45d <- max(dd(in45points,currentcenter))*resolution
    whichmax45d <- in45points[which.max(dd(in45points,currentcenter)),]
    }
    in90points <- INradiuspoints[which(abs(anglesdeg-90)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in90points)==0){max90d <- NA; whichmax90d <- c(NA, NA)} else{
    max90d <- max(dd(in90points,currentcenter))*resolution
    whichmax90d <- in90points[which.max(dd(in90points,currentcenter)),]
    }
    in135points <- INradiuspoints[which(abs(anglesdeg-135)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in135points)==0){max135d <- NA; whichmax135d <- c(NA, NA)} else{
    max135d <- max(dd(in135points,currentcenter))*resolution
    whichmax135d <- in135points[which.max(dd(in135points,currentcenter)),]
    }
    in180points <- INradiuspoints[which(abs(anglesdeg-180)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in180points)==0){max180d <- NA; whichmax180d <- c(NA, NA)} else{
    max180d <- max(dd(in180points,currentcenter))*resolution
    whichmax180d <- in180points[which.max(dd(in180points,currentcenter)),]
    }
    in225points <- INradiuspoints[which(abs(anglesdeg-225)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in225points)==0){max225d <- NA; whichmax225d <- c(NA, NA)} else{
    max225d <- max(dd(in225points,currentcenter))*resolution
    whichmax225d <- in225points[which.max(dd(in225points,currentcenter)),]
    }
    in270points <- INradiuspoints[which(abs(anglesdeg-270)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in270points)==0){max270d <- NA; whichmax270d <- c(NA, NA)} else{
    max270d <- max(dd(in270points,currentcenter))*resolution
    whichmax270d <- in270points[which.max(dd(in270points,currentcenter)),]
    }
    in315points <- INradiuspoints[which(abs(anglesdeg-315)<2,arr.ind = TRUE),,drop=FALSE]
    if(nrow(in315points)==0){max315d <- NA; whichmax315d <- c(NA, NA)} else{
    max315d <- max(dd(in315points,currentcenter))*resolution
    whichmax315d <- in315points[which.max(dd(in315points,currentcenter)),]
    }
    
    if(sum(is.na(whichmax0d))==2){grad0 <- NA} else{
    grad0 <- (zPLEVEL[whichmax0d[1],whichmax0d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max0d
    }
    if(sum(is.na(whichmax45d))==2){grad45 <- NA} else{
    grad45 <- (zPLEVEL[whichmax45d[1],whichmax45d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max45d
    }
    if(sum(is.na(whichmax90d))==2){grad90 <- NA} else{
    grad90 <- (zPLEVEL[whichmax90d[1],whichmax90d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max90d
    }
    if(sum(is.na(whichmax135d))==2){grad135 <- NA} else{
    grad135 <- (zPLEVEL[whichmax135d[1],whichmax135d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max135d
    }
    if(sum(is.na(whichmax180d))==2){grad180 <- NA} else{
    grad180 <- (zPLEVEL[whichmax180d[1],whichmax180d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max180d
    }
    if(sum(is.na(whichmax225d))==2){grad225 <- NA} else{
    grad225 <- (zPLEVEL[whichmax225d[1],whichmax225d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max225d
    }
    if(sum(is.na(whichmax270d))==2){grad270 <- NA} else{
    grad270 <- (zPLEVEL[whichmax270d[1],whichmax270d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max270d
    }
    if(sum(is.na(whichmax135d))==2){grad315 <- NA} else{
    grad315 <- (zPLEVEL[whichmax315d[1],whichmax315d[2]]-zPLEVEL[currentcenter[1],currentcenter[2]])/max315d
    }
    
    vectofgrads <- c(grad0,grad45,grad90,grad135,grad180,grad225,grad270,grad315)
    vectofgrads <- vectofgrads[1:ZGradSymmNumberOfDirections]
    if(sum(is.na(vectofgrads))==length(vectofgrads)){isZPLEVELSymm <- FALSE} else{
      if(ZGradSymmMeasure=='mean'){isZPLEVELSymm <- mean(vectofgrads,na.rm = TRUE)>ZGradSymmthreshold}
      if(ZGradSymmMeasure=='min'){isZPLEVELSymm <- min(vectofgrads, na.rm = TRUE)>ZGradSymmthreshold}
    }
     
  }
  
  return(isZPLEVELSymm)
    
}
