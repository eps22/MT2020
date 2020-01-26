ThermalSymmB <- function(compfolder, currentimeidx, consecutivecenters, resolution, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, zerovortradius2){
  
    zfieldnc <- paste(compfolder,'/output/outputfile-z.nc',sep='')
 
    plevlow <- Blowerpressurelevel
    plevhigh <- Bupperpressurelevel
    if(plevlow<plevhigh){cat('Selected lower pressure level is above the upper one and may lead to inconsistent medicanes detections')}
    h=+1

    library(ncdf4)

    zheight <- nc_open(zfieldnc)
    zz <- ncvar_get(zheight,'height',start = c(1,1,1,currentimeidx),count=c(-1,-1,-1,1))
    
    #Get requested plevels indices
    plevs <- zheight$dim$interp_level$vals
    hpaHIGHidx <- which(plevs==plevlow, arr.ind = TRUE)
    hpaLOWidx <- which(plevs==plevhigh, arr.ind = TRUE)
    
    if(length(hpaHIGHidx)==0){stop('Error: selected lower pressure level for the calculation of thermal symmetry is not a valid pressure level.')}
    if(length(hpaLOWidx)==0){stop('Error: selected higher pressure level for the calculation of thermal symmetry is not a valid pressure level.')}
    
    z900 <- zz[,,hpaHIGHidx]
    z600 <- zz[,,hpaLOWidx] 
    z600m900 <- z600-z900    
    z500 <- zz
  
    currentcenter <- consecutivecenters[1,]
    nextcenter <- consecutivecenters[2,]
    movementdirectionx <- as.numeric(nextcenter[1] - currentcenter[1])
    movementdirectiony <- as.numeric(nextcenter[2] - currentcenter[2])
    movementdirectionMOD <- sqrt(movementdirectionx^2 + movementdirectiony^2)
  
    if(movementdirectionMOD==0){
      B <- NA
    } else{
        INradiuspoints <- which((dd(z900,currentcenter)*resolution)<zerovortradius2,arr.ind=TRUE) #c(155,86)
        angles <- c()
        rightpoints <- c()
        leftpoints <- c()
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
      
       #Divide by motion direction
       pointsrightidxs <- which(anglesdeg>180,arr.ind = TRUE) #points right to the direction of motion
       pointsleftidxs <- which(anglesdeg<180  & anglesdeg!=0,arr.ind = TRUE) #points left to the direction of motion
       rightpoints <- INradiuspoints[pointsrightidxs,]
       leftpoints <- INradiuspoints[pointsleftidxs,]
     
       z600m900right <- mean(z600m900[rightpoints])
       z600m900left <- mean(z600m900[leftpoints])
        
       B0=h*(z600m900right - z600m900left)
       
       #Divide by a 45º rotated line respect to the motion direction
       pointsrightidxs45 <- which(anglesdeg>135 & anglesdeg<315,arr.ind = TRUE) #points left to the direction of motion      
       pointsleftidxs45 <- c(which(anglesdeg<135  & anglesdeg!=0,arr.ind = TRUE),which(anglesdeg>315,arr.ind = TRUE )) #points right to the direction of motion
       rightpoints45 <- INradiuspoints[pointsrightidxs45,]
       leftpoints45 <- INradiuspoints[pointsleftidxs45,]
       
       z600m900right45 <- mean(z600m900[rightpoints45])
       z600m900left45 <- mean(z600m900[leftpoints45])
       
       B45=h*(z600m900right45 - z600m900left45)
      
       #Divide by the perpendicular to the motion direction
       pointsleftidxsP <- c(which(anglesdeg>270,arr.ind = TRUE),which(anglesdeg<90 & anglesdeg!=0,arr.ind = TRUE )) #points right to the direction of motion
       pointsrightidxsP <- which(anglesdeg>90 & anglesdeg<270,arr.ind = TRUE) #points left to the direction of motion
       rightpointsP <- INradiuspoints[pointsrightidxsP,]
       leftpointsP <- INradiuspoints[pointsleftidxsP,]

       z600m900rightP <- mean(z600m900[rightpointsP])
       z600m900leftP <- mean(z600m900[leftpointsP])
       
       BP=h*(z600m900rightP - z600m900leftP)
        
       #Divide by a -45º rotated line respect to the motion direction
       pointsrightidxsm45 <- which(anglesdeg>45 & anglesdeg<225,arr.ind = TRUE) #points left to the direction of motion
       pointsleftidxsm45 <- c(which(anglesdeg<45  & anglesdeg!=0,arr.ind = TRUE),which(anglesdeg>225,arr.ind = TRUE )) #points right to the direction of motion
       rightpointsm45 <- INradiuspoints[pointsrightidxsm45,]
       leftpointsm45 <- INradiuspoints[pointsleftidxsm45,]
      
       z600m900rightm45 <- mean(z600m900[rightpointsm45])
       z600m900leftm45 <- mean(z600m900[leftpointsm45])
      
      Bm45=h*(z600m900rightm45 - z600m900leftm45)
      
      if(Bmultiplemeasure=='mean'){
        if(Bdirections==1){
          B <- abs(B0)
        } else if(Bdirections==2){
          B <- mean(c(abs(B0), abs(BP))[!is.na(c(abs(B0), abs(BP)))])
        } else if(Bdirections==3){
          B <- mean(c(abs(B0), abs(BP), abs(B45))[!is.na(c(abs(B0), abs(BP), abs(B45)))])
        } else if(Bdirections==4){
          B <- mean(c(abs(B0), abs(BP), abs(B45), abs(Bm45))[!is.na(c(abs(B0), abs(BP), abs(B45), abs(Bm45)))])
        }
      } else if(Bmultiplemeasure=='max'){
        if(Bdirections==1){
          B <- abs(B0)
        } else if(Bdirections==2){
          B <- max(c(abs(B0), abs(BP))[!is.na(c(abs(B0), abs(BP)))])
        } else if(Bdirections==3){
          B <- max(c(abs(B0), abs(BP), abs(B45))[!is.na(c(abs(B0), abs(BP), abs(B45)))])
        } else if(Bdirections==4){
          B <- max(c(abs(B0), abs(BP), abs(B45), abs(Bm45))[!is.na(c(abs(B0), abs(BP), abs(B45), abs(Bm45)))])
        }
      }
    
  }

  return(B)
  
  
} 
