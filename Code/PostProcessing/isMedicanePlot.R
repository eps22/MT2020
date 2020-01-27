#PlotComplete <- function(directory1,directory2, dot.colors, dot.sizes){


### EDITABLE PART

directory <- '/Datos/FindMedicanes-V13'

#Event <- 'Rolf'
#Event <- 'Numa'
#Event <- 'Cornelia'
Event <- 'Celeno'
#Event <- 'Qendresa'
#Event <- 'Leucosia'

same <- FALSE

fix <- 'chem' #'chem' or 'nud'

chemnochem <- 'WRF' #'WRF' or 'WRF-chem'. Only used when fix='chem'
nonudspnud <- 'nonud' #'nonud' or 'spnud'. Only used when fix='nud'

exclude <- c()
#include <- c(1:9)

scale_lower_limit <- 995 #NULL #For a default handling of the scale limits, please set to NULL  
scale_upper_limit <- 1010 #NULL #For a default handling of the scale limits, please set to NULL  

### END OF EDITABLE PART





{
  #colors: one of {'categoryK', 'minslp', 'energy'}
  #sizes: one of {numeric, 'categoryK', 'minslp', 'energy'}
  dot.colors <- 'minslp'
  dot.sizes <- 2
  
  inbetweenfactor <- 8
  displacement <- -2.0
  
  if(fix=='nud'){
    dir1identifier <- 'no aerosols'
    dir2identifier <- 'aerosols'
  } else{
    dir1identifier <- 'no nudging'
    dir2identifier <- 'spectral nudging'
  }
  
  library(varhandle)
  library(ncdf4)
  library(s2dverification)
  library(RColorBrewer)
  Spectral <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))
  
  if(fix=='chem'){
    if(chemnochem=='WRF'){
      chemnochem2 <- '-'
    } else{chemnochem2 <- '-chem-'}
    if(same==TRUE){
      dirs1 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,chemnochem2,nonudspnud,sep=""),list.files(paste(directory,Event,sep='/')))]
      dirs2 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,chemnochem2,nonudspnud,sep=""),list.files(paste(directory,Event,sep='/')))]
    } else{
      dirs1 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,chemnochem2,'nonud',sep=""),list.files(paste(directory,Event,sep='/')))]
      dirs2 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,chemnochem2,'spnud',sep=""),list.files(paste(directory,Event,sep='/')))]
    }
  } 
  if(fix=='nud'){
    if(chemnochem=='WRF'){
      chemnochem2 <- '-'
    } else{chemnochem2 <- '-chem-'}
    if(same==TRUE){
      dirs1 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,chemnochem2,nonudspnud,sep=""),list.files(paste(directory,Event,sep='/')))]
      dirs2 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,chemnochem2,nonudspnud,sep=""),list.files(paste(directory,Event,sep='/')))]
    } else{
      dirs1 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,nonudspnud,sep="-"),list.files(paste(directory,Event,sep='/')))]
      dirs2 <- list.files(paste(directory,Event,sep='/'))[grep(paste(Event,'chem',nonudspnud,sep="-"),list.files(paste(directory,Event,sep='/')))]
    }
  }
  
  if(length(exclude)>0){
    dirs1 <- dirs1[-exclude]
    dirs2 <- dirs2[-exclude]
  }
  
  if(Event=='Leucosia'){
    medicaneyear <- '1982'              #Leucosia
    respecttoday <-  26                 #Leucosia
    respecttohour <- 00                 #Leucosia
    respectdate <- '1982-01-26 00h'     #Leucosia
    initialdate <- '1982-01-23 12:00:00' #Leucosia
    finaldate <- '1982-01-28 00:00:00'  #Leucosia
  }
  if(Event=='Celeno'){
    medicaneyear <- '1995'              #Celeno
    respecttoday <-  15                 #Celeno
    respecttohour <- 00                 #Celeno
    respectdate <- '1995-01-14 00h'     #Celeno
    initialdate <- '1995-01-13 12:00:00' #Celeno
    finaldate <- '1995-01-16 22:00:00'  #Celeno
  }
  if(Event=='Cornelia'){
    medicaneyear <- '1996'              #Cornelia
    respecttoday <- 09                  #Cornelia
    respecttohour <- 00                 #Cornelia
    respectdate <- '1996-10-09 00h'     #Cornelia
    initialdate <- '1996-10-06 12:00:00' #Cornelia
    finaldate <- '1996-10-10 12:00:00'  #Cornelia
  }
  #if(Event=='Cornelia'){
  #  medicaneyear <- '1996'                #EarlyCornelia
  #  respecttoday <- 05                    #EarlyCornelia
  #  respecttohour <- 00                   #EarlyCornelia
  #  respectdate <- '1996-10-05 00h'       #EarlyCornelia
  #  initialdate <- '1996-10-02 00:00:00'  #EarlyCornelia
  #  finaldate <- '1996-10-07 00:00:00'    #EarlyCornelia
  #}
  if(Event=='Rolf'){
    medicaneyear <- '2011'                #Rolf
    respecttoday <-  08                   #Rolf
    respecttohour <- 00                   #Rolf
    respectdate <- '2011-11-08 00h'       #Rolf
    initialdate <- '2011-11-05 12:00:00'
    finaldate <- '2011-11-10 00:00:00'    #Rolf
  }  
  if(Event=='Qendresa'){
    medicaneyear <- '2014'              #Qendresa
    respecttoday <-  08                 #Qendresa
    respecttohour <- 00                 #Qendresa
    respectdate <- '2014-11-07 12h'     #Qendresa
    initialdate <- '2014-11-06 00:00:00'  #Qendresa
    finaldate <- '2014-11-09 00:00:00'  #Qendresa
  }
  if(Event=='Numa'){
    medicaneyear <- '2017'              #Numa
    respecttoday <-  18                 #Numa
    respecttohour <- 00                 #Numa
    respectdate <- '2017-11-18 00h'
    initialdate <- '2017-11-13 12:00:00'
    finaldate <- '2017-11-20 00:00:00'  #Numa
  }
  
  
  ismedicanedf <- c()
  slpmindf <- c()
  rut <- c()
  
  dir <- dirs1[1]
  file <- paste(directory,'/',dir,'/',dir,'-trackdf.csv',sep="")
  ncinslp <- nc_open(paste(directory,'/',Event,'/',dir,'/output/outputfile-slp.nc',sep=""))
  absmaxtime <- ncinslp$dim$Time$len
  numofdates <- ncinslp$dim$Time$vals + 1
  initdate <- strsplit((strsplit(ncinslp$dim$Time$units, ' ')[[1]])[3], '-')[[1]]
  inithour <- strsplit((strsplit(ncinslp$dim$Time$units, ' ')[[1]])[4],':')[[1]]
  init <- c(initdate, inithour)
  alldates <- seq.POSIXt(do.call(ISOdate,as.list(init)), by="hour", length.out = max(numofdates))
  
  for(dir in dirs1){
    
    file <- paste(directory,'/',Event,'/',dir,'/',dir,'-trackingdf.csv',sep="")
    maxtime <- nc_open(paste(directory,'/',Event,'/',dir,'/output/outputfile-slp.nc',sep=""))$dim$Time$len
    
    if(file.exists(file)){
      data0 <- read.csv(file,header=TRUE,row.names = 1)
      data <- c()
      for(i in 1:maxtime){
        if(i %in% data0$timestep){
          selectindex <- which(data0$timestep==i,arr.ind = TRUE)
          if(length(selectindex)==1){
            data <- rbind(data, as.numeric(data0[selectindex,c('timestep','x','y','center.slp.hPa')]))
          } else if(length(selectindex)>1){
            if(nrow(data)>0){
              candidatesslps <- as.numeric(data0[selectindex,'center.slp.hPa'])
              previousdatarowslp <- tail(data[!is.na(data[,4]),],n=1)[4]
              selectindex <- selectindex[which.min(abs(candidatesslps-previousdatarowslp))]
              data <- rbind(data, as.numeric(data0[selectindex,c('timestep','x','y','center.slp.hPa')]))
            } else{
              data <- rbind(data,c(i, NA, NA, NA))
            }
          } else{
            data <- rbind(data,c(i, NA, NA, NA))
          }
        } else{
          data <- rbind(data,c(i, NA, NA, NA))
        }
      }
      data <- data.frame(data)
      colnames(data) <- c('timestep','x','y','center.slp.hPa')
    } else{data <- data.frame(timestep=1:maxtime,x=rep(NA,maxtime),y=rep(NA,maxtime),center.slp.hPa=rep(NA,maxtime))}
    
    #Is medicane or not -by means of the Hart parameters- for the pch attribute
    ismedicane <- !is.na(data$center.slp.hPa)
    
    if(!is.null(ismedicanedf)){
      if(length(ismedicane) < nrow(ismedicanedf)){
        ismedicane <- c(rep(NA,(nrow(ismedicanedf)-length(ismedicane))), ismedicane)
      }
    }
    ismedicanedf <- cbind(ismedicanedf, ismedicane)
    
    #Slp minimum values for the color attribute
    
    slpmin <- data$center.slp.hPa
    if(!is.null(slpmindf)){
      if(length(slpmin) < nrow(slpmindf)){
        slpmin <- c(rep(NA,(nrow(slpmindf)-length(slpmin))), slpmin)
      }
    }
    slpmindf <- cbind(slpmindf, slpmin)
    
    rut <- c(rut, absmaxtime-maxtime)
    
  }
  
  casenames <- sapply(dirs1, function(x) strsplit(strsplit((strsplit(x,'to')[[1]][1]),paste('-',medicaneyear,'-',sep=''))[[1]][2],'-')[[1]][2], USE.NAMES=FALSE)
  day <- sapply(casenames, function(x) strsplit(x,'_')[[1]][1], USE.NAMES = FALSE)
  hour <- sapply(casenames, function(x) strsplit(x,'_')[[1]][2], USE.NAMES = FALSE)
  colnames(ismedicanedf) <- casenames
  colnames(slpmindf) <- casenames
  
  
  datafr <- data.frame(x=1:nrow(ismedicanedf))
  for(cc in 1:ncol(ismedicanedf)){
    vec <- ismedicanedf[,cc]*(cc)
    vec[vec==0] <- NA
    datafr <- cbind(datafr,vec)
  }
  colnames(datafr) <- c('x',paste('v',1:ncol(ismedicanedf),sep=""))
  
  
  ismedicanedf2 <- c()
  slpmindf2 <- c()
  rut2 <- c()
  
  dir <- dirs2[1]
  file <- paste(directory,'/',Event,'/',dir,'/',dir,'-trackingdf.csv',sep="")
  absmaxtime <- nc_open(paste(directory,'/',Event,'/',dir,'/output/outputfile-slp.nc',sep=""))$dim$Time$len
  
  for(dir in dirs2){
    
    file <- paste(directory,'/',Event,'/',dir,'/',dir,'-trackingdf.csv',sep="")
    maxtime <- nc_open(paste(directory,'/',Event,'/',dir,'/output/outputfile-slp.nc',sep=""))$dim$Time$len
    
    if(file.exists(file)){
      data0 <- read.csv(file,header=TRUE,row.names = 1)
      data <- c()
      for(i in 1:maxtime){
        if(i %in% data0$timestep){
          selectindex <- which(data0$timestep==i,arr.ind = TRUE)
          if(length(selectindex)==1){
            data <- rbind(data, as.numeric(data0[selectindex,c('timestep','x','y','center.slp.hPa')]))
          } else if(length(selectindex)>1){
            if(nrow(data)>0){
              candidatesslps <- as.numeric(data0[selectindex,'center.slp.hPa'])
              previousdatarowslp <- tail(data[!is.na(data[,4]),],n=1)[4]
              selectindex <- selectindex[which.min(abs(candidatesslps-previousdatarowslp))]
              data <- rbind(data, as.numeric(data0[selectindex,c('timestep','x','y','center.slp.hPa')]))
            } else{
              data <- rbind(data,c(i, NA, NA, NA))
            }
          } else{
            data <- rbind(data,c(i, NA, NA, NA))
          }
        } else{
          data <- rbind(data,c(i, NA, NA, NA))
        }
      }
      data <- data.frame(data)
      colnames(data) <- c('timestep','x','y','center.slp.hPa')
    } else{data <- data.frame(timestep=1:maxtime,x=rep(NA,maxtime),y=rep(NA,maxtime),center.slp.hPa=rep(NA,maxtime))}
    
    #Is medicane or not -by means of the Hart parameters- for the pch attribute
    ismedicane <- !is.na(data$center.slp.hPa)
    
    if(!is.null(ismedicanedf2)){
      if(length(ismedicane) < nrow(ismedicanedf2)){
        ismedicane <- c(rep(NA,(nrow(ismedicanedf2)-length(ismedicane))), ismedicane)
      }
    }
    ismedicanedf2 <- cbind(ismedicanedf2, ismedicane)
    
    #Slp minimum values for the color attribute
    
    slpmin <- data$center.slp.hPa
    if(!is.null(slpmindf2)){
      if(length(slpmin) < nrow(slpmindf2)){
        slpmin <- c(rep(NA,(nrow(slpmindf2)-length(slpmin))), slpmin)
      }
    }
    slpmindf2 <- cbind(slpmindf2, slpmin)
    
    rut2 <- c(rut2, absmaxtime-maxtime)
    
  }
  
  colnames(ismedicanedf2) <- casenames
  colnames(slpmindf2) <- casenames
  
  
  datafr2 <- data.frame(x=1:nrow(ismedicanedf2))
  for(cc in 1:ncol(ismedicanedf2)){
    vec <- ismedicanedf2[,cc]*(cc)
    vec[vec==0] <- NA
    datafr2 <- cbind(datafr2,vec)
  }
  colnames(datafr2) <- c('x',paste('v',1:ncol(ismedicanedf2),sep=""))
  
  #Crop time to the requested interval
  id <- which(as.character(alldates)==initialdate, arr.ind = TRUE)
  fd <- which(as.character(alldates)==finaldate, arr.ind = TRUE)
  
  datafr <- datafr[id:fd,]
  datafr$x <- 1:length(id:fd)
  datafr2 <- datafr2[id:fd,]
  datafr2$x <- 1:length(id:fd)
  slpmindf <- slpmindf[id:fd,]
  slpmindf2 <- slpmindf2[id:fd,]
  
  # plot
  xlabsseq <- id:fd
  xlabscrop <- alldates[xlabsseq]
  xlabscropall <- sapply(as.character(xlabscrop) , function(x){strsplit(x,':')[[1]][1]}, USE.NAMES = FALSE)
  xlabs <- xlabscropall[seq(1,length(xlabscropall),by=12)]
  xlabsat <- seq(1,length(xlabscropall),by=12)
  
  initdifference <- ((as.numeric(day)-1)*24 + as.numeric(hour)) - ((respecttoday-1)*24 + (respecttohour)) #respect to 2011-11-08 00:00:00
  leglabdays <- as.integer(initdifference/24)
  leglabhours <- (abs(initdifference/24) - as.integer(abs(initdifference/24)))*24
  leglabhours[leglabhours==0] <- "00"
  ylabs <- mapply(function(x,y,z){if(x == 0){paste('+',as.character(y),"d ",as.character(z),"h",sep="")} else{paste("-",as.character(abs(y)),"d ",as.character(z),"h",sep="")}},initdifference,leglabdays,leglabhours,USE.NAMES = FALSE)
  
  
  library(ggplot2)
  
  if(fix=='chem'){
    outname <- paste('~/Escritorio/',paste(Event,chemnochem,'ismedicane.pdf',sep='-'),sep='')
  } 
  if(fix=='nud'){
    outname <-  paste('~/Escritorio/',paste(Event,nonudspnud,'ismedicane.pdf',sep='-'),sep='')
  }
  
  #pdf(outname, 16, 8)
  
  p <- ggplot(data=datafr) + ylim(0,ncol(ismedicanedf)+1) #+ ggtitle(paste(directory1,'vs',directory2,sep=" "))
  
  
  starts <- c()
  ends <- c()
  starts2 <- c()
  ends2 <- c()
  
  #Add simulations to the plot
  for(i in 2:ncol(datafr)){
    #print(i)
    
    colors <- slpmindf[,i-1]
    sizes <- dot.sizes 
    
    colors2 <- slpmindf2[,i-1]
    sizes2 <- dot.sizes
    
    vv <- paste('v',as.character(i-1),sep="")
    
    #Display a dotted grey line in that dates -x points- where the simulation has no started yet
    #Directory1
    rutt <- rut-id
    rutt[rutt<0] <- 0
    notstarted <- rep(NA, nrow(datafr))
    if(rutt[i-1]!=0){
      notstarted[1:(rutt[i-1]+1)] <- i-1
    }
    #Directory2
    rutt2 <- rut2-id
    rutt2[rutt2<0] <- 0
    notstarted2 <- rep(NA, nrow(datafr2))
    if(rutt2[i-1]!=0){
      notstarted2[1:(rutt2[i-1]+1)] <- i-1
    }
    
    #Display also the non-medicane points
    #Directory1
    ismedicanepoints <- datafr[,vv]
    #Directory2
    ismedicanepoints2 <- datafr2[,vv]
    
    
    abc <- !is.na(ismedicanepoints)
    os <- c()
    qs <- c()
    comps <- c()
    for(o in 1:(length(abc)-1)){
      for(q in 1:(length(abc)-1-o)){
        comp <- sum(abc[(o+1):(o+q)])-sum(!abc[(o+1):(o+q)])
        #print(paste(o,q,comp))
        os <- c(os,o)
        qs <- c(qs,q)
        comps <- c(comps,comp)
      }
    }
    
    initial <- os[which(comps==max(comps),arr.ind = TRUE)]
    if(length(initial)>1){initial <- min(initial)}
    final <- initial + qs[which(comps==max(comps),arr.ind = TRUE)]
    if(length(final)>1){final <- max(final)}
    
    cyclonestarttimeidx <- initial
    cycloneendtimeidx <- final
    
    starts <- c(starts, cyclonestarttimeidx)
    ends <- c(ends, cycloneendtimeidx)
    
    
    abc2 <- !is.na(ismedicanepoints2)
    os2 <- c()
    qs2 <- c()
    comps2 <- c()
    for(o2 in 1:(length(abc2)-1)){
      for(q2 in 1:(length(abc2)-1-o2)){
        comp2 <- sum(abc2[(o2+1):(o2+q2)])-sum(!abc2[(o2+1):(o2+q2)])
        #print(paste(o,q,comp))
        os2 <- c(os2,o2)
        qs2 <- c(qs2,q2)
        comps2 <- c(comps2,comp2)
      }
    }
    initial2 <- os2[which(comps2==max(comps2),arr.ind = TRUE)]
    if(length(initial2)>1){initial2 <- min(initial2)}
    final2 <- initial2 + qs2[which(comps2==max(comps2),arr.ind = TRUE)]
    if(length(final2)>1){final2 <- max(final2)}
    
    cyclonestarttimeidx2 <- initial2
    cycloneendtimeidx2 <- final2
    
    #not medicane
    nm <- is.na(datafr[,vv])
    nm[!is.na(notstarted)] <- NA
    nm[nm==FALSE] <- NA
    nm[nm==TRUE] <- 1
    nm <- nm*(i-1)
    
    nm2 <- is.na(datafr2[,vv])
    nm2[!is.na(notstarted2)] <- NA
    nm2[nm2==FALSE] <- NA
    nm2[nm2==TRUE] <- 1
    nm2 <- nm2*(i-1)
    
    startvector <- rep(NA, length(abc))
    startvector[cyclonestarttimeidx] <- i-1
    endvector <- rep(NA, length(abc))
    endvector[cycloneendtimeidx] <- i-1
    insidepoints <- rep(NA,length(abc))
    insidepoints[cyclonestarttimeidx:cycloneendtimeidx] <- i-1
    
    startvector2 <- rep(NA, length(abc2))
    startvector2[cyclonestarttimeidx2] <- i-1
    endvector2 <- rep(NA, length(abc2))
    endvector2[cycloneendtimeidx2] <- i-1
    insidepoints2 <- rep(NA,length(abc2))
    insidepoints2[cyclonestarttimeidx2:cycloneendtimeidx2] <- i-1
    
    starts2 <- c(starts2, cyclonestarttimeidx2)
    ends2 <- c(ends2, cycloneendtimeidx2)
    
    #Directory1
    if(sum(!is.na(notstarted))!=0){
      p <- p + geom_line(aes_string(x=factor(datafr$x), y=notstarted*inbetweenfactor, group=1), size = 1, color = 1, linetype = "dashed", alpha = 0.22)
    }
    #Directory2
    if(sum(!is.na(notstarted2))!=0){
      p <- p + geom_line(aes_string(x=factor(datafr2$x), y=notstarted2*inbetweenfactor+displacement, group=1), size = 1, color = 1, linetype = "dashed", alpha = 0.22)
    }
    
    p <- p + geom_point(aes_string(x=factor(datafr$x), y=datafr[,vv]*inbetweenfactor, colour = colors), size=sizes)#, color=rgb(0.2,0.7,0.1,0.5))
    p <- p + geom_point(aes_string(x=factor(datafr2$x), y=datafr2[,vv]*inbetweenfactor+displacement,  colour = colors2), size=sizes2)#, color=rgb(0.2,0.7,0.1,0.5))
    
    p <- p + geom_point(aes_string(x=factor(datafr$x), y=nm*inbetweenfactor, colour = colors), shape=4, size=sizes, alpha=0.22)#, color=rgb(0.2,0.7,0.1,0.5))
    p <- p + geom_point(aes_string(x=factor(datafr2$x), y=nm2*inbetweenfactor+displacement,  colour = colors2), shape=4, size=sizes2, alpha=0.22)#, color=rgb(0.2,0.7,0.1,0.5))
    
    p <- p + geom_rect(aes_string(xmin=cyclonestarttimeidx+0.5, xmax=cycloneendtimeidx+0.5, ymin=(inbetweenfactor*(i-1))-1, ymax=(inbetweenfactor*(i-1))+1),color='grey',size=0.1,fill=NA,alpha=0.04)
    p <- p + geom_rect(aes_string(xmin=cyclonestarttimeidx2+0.5, xmax=cycloneendtimeidx2+0.5, ymin=((inbetweenfactor*(i-1))+displacement)-1, ymax=((inbetweenfactor*(i-1))+displacement)+1),color='grey',size=0.1,fill=NA,alpha=0.08)
    
  }
  
  if(!is.null(scale_lower_limit) & !is.null(scale_upper_limit)){
    p <- p + scale_colour_gradientn(name="Min SLP (hPa)",colours = Spectral, limits=c(scale_lower_limit,scale_upper_limit), oob=scales::squish) 
  } else{
    p <- p + scale_colour_gradientn(name="Min SLP (hPa)",colours = Spectral) 
  }
  
  if(displacement<0){
    secondyaxislabs <- rep(paste(dir1identifier,rep('\n',1),dir2identifier,sep=""),ncol(ismedicanedf))
    secondyaxislabs1 <- rep(dir1identifier,ncol(ismedicanedf))
    secondyaxislabs2 <- rep(dir2identifier,ncol(ismedicanedf))
    secondyaxisbreaks <-  (seq(1,ncol(ismedicanedf),by=1)*inbetweenfactor)+(displacement/2)
  } else{
    secondyaxislabs <- rep(paste(dir2identifier,'\n',dir1identifier,sep=""),ncol(ismedicanedf))
    secondyaxisbreaks <- c(rbind(seq(1,ncol(ismedicanedf),by=1)*2, seq(1,ncol(ismedicanedf),by=1)*2+displacement)) 
  }
  if(fix=='chem'){
    if(same==TRUE){
      p <- p + theme_light() + ggtitle(label = paste(Event,' : ',chemnochem,'-',nonudspnud,sep='')) + 
        theme(plot.title=element_text(color = "black", size = 10, face = "bold"), panel.border = element_blank(), legend.text = element_text(family='serif'), legend.title = element_text(family='serif'), axis.title=element_text(size=10,face="bold",family='serif'), axis.text.x = element_text(family='serif'), axis.text.y = element_text(family='serif')) + 
        scale_x_discrete(name="Observation time", breaks=xlabsat, labels=xlabs) +
        scale_y_continuous(name=paste("Run-up time (from", respectdate,")"), breaks=(seq(1,ncol(ismedicanedf),by=1)*inbetweenfactor)+(displacement/2), labels=ylabs) +
        scale_shape_manual(name="No medicane", labels=c('|B| > 20',expression(paste('V'[' T']^' L', ' <  0')),expression(paste('V'[' T']^' U', ' <  0')),expression(paste('V'[' T']^' L', ' < V'[' T']^' U'))), values = c(2,3,4,5)) + 
        scale_size(name = "category K")
    } else{
      p <- p + theme_light() + ggtitle(label = paste(Event,chemnochem,sep=' : ')) + 
        theme(plot.title=element_text(color = "black", size = 10, face = "bold"), panel.border = element_blank(), legend.text = element_text(family='serif'), legend.title = element_text(family='serif'), axis.title=element_text(size=10,face="bold",family='serif'), axis.text.x = element_text(family='serif'), axis.text.y = element_text(family='serif')) + 
        scale_x_discrete(name="Observation time", breaks=xlabsat, labels=xlabs) +
        scale_y_continuous(name=paste("Run-up time (from", respectdate,")"), breaks=(seq(1,ncol(ismedicanedf),by=1)*inbetweenfactor)+(displacement/2), labels=ylabs,
                           sec.axis = dup_axis(~., name = "", breaks=secondyaxisbreaks, labels = secondyaxislabs)) +
        scale_shape_manual(name="No medicane", labels=c('|B| > 20',expression(paste('V'[' T']^' L', ' <  0')),expression(paste('V'[' T']^' U', ' <  0')),expression(paste('V'[' T']^' L', ' < V'[' T']^' U'))), values = c(2,3,4,5)) + 
        scale_size(name = "category K")
    }
  } 
  if(fix=='nud'){
    if(same==TRUE){
      p <- p + theme_light() + ggtitle(label = paste(Event,' : ',chemnochem,'-',nonudspnud,sep='')) + 
        theme(plot.title=element_text(color = "black", size = 10, face = "bold"), panel.border = element_blank(), legend.text = element_text(family='serif'), legend.title = element_text(family='serif'), axis.title=element_text(size=10,face="bold",family='serif'), axis.text.x = element_text(family='serif'), axis.text.y = element_text(family='serif')) + 
        scale_x_discrete(name="Observation time", breaks=xlabsat ,labels=xlabs) + #breaks=seq(1,length(datacomp$Date),by=12)
        scale_y_continuous(name=paste("Run-up time (from", respectdate,")"), breaks=(seq(1,ncol(ismedicanedf),by=1)*inbetweenfactor)+(displacement/2), labels=ylabs) +
        scale_shape_manual(name="No medicane", labels=c('|B| > 20',expression(paste('V'[' T']^' L', ' <  0')),expression(paste('V'[' T']^' U', ' <  0')),expression(paste('V'[' T']^' L', ' < V'[' T']^' U'))), values = c(2,3,4,5)) + 
        scale_size(name = "category K")
    } else{
      p <- p + theme_light() + ggtitle(label = paste(Event,nonudspnud,sep=' : ')) + 
        theme(plot.title=element_text(color = "black", size = 10, face = "bold"), panel.border = element_blank(), legend.text = element_text(family='serif'), legend.title = element_text(family='serif'), axis.title=element_text(size=10,face="bold",family='serif'), axis.text.x = element_text(family='serif'), axis.text.y = element_text(family='serif')) + 
        scale_x_discrete(name="Observation time", breaks=xlabsat ,labels=xlabs) + #breaks=seq(1,length(datacomp$Date),by=12)
        scale_y_continuous(name=paste("Run-up time (from", respectdate,")"), breaks=(seq(1,ncol(ismedicanedf),by=1)*inbetweenfactor)+(displacement/2), labels=ylabs,
                           sec.axis = dup_axis(~., name = "", breaks=secondyaxisbreaks, labels = secondyaxislabs)) +
        scale_shape_manual(name="No medicane", labels=c('|B| > 20',expression(paste('V'[' T']^' L', ' <  0')),expression(paste('V'[' T']^' U', ' <  0')),expression(paste('V'[' T']^' L', ' < V'[' T']^' U'))), values = c(2,3,4,5)) + 
        scale_size(name = "category K")
    }
  }
  
  
  
  print(p)
  #dev.off()
}

#}
