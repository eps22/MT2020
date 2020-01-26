listodf <- function(trackinglist){
  
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
  
  return(trackingdf)
  
}