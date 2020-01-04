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
