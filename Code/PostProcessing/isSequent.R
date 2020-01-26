isSequent <- function(point1, point2, resolution, dthreshold, dtthreshold){
  
  time1 <- point1[1]
  time2 <- point2[1]
  
  if(time2<=(time1+dtthreshold)){timesequent <- TRUE} else{timesequent <- FALSE}
  
  source('./PostProcessing/dd.R')
  
  point1position <- c(point1[2],point1[3])
  point2position <- c(point2[2],point2[3])
  
  if((dd(point1position, point2position)*resolution)<dthreshold){spatialsequent <- TRUE} else{spatialsequent <- FALSE}
  
  if(timesequent && spatialsequent){
    return(TRUE)
  } else{
    return(FALSE)
  }
  
}
