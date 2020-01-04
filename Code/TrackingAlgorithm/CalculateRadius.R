CalculateRadius <- function(resolution, uv10curl, currentcenter, threshold='mean', npoints=30){

    if((dim(uv10curl)[1]-currentcenter[1])<(npoints+1)){
      valuesright <- uv10curl[currentcenter[1]:dim(uv10curl)[1],currentcenter[2]]
    } else{
      valuesright <- uv10curl[currentcenter[1]:(currentcenter[1]+npoints),currentcenter[2]]
    }
    if(threshold=='mean'){ 
    valuesright[valuesright>mean(uv10curl)] <- 1
    valuesright[valuesright<mean(uv10curl)] <- -1
    } #if threshold='zero', then values are not changed
    signchange1 <- which(diff(sign(valuesright))!=0)
    if(length(signchange1)==0){dright <- 1e10} else{dright <- (min(signchange1)+0.5)*resolution}
    

    if(currentcenter[1]<(npoints+1)){
      valuesleft <- uv10curl[currentcenter[1]:1,currentcenter[2]]
    } else{
      valuesleft <- uv10curl[currentcenter[1]:(currentcenter[1]-npoints),currentcenter[2]]
    }
    if(threshold=='mean'){ 
    valuesleft[valuesleft>mean(uv10curl)] <- 1
    valuesleft[valuesleft<mean(uv10curl)] <- -1
    }
    signchange2 <- which(diff(sign(valuesleft))!=0)
    if(length(signchange2)==0){dleft <- 1e10} else{dleft <- (min(signchange2)+0.5)*resolution}


    if((dim(uv10curl)[2]-currentcenter[2])<(npoints+1)){
      valuesup <- uv10curl[currentcenter[1],currentcenter[2]:(dim(uv10curl)[2])]
    } else{
      valuesup <- uv10curl[currentcenter[1],currentcenter[2]:(currentcenter[2]+npoints)]
    }
    if(threshold=='mean'){ 
    valuesup[valuesup>mean(uv10curl)] <- 1
    valuesup[valuesup<mean(uv10curl)] <- -1
    }
    signchange3 <- which(diff(sign(valuesup))!=0)
    if(length(signchange3)==0){dup <- 1e10} else{dup <- (min(signchange3)+0.5)*resolution}
    

    if(currentcenter[2]<(npoints+1)){
      valuesdown <- uv10curl[currentcenter[1],currentcenter[2]:1]
    } else{
      valuesdown <- uv10curl[currentcenter[1],currentcenter[2]:(currentcenter[2]-npoints)]
    }
    if(threshold=='mean'){ 
    valuesdown[valuesdown>mean(uv10curl)] <- 1
    valuesdown[valuesdown<mean(uv10curl)] <- -1
    }
    signchange4 <- which(diff(sign(valuesdown))!=0)
    if(length(signchange4)==0){ddown <- 1e10} else{ddown <- (min(signchange4)+0.5)*resolution}
    

    vectorofds <- c(dright, dleft, dup, ddown)
    zerovortradius <- round(mean(vectorofds),1)
    
    maxdeviation <- max(abs(vectorofds-zerovortradius))
    if(maxdeviation>300){zerovortradius <- 1e10} #If, in any of the four directions, the calculated radius is higher or lower than the mean in more than 300 km, then
                                                 #there exists an asymmetry that discards the point as medicane center candidate. 

    return(zerovortradius)
}
