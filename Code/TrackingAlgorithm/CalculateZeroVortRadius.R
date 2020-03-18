CalculateZeroVortRadius <- function(resolution, uv10curl, currentcenter, threshold='mean', npoints=30, maxasymm = 300, minsymmdirs=6){
  
  #RIGHT DIRECTION
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
  
  #LEFT DIRECTION
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
  
  #UP DIRECTION
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
  
  #DOWN DIRECTION
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
  
  #UP RIGHT DIRECTION
  if((dim(uv10curl)[1]-currentcenter[1])>(npoints+1)){
    if((dim(uv10curl)[2]-currentcenter[2])>(npoints+1)){
      valuesupright <- c()
      for(i in 1:npoints){
        valuesupright <- c(valuesupright,uv10curl[(currentcenter[1]+i),(currentcenter[2]+i)])  
      }
    } else{
      limmuprighty <- dim(uv10curl)[2]-currentcenter[2]
      valuesupright <- c()
      for(i in 1:limmuprighty){
        valuesupright <- c(valuesupright,uv10curl[(currentcenter[1]+i),(currentcenter[2]+i)])  
      }
    }
  } else{
    limmuprightx <-  dim(uv10curl)[1]-currentcenter[1]
    if((dim(uv10curl)[2]-currentcenter[2])>(npoints+1)){
      valuesupright <- c()
      for(i in 1:limmuprightx){
        valuesupright <- c(valuesupright,uv10curl[(currentcenter[1]+i),(currentcenter[2]+i)])  
      }
    } else{
      limmuprighty <- dim(uv10curl)[2]-currentcenter[2]
      limmupright <- min(c(limmuprightx,limmuprighty))
      valuesupright <- c()
      for(i in 1:limmupright){
        valuesupright <- c(valuesupright,uv10curl[(currentcenter[1]+i),(currentcenter[2]+i)])  
      }
    }
  }
  if(threshold=='mean'){ 
    valuesupright[valuesupright>mean(uv10curl)] <- 1
    valuesupright[valuesupright<mean(uv10curl)] <- -1
  }
  signchange5 <- which(diff(sign(valuesupright))!=0)
  if(length(signchange5)==0){dupright <- 1e10} else{dupright <- (min(signchange5)+0.5)*resolution}
  
  #UP LEFT DIRECTION
  if(currentcenter[1]>(npoints+1)){
    if((dim(uv10curl)[2]-currentcenter[2])>(npoints+1)){
      valuesupleft <- c()
      for(i in 1:npoints){
        valuesupleft <- c(valuesupleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]+i)])  
      }
    } else{
      limmuplefty <- dim(uv10curl)[2]-currentcenter[2]
      valuesupleft <- c()
      for(i in 1:limmuplefty){
        valuesupleft <- c(valuesupleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]+i)])  
      }
    }
  } else{
    limmupleftx <-  currentcenter[1]-1
    if((dim(uv10curl)[2]-currentcenter[2])>(npoints+1)){
      valuesupleft <- c()
      for(i in 1:limmupleftx){
        valuesupleft <- c(valuesupleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]+i)])  
      }
    } else{
      limmuplefty <- dim(uv10curl)[2]-currentcenter[2]
      limmupleft <- min(c(limmupleftx,limmuplefty))
      valuesupleft <- c()
      for(i in 1:limmupleft){
        valuesupleft <- c(valuesupleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]+i)])  
      }
    }
  }
  signchange6 <- which(diff(sign(valuesupleft))!=0)
  if(length(signchange6)==0){dupleft <- 1e10} else{dupleft <- (min(signchange6)+0.5)*resolution}
  
  #DOWN RIGHT DIRECTION
  if((dim(uv10curl)[1]-currentcenter[1])>(npoints+1)){
    if(currentcenter[2]>(npoints+1)){
      valuesdownright <- c()
      for(i in 1:npoints){
        valuesdownright <- c(valuesdownright,uv10curl[(currentcenter[1]+i),(currentcenter[2]-i)])  
      }
    } else{
      limmdownrighty <- currentcenter[2]-1
      valuesdownright <- c()
      for(i in 1:limmdownrighty){
        valuesdownright <- c(valuesdownright,uv10curl[(currentcenter[1]+i),(currentcenter[2]-i)])  
      }
    }
  } else{
    limmdownrightx <-  dim(uv10curl)[1]-currentcenter[1]
    if(currentcenter[2]>(npoints+1)){
      valuesdownright <- c()
      for(i in 1:limmdownrightx){
        valuesdownright <- c(valuesdownright,uv10curl[(currentcenter[1]+i),(currentcenter[2]-i)])  
      }
    } else{
      limmdownrighty <- currentcenter[2]-1
      limmdownright <- min(c(limmdownrightx,limmdownrighty))
      valuesdownright <- c()
      for(i in 1:limmdownright){
        valuesdownright <- c(valuesdownright,uv10curl[(currentcenter[1]+i),(currentcenter[2]-i)])  
      }
    }
  }
  if(threshold=='mean'){ 
    valuesdownright[valuesdownright>mean(uv10curl)] <- 1
    valuesdownright[valuesdownright<mean(uv10curl)] <- -1
  }
  signchange7 <- which(diff(sign(valuesdownright))!=0)
  if(length(signchange7)==0){ddownright <- 1e10} else{ddownright <- (min(signchange7)+0.5)*resolution}
  
  #DOWN LEFT DIRECTION
  if(currentcenter[1]>(npoints+1)){
    if(currentcenter[2]>(npoints+1)){
      valuesdownleft <- c()
      for(i in 1:npoints){
        valuesdownleft <- c(valuesdownleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]-i)])  
      }
    } else{
      limmdownlefty <- currentcenter[2]-1
      valuesdownleft <- c()
      for(i in 1:limmdownlefty){
        valuesdownleft <- c(valuesdownleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]-i)])  
      }
    }
  } else{
    limmdownleftx <-  currentcenter[1]-1
    if(currentcenter[2]>(npoints+1)){
      valuesdownleft <- c()
      for(i in 1:limmdownleftx){
        valuesdownleft <- c(valuesdownleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]-i)])  
      }
    } else{
      limmdownlefty <- currentcenter[2]-1
      limmdownleft <- min(c(limmdownleftx,limmdownlefty))
      valuesdownleft <- c()
      for(i in 1:limmdownleft){
        valuesdownleft <- c(valuesdownleft,uv10curl[(currentcenter[1]-i),(currentcenter[2]-i)])  
      }
    }
  }
  if(threshold=='mean'){ 
    valuesdownleft[valuesdownleft>mean(uv10curl)] <- 1
    valuesdownleft[valuesdownleft<mean(uv10curl)] <- -1
  }
  signchange8 <- which(diff(sign(valuesdownleft))!=0)
  if(length(signchange8)==0){ddownleft <- 1e10} else{ddownleft <- (min(signchange8)+0.5)*resolution}
  
  
  
  vectorofds <- c(dright, dleft, dup, ddown,dupright,dupleft,ddownright,ddownleft)
  #print(sum(vectorofds!=1e10))
  if(sum(vectorofds!=1e10)>=minsymmdirs){
    zerovortradius <- round(mean(vectorofds[vectorofds!=1e10]),1)
    maxdeviation <- max(abs(vectorofds[vectorofds!=1e10]-zerovortradius))
    #print(maxdeviation)
    if(maxdeviation>maxasymm){zerovortradius <- 1e10} #If, in any of the four directions, the calculated radius is higher or lower than the mean in more than 300 km, then
    #there exists an asymmetry that discards the point as medicane center candidate. 
  } else{zerovortradius <- 1e10}
  
  return(zerovortradius)
}
