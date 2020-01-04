conditionalSLP <- function(currentimestampslpfield,twoprevioustimestampcenter,previoustimestampcenter,currentimestampcentercandidate,previnstspeed,timestep,resolution){
 #Condition 1: deltat <= 24h
  if(timestep > 24){
    cond1 <- FALSE
  } else{
    cond1 <- TRUE
  }
  
  #Condition 3: V=deltad/deltat < 40 m·s-1
  deltad <- (sqrt((previoustimestampcenter[1] - currentimestampcentercandidate[1])^2 + (previoustimestampcenter[2] - currentimestampcentercandidate[2])^2))*resolution
  
  instspeedkmh <- (deltad/timestep) 
  instspeedms <- instspeedkmh/3.6
  if(instspeedms*3.6 > 40){ #km·h-1 to m·s-1
    cond3 <- FALSE
  } else{
    cond3 <- TRUE
  }
  
  #Condition 4: deltad < deltadmax | deltadmax=max(500km,3*deltat*Vprev)
  deltadmaxx <- max(c(500/3.6,3*timestep*previnstspeed))
  if(deltad > deltadmaxx){
    cond4 <- FALSE
  } else{
    cond4 <- TRUE
  }
  
  #Condition 5: angle restriction for the motion direction change
  angles <- c(Inf,135,90,75,60,45)
  if(instspeedms*3.6 < 10.0){
    anglelimit <- angles[1]
  }
  if(instspeedms*3.6 >= 10.0 && instspeedms*3.6 < 12.5){
    anglelimit <- angles[2]
  }
  if(instspeedms*3.6 >= 12.5 && instspeedms*3.6 < 15.0){
    anglelimit <- angles[3]
  }
  if(instspeedms*3.6 >= 15.0 && instspeedms*3.6 < 17.5){
    anglelimit <- angles[3]
  }
  if(instspeedms*3.6 >= 17.5 && instspeedms*3.6 < 20.0){
    anglelimit <- angles[3]
  }
  if(instspeedms*3.6 >= 20.0 && instspeedms*3.6 < 22.5){
    anglelimit <- angles[4]
  }
  if(instspeedms*3.6 >= 22.5 && instspeedms*3.6 < 25.0){
    anglelimit <- angles[4]
  }
  if(instspeedms*3.6 >= 25.0 && instspeedms*3.6 < 27.5){
    anglelimit <- angles[5]
  }
  if(instspeedms*3.6 >= 27.5 && instspeedms*3.6 < 30.0){
    anglelimit <- angles[5]
  }
  if(instspeedms*3.6 >= 30.0){
    anglelimit <- angles[6]
  }
  
  vector1x <- as.numeric(previoustimestampcenter[1]-twoprevioustimestampcenter[1])
  vector1y <- as.numeric(previoustimestampcenter[2]-twoprevioustimestampcenter[2])
  vector1mod <- sqrt(vector1x^2 + vector1y^2)
  vector2x <- as.numeric(currentimestampcentercandidate[1]-previoustimestampcenter[1])
  vector2y <- as.numeric(currentimestampcentercandidate[2]-previoustimestampcenter[2])
  vector2mod <- sqrt(vector2x^2 + vector2y^2)
  v1v2dotproduct <- vector1x*vector2x + vector1y*vector2y
  if(vector1mod==0 | vector2mod==0){
    angle <- 0
  } else{
    cosangle <- v1v2dotproduct/(vector1mod*vector2mod)
    zcross <- vector1x*vector2y - vector2x*vector1y
    if(zcross < 0){
      angle <- 360 - acos(round(cosangle,8))*360/(2*pi) #rad to deg
    } 
    if(zcross >= 0){
      angle <- acos(round(cosangle,8))*360/(2*pi) #rad to deg
    }
    
  }
 
  
  if(angle > anglelimit){
    cond5 <- FALSE
  } else{
    cond5 <- TRUE
  }
  
  if(cond1 == TRUE && cond3 == TRUE && cond4 == TRUE){ 
    return(TRUE)
  }
  if(cond1 == FALSE | cond3 == FALSE | cond4 == FALSE){
    return(FALSE)
  }
  
  
}
