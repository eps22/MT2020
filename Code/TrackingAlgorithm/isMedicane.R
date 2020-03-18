isMedicane <- function(compfolder, center, currentimeidx, timestep, resolution, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, PressureVertLevelDimName, SLPVarName, ZVarName){
  
  
  track <- TrackCenterCandidate(compfolder,center,currentimeidx,resolution,timestep,SLPVarName)
  
  if(sum(track$centers[1,]==track$centers[2,])==2){
    abc <- list()
    abc$logical <- FALSE
    abc$triad <- NA
  } else{
    abc <- FulfillsHart(compfolder, currentimeidx, track$centers, resolution, HartConditionsTocheck, Blowerpressurelevel, Bupperpressurelevel, Bmultiplemeasure, Bdirections, Bthreshold, LTWlowerpressurelevel, LTWupperpressurelevel, UTWlowerpressurelevel, UTWupperpressurelevel, zerovortradius2, PressureVertLevelDimName, ZVarName)
  }
  
  return(abc)
  
}
