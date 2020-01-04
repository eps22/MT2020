
args <- commandArgs(trailingOnly = FALSE)

inputdir <- args[2]
resol <- args[3]
ts <- args[4]
if(length(args)==6){
  ncores <- args[5]
} else{
  ncores <- NULL 
}


source('./TrackingAlgorithm/ReadNamelist.R')
parlist <- ReadNamelist()

setwd(paste(getwd(),'TrackingAlgorithm',sep='/'))
source('./ParallelMasterTracker.R')

if(is.null(ncores)){
   ParallelMasterTracker(inputdir, Resolution = as.numeric(resol), timestep = as.numeric(ts), ParamsList = parlist)
} else{
   ParallelMasterTracker(inputdir, Resolution = as.numeric(resol), timestep = as.numeric(ts), ncores = as.numeric(ncores), ParamsList = parlist)
}
