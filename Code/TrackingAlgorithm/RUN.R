
args <- commandArgs(trailingOnly = TRUE)

inputdir <- args[1]
if(length(args)==2){
  ncores <- args[2]
} else{
  ncores <- NULL 
}

source('./TrackingAlgorithm/ReadNamelist.R')
parlist <- ReadNamelist()

setwd(paste(getwd(),'TrackingAlgorithm',sep='/'))
source('./ParallelMasterTracker.R')

if(is.null(ncores)){
   ParallelMasterTracker(inputdir, ParamsList = parlist)
} else{
   ParallelMasterTracker(inputdir, ncores = as.numeric(ncores), ParamsList = parlist)
}
