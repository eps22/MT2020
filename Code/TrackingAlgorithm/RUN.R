
args <- commandArgs(trailingOnly = FALSE)

inputdir <- args[2]
if(length(args)==4){
  ncores <- args[3]
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
