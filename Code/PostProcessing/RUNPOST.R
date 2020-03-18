
args <- commandArgs(trailingOnly = FALSE)

folder <- args[2]
adjust <- as.logical(args[3])
complete <- as.logical(args[4])
connect <- as.logical(args[5])

source('./PostProcessing/PlotTrack.R')

if(length(args)==6){ 
  PlotTrack(folder=folder, connect=connect, adjust=adjust, complete = complete)
}

if(length(args)==8){ 
  dthreshold <- args[6]
  dtthreshold <- args[7]
  PlotTrack(folder=folder, adjust=adjust, complete = complete, connect=connect, dthreshold=as.numeric(dthreshold), dtthreshold=as.numeric(dtthreshold))
}

if(length(args)==9){ 
  SLPexpandN <- args[6]
  dthreshold <- args[7]
  dtthreshold <- args[8]
  PlotTrack(folder=folder, adjust=adjust, complete = complete, connect=connect, SLPexpandN=as.numeric(SLPexpandN), dthreshold=as.numeric(dthreshold), dtthreshold=as.numeric(dtthreshold))
}


