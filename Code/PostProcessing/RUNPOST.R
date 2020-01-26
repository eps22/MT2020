
args <- commandArgs(trailingOnly = FALSE)

folder <- args[2]
adjust <- as.logical(args[3])
complete <- as.logical(args[4])
connect <- as.logical(args[5])

source('./PostProcessing/PlotTrack.R')

if(length(args)==6){ #resolution and dthreshold are not given, and then connection should not be made
  PlotTrack(folder=folder, connect=connect, adjust=adjust, complete = complete)
}

if(length(args)==10){ #resolution and dthreshold are given, and then connection should be made
  SLPexpandN <- args[6]
  resolution <- args[7]
  dthreshold <- args[8]
  dtthreshold <- args[9]
  PlotTrack(folder=folder, adjust=adjust, complete = complete, connect=connect, SLPexpandN=as.numeric(SLPexpandN), resolution=as.numeric(resolution), dthreshold=as.numeric(dthreshold), dtthreshold=as.numeric(dtthreshold))
}


