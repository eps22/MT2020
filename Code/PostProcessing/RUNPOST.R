
args <- commandArgs(trailingOnly = FALSE)

folder <- args[2]
adjust <- as.logical(args[3])

source('./PostProcessing/PlotTrack.R')

if(length(args)==4){ #resolution and dthreshold are not given, and then connection should not be made
  PlotTrack(folder=folder, connect=FALSE, adjust=adjust)
}

if(length(args)>4){ #resolution and dthreshold are given, and then connection should be made
resolution <- args[4]
dthreshold <- args[5]
dtthreshold <- args[6]
PlotTrack(folder=folder, connect=TRUE, adjust=adjust, resolution=as.numeric(resolution), dthreshold=as.numeric(dthreshold), dtthreshold=as.numeric(dtthreshold))
}


