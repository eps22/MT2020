
args <- commandArgs(trailingOnly = FALSE)

folder <- args[2]

source('./PostProcessing/PlotTrack.R')

if(length(args)==3){ #resolution and dthreshold are not given, and then connection should not be made
  PlotTrack(folder=folder, connect=FALSE)
}

if(length(args)>3){ #resolution and dthreshold are given, and then connection should be made
resolution <- args[3]
dthreshold <- args[4]
dtthreshold <- args[5]
PlotTrack(folder=folder, connect=TRUE, resolution=as.numeric(resolution), dthreshold=as.numeric(dthreshold), dtthreshold=as.numeric(dtthreshold))
}


