#!/bin/bash

localdir=$(pwd)

ln -sf $1 inputfile 

python3 pinterpy/interpy.py pinterpy/interp-namelist

mkdir outfolder
mkdir outfolder/output
mv outputfile-*.nc outfolder/output
find . -type l | xargs rm

R < TrackingAlgorithm/RUN.R $localdir/outfolder $2 $3 $4 --no-save

mv outfolder-track.RData $1-track.RData
mv outfolder $1
mv $1-track.RData $1
