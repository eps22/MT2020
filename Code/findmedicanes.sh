#!/bin/bash

localdir=$(pwd)

PASSED=$1

if [[ -d $PASSED ]]; then
    if [ -d "$PASSED/output" ]; then
       mkdir -p outfolder/output
       ln -sf "$PASSED"/output/outputfile-* "outfolder/output"
    elif [ -f "$PASSED/outputfile-slp.nc" ]; then
       mkdir -p outfolder/output
       ln -sf "$PASSED"/outputfile-* outfolder/output
    else 
       echo "$PASSED is not a valid input for the algorithm. Please check the README file."
       exit 1
    fi
elif [[ -f $PASSED ]]; then
    ln -sf "$PASSED" "inputfile"
    python3.7 pinterpy/interpy.py pinterpy/interp-namelist
    mkdir -p outfolder/output
    mv outputfile-*.nc outfolder/output
    find . -type l | xargs rm
else
    echo "$PASSED is not a valid input for the algorithm. Please check the README file."
    exit 1
fi

Rscript TrackingAlgorithm/RUN.R $localdir/outfolder ${@:2}

mv outfolder-track.RData track.RData
echo "PWD: "$localdir/$'   \n'"INPUTDIR: ""$PASSED" >> outfolder/inputinfo
mv track.RData outfolder
abc="completed$(ls -d completed* | wc -l)"
mv outfolder $abc
mv $abc/track.RData $abc/$abc-track.RData


