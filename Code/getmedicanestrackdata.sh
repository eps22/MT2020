#!/bin/bash

folder=$1
ttype=$2

echo "$ttype"

if [ $ttype == "reduced" ]
  then
    R -e 'source("./PostProcessing/getTrackData.R"); getTrackData("'$folder'", complete=FALSE, save=TRUE)'
fi

if [ $ttype == "complete" ]
  then
    R -e 'source("./PostProcessing/getTrackData.R"); getTrackData("'$folder'", complete=TRUE, save=TRUE)'
fi
