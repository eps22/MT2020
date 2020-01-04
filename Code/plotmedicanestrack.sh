#!/bin/bash


if [ $# -eq 1 ]
  then 
    R < PostProcessing/RUNPOST.R $1 --no-save
fi

if [ $# -eq 4 ]
  then
    R < PostProcessing/RUNPOST.R $1 $2 $3 $4 --no-save
fi
