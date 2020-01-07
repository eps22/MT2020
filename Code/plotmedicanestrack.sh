#!/bin/bash


if [ $# -eq 2 ]
  then 
    R < PostProcessing/RUNPOST.R $1 $2 --no-save
fi

if [ $# -eq 5 ]
  then
    R < PostProcessing/RUNPOST.R $1 $2 $3 $4 $5 --no-save
fi
