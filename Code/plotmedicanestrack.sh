#!/bin/bash


if [ $# -eq 4 ]
  then 
    R < PostProcessing/RUNPOST.R $1 $2 $3 $4 --no-save
fi

if [ $# -eq 6 ]
  then
    R < PostProcessing/RUNPOST.R $1 $2 $3 $4 $5 $6 --no-save
fi

if [ $# -eq 7 ]
  then
    R < PostProcessing/RUNPOST.R $1 $2 $3 $4 $5 $6 $7 --no-save
fi
