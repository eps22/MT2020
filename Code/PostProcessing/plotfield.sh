#!/bin/bash

args="$1"

echo $args

R < plotfield.R "$@" --no-save

