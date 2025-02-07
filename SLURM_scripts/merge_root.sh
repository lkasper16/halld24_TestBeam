#!/bin/bash

# This script merges the root files for a given run number
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <run_number> <tag>"
    exit 1
fi

run_number=$1
tag=$2

# Set the path to the root files
path="/work/halld2/home/nseptian/halld24_TestBeam/"

# Set the output file name
output="/work/halld2/home/nseptian/halld24_TestBeam/RootOutput"

# Merge the root files
hadd -f $output/out_hits_Run_00${run_number}_$tag.root $path/out_hits_Run_00${run_number}*.root

