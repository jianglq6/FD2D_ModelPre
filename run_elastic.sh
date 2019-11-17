#!/usr/bin/bash

echo "running: `date`"
currentdir=`pwd`

mkdir -p OUTPUT

#clean output files
rm -rf OUTPUT/*

#cd currentdir

#stores setup
cp DATA/Par_file OUTPUT/

./bin/FD2D_modeling 

#store output
cp DATA/*STATIONS* OUTPUT

echo 
echo "see results in directory: OUTPUT/"
echo
echo "done"
echo `date`
