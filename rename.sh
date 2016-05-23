#!/bin/bash

# This script searches for all instances of "trajOpt" and 
# replaces it with "optimTraj"

# lower case
files=$(grep -r -l --exclude-dir=".git" --exclude="*.sh" 'trajOpt' .)
for file in $files
do
	sed -i.bak 's/trajOpt/optimTraj/g' $file
done
find . -name "*.bak" -type f -delete

# Upper Case
files=$(grep -r -l --exclude-dir=".git" --exclude="*.sh" 'TrajOpt' .)
for file in $files
do
	sed -i.bak 's/TrajOpt/OptimTraj/g' $file
done
find . -name "*.bak" -type f -delete

