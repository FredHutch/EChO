#!/usr/bin/env bash

set -ue

if [ $# -lt 4 ] 
then
	echo "
	
	Enhanced Chromatin Occupancy (EChO): Find local minima in average fragment size profiles from CUT&RUN data
	
	Requires bedtools, perl, and R to be in your path


	Usage: EChO_1.0.sh <region bed file> <fragment bed file> ["foci" | "matrix"] <output prefix>


	Foci mode generates average fragment size profiles from a set of regions and identifies local minima in those profiles (foci)
	
	Matrix mode similarly generates average fragment size profiles and generates a matrix of average fragment size values surrounding every focus
	NOTE: Matrix mode requires the fourth column of the region bed file to represent the region focus as defined by EChO Foci mode

	"
	exit 1
fi
#module load bedtools

path=`dirname $0`

mode=`echo $3`

columns=`awk '{print NF}' $1 | head -n 1`

if [[ $mode == "foci" ]]
then
	bedtools intersect -wao -a $1 -b $2 | perl $path/EChO_1.0.pl - $4 line $path/EChO_1.0.R $columns | sed '/^\s*$/d' > $4.EChO.bed
elif [[ $mode == "matrix" ]]
then	
	bedtools intersect -wao -a $1 -b $2 | perl $path/EChO.matrix_1.0.pl - $4 line $path/EChO.matrix_1.0.R $columns | sed '/^\s*$/d' > $4.EChO.matrix
else
	echo "**Third entry must be \"foci\" for foci mode or \"matrix\" for matrix mode**


	Enhanced Chromatin Occupancy (EChO): Find local minima in average fragment size profiles from CUT&RUN data
	
	Requires bedtools, perl, and R to be in your path


	Usage: EChO_1.0.sh <region bed file> <fragment bed file> ["foci" | "matrix"] <output prefix>


	Foci mode generates average fragment size profiles from a set of regions and identifies local minima in those profiles (foci)
	
	Matrix mode similarly generates average fragment size profiles and generates a matrix of average fragment size values surrounding every focus
	NOTE: Matrix mode requires the fourth column of the region bed file to represent the region focus as defined by EChO Foci mode

	"
	exit 1
fi
