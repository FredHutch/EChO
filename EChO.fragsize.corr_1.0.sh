#!/usr/bin/env bash

set -ue

if [ $# -lt 7 ]
then
	echo "

	Compare fragment size profiles of overlapping regions from two CUT&RUN experiments

	Usage: bash EChO.fragsize.corr.sh <first region bed file> <first summit bed file> <first fragment bed file> <second region bedfile> <second summit bed file> <second fragment bed file> <output prefix>
	
	This script computes the true correlation, maximum correlation, and cross-correlation shift for every pair of overlapping regions
	
	Output: <input prefix>.crosscor.txt
	"
	exit
fi

password=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13; echo ''`

module load bedops

bedops -i $1 $4 | bedtools intersect -wao -a - -b $1 | cut -f 1,2,3,4,5,6 | bedtools intersect -wao -a - -b $2 | awk '$7 != "." {print $0}' | bedtools groupby -g 1,2,3 -c 4,5,6,12 -o first,first,first,max | bedtools intersect -wao -a - -b $4 | cut -f 1,2,3,4,5,6,7,8,9,10 | bedtools intersect -wao -a - -b $5 | awk '$11 != "." {print $0}' | bedtools groupby -g 1,2,3 -c 4,5,6,7,8,9,10,16 -o first,first,first,first,first,first,first,min | awk 'BEGIN{s=1}; {print $0"\t"s; s++}' | awk '$3-$2 >= 401 {print $0}' > $password.comb.bed

cut -f 4,5,6,7,12 $password.comb.bed | bedtools intersect -wao -a - -b $3 > $password.first.frags
cut -f 8,9,10,11,12 $password.comb.bed | bedtools intersect -wao -a - -b $6 > $password.second.frags

path=`dirname $0`
len=`wc -l $password.comb.bed | awk '{print $1}'`
for i in $(seq 1 $len)
do
#	center=`sed -n '$i p' $password.comb.bed | cut -f 1,2,3 | awk '{if(($3-$2)%2==0){print (($3-$2)/2)+$2}else{print (($3-$2+1)/2)+$2}}'` 
	num=$(($i + 0))
	awk -v which="$num" '$12==which {print $0}' $password.comb.bed > $password.int.bed
	awk -v which="$num" '$5==which {if(($8-$7)%2==0){fragcen=(($8-$7)/2)+$7}else{fragcen=(($8-$7+1)/2)+$7}; print fragcen"\t"$9"\t"$4}' $password.first.frags > $password.first.bed
	awk -v which="$num" '$5==which {if(($8-$7)%2==0){fragcen=(($8-$7)/2)+$7}else{fragcen=(($8-$7+1)/2)+$7}; print fragcen"\t"$9"\t"$4}' $password.second.frags > $password.second.bed
#	sed -n "$num p" $password.comb.bed | cut -f 1,2,3 > $password.int.bed
#	sed -n "$num p" $password.comb.bed | cut -f 4,5,6,7 | bedtools intersect -wao -a - -b $3 | awk '{if(($7-$6)%2==0){fragcen=(($7-$6)/2)+$6}else{fragcen=(($7-$6+1)/2)+$6}; print fragcen"\t"$8"\t"$4}' > $password.first.bed
#	sed -n "$num p" $password.comb.bed | cut -f 8,9,10,11 | bedtools intersect -wao -a - -b $6 | awk '{if(($7-$6)%2==0){fragcen=(($7-$6)/2)+$6}else{fragcen=(($7-$6+1)/2)+$6}; print fragcen"\t"$8"\t"$4}' > $password.second.bed
	lenfirst=`wc -l $password.first.bed | awk '{print $1}'`
	lensecond=`wc -l $password.second.bed | awk '{print $1}'`
	if [ $lenfirst -gt 10 ] && [ $lensecond -gt 10 ]
	then
		Rscript $path/EChO.fragsize.corr_1.0.R --first=$password.first.bed --second=$password.second.bed --int=$password.int.bed --output=$7
	fi
done

#I need to process the -wao output with a perl script before this step to separate fragments  by region


rm $password.comb.bed
rm $password.int.bed
rm $password.first.bed
rm $password.second.bed
rm $password.first.frags
rm $password.second.frags
