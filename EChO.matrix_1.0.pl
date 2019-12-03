#!/usr/bin/perl

#This script is designed to take an input generated as follows:
#
# 1) Generate an enriched region .bed file of the following structure:
#
# <chr> <start> <stop> <single bp coordinate of summit>
#
# 2) use bedtools to intersect the above .bed file with a paired-end fragment .bed file:
#
# bedtools intersect -wao -a <regions .bed file> -b <fragments .bed file>
#
# 3) Use the resulting file as an input to this script:
#
# perl scripts/fragsize.matrix.pl <intersected .bed file> <output name>

use strict;
#use Statistics::R;

my @newoutput;
my %points;
my $chr=0;
my $start=0;
my $stop=0;
my $center=0;
#my $center2=0;
#my $center3=0;
my $i=1;
my $j=0;
my $bedlen = $ARGV[4]-1;
#if(@ARGV < 3){
#	die "Three argruments required: intersected bed file, output name, and span value\n";
#}
open(BUFFER, $ARGV[0]);
while(my $line = <BUFFER>){
	chomp($line);
	my @fields = split("\t", $line);
	if($chr ne $fields[0] || $start != $fields[1] || $stop != $fields[2] || $center != $fields[3]){
		unless(!keys %points){
			open(OUT, ">$ARGV[1].test.txt");
			foreach my $key (keys %points){
				print OUT "$points{$key}{offset}\t$points{$key}{len}\n";
			}
			close(OUT);
			$j=0;
			my $output = `Rscript $ARGV[3] --frame=$ARGV[1].test.txt --center=$center --line=$i --output=$ARGV[1] --span=$ARGV[2]`;
			push(@newoutput, $output);
			$i++;
		}
		%points = ();
		$chr = $fields[0];
		$start = $fields[1];
		$stop = $fields[2];
		$center = $fields[3];
	}
	my $fragstart = $fields[$bedlen+2];
	my $fragstop = $fields[$bedlen+3];
	my $fragcenter;
	if(($fragstop - $fragstart)%2 == 1){
		$fragcenter = (($fragstop - $fragstart + 1)/2)+$fragstart;
	}else{
		$fragcenter = (($fragstop-$fragstart)/2)+$fragstart;
	}
	my $offset = $fragcenter - $center;
	my $fraglength = $fields[$bedlen+3]-$fields[$bedlen+2];
	my $index = $j;
	$points{$index}{offset} = $offset;
	$points{$index}{len} = $fraglength;
	$j++;
}
open(OUT, ">$ARGV[1].test.txt");
foreach my $key (keys %points){
	print OUT "$points{$key}{offset}\t$points{$key}{len}\n";
}
close(OUT);
my $output = `Rscript $ARGV[3] --frame=$ARGV[1].test.txt --center=$center --line=$i --output=$ARGV[1] --span=$ARGV[2]`;
push(@newoutput, $output);
foreach(@newoutput){
	print "$_\n";
}
#close(OUT2);
