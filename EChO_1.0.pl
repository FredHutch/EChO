#!/usr/bin/perl

#This script reads a bed file from EChO_1.0.sh with fragments mapped
#onto it and reports offset and length values for each fragment mapped.

use strict;

my @newoutput;
my %points;
my $chr=0;
my $start=0;
my $stop=0;
my $center;
my $i = 1;
my $bedlen = $ARGV[4]-1;
open(BUFFER, $ARGV[0]);
while(my $line = <BUFFER>){
	chomp($line);
	my @fields = split("\t", $line);
	if($chr ne $fields[0] || $start != $fields[1] || $stop != $fields[2]){
		unless(keys %points < 10){
			open(OUT, ">$ARGV[1].test.txt");
			foreach my $key (keys %points){
				print OUT "$points{$key}{offset}\t$points{$key}{len}\n";
			}
			close(OUT);
			my $output = `Rscript $ARGV[3] --frame=$ARGV[1].test.txt --center=$center --chr=$chr --output=$ARGV[1] --span=$ARGV[2]`;
			push(@newoutput, $output);
		}
		%points = ();
		$i = 1;
		$chr = $fields[0];
		$start = $fields[1];
		$stop = $fields[2];
		if(($stop - $start)%2 == 1){
			$center = (($stop - $start + 1)/2)+$start;
		}else{
			$center = (($stop-$start)/2)+$start;
		}
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
	my $fraglength = $fields[$bedlen+3] - $fields[$bedlen+2];
	my $index = $i;
	$i++;
	$points{$index}{offset} = $offset;
	$points{$index}{len} = $fraglength;
}
foreach(@newoutput){
	print "$_\n";
}
