#!/usr/bin/perl -w
############################### This program is free software; you can redistribute it and/or modify
############################### it under the terms of the GNU General Public License as published by
############################### the Free Software Foundation; either version 2 of the License, or
############################### (at your option) any later version.
###############################
############################### This program is distributed in the hope that it will be useful,
############################### but WITHOUT ANY WARRANTY; without even the implied warranty of
############################### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
############################### GNU General Public License for more details.
###############################
############################### You should have received a copy of the GNU General Public License
############################### along with this program; if not, write to the Free Software
############################### Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
###############################
############################### Copyright 2019 Damon P. Little
###############################
############################### INCLUDE
use strict;
my $i = '';
for(my $k = $#ARGV; $k >= 0; $k--){
	if($ARGV[$k] eq '-i'){
		if(-e $ARGV[$k+1]){
			$i = $ARGV[$k+1];
		}
	}
}
if(length($i)){
	open(INFILE, '<', $i) or die("could not open $i!\n");
	my $alignment = ();
	my $sequence = '';
	while(my $line = <INFILE>){
		chomp($line);
		if($line =~ m/^>/){
			if(length($sequence)){
				process($sequence);
			}
			$sequence = '';
		} else {
			$sequence .= $line;
		}
	}
	process($sequence);
	close(INFILE);
	sub process {
		my $sequence = uc($_[0]);
		$sequence =~ tr/ACDEFGHIKLMNPQRSTVWYX\-//cd;
		my @bases = split(//, $sequence);
		my $k = $#{$alignment}+1;
		for(my $j = $#bases; $j >= 0; $j--){
			$alignment->[$k][$j] = $bases[$j];
		}
		return(0);
	}
	my $count = 0;
	my $mean = 0;
	my $M2 = 0;
	for(my $k = $#{$alignment}; $k > 0; $k--){
		for(my $kk = $k-1; $kk >= 0; $kk--){
			my $distance = 0;
			for(my $j = $#{$alignment->[$k]}; $j >= 0; $j--){
				if($alignment->[$k][$j] ne $alignment->[$kk][$j]){
					$distance++;
				}
			}
			$count++;
			my $nextMean = $mean+($distance-$mean)/$count;
			$M2 = $M2+($distance-$mean)*($distance-$nextMean);
			$mean = $nextMean;
		}
	}
	if($count > 1){
		print(int($mean+0.5) . ' ' . int(($M2/($count-1))+0.5) . "\n");
	} else {
		print(int($mean+0.5) . " 0\n");
	}
} else {
	print(STDERR "\neditAA.pl a script to calculate the mean and sample variance of the edit\n");
	print(STDERR "distance from an amino acid FASTA file.\n\n");
	print(STDERR "USAGE: editAA.pl -i file\n");
	print(STDERR "-i\tspecifies an aligned FASTA file (required)\n\n");
}
exit(0);
