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
	my $score = 0;
	my $realColumns = 0;
	my @aa = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');
	for(my $j = $#{$alignment->[0]}; $j >= 0; $j--){
		my %column = ();
		for(my $k = $#{$alignment}; $k >= 0; $k--){
			if(!exists($column{$alignment->[$k][$j]})){
				$column{$alignment->[$k][$j]} = 1;
			}
		}
		my $count = 0;
		for(my $k = $#aa; $k >= 0; $k--){
			if(exists($column{$aa[$k]})){
				$count++;
			}
		}
		if($count != 0){
			$realColumns++;
			if(exists($column{'-'})){
				$count++;
			}
		}
		$score += $count;
	}
	if($realColumns > 1){
		print(int(($score/$realColumns)+0.5) . "\n");
	} else {
		print("0\n");
	}
} else {
	print(STDERR "\nalignAA.pl a script to calculate the Nixon/Little alignment score\n");
	print(STDERR "from an amino acid FASTA file.\n\n");
	print(STDERR "USAGE: alignAA.pl -i file\n");
	print(STDERR "-i\tspecifies an aligned FASTA file (required)\n\n");
}
exit(0);
