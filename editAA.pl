#!/usr/bin/perl
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
