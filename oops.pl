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
############################### Copyright 2019 Damon P. Little and Apurva Narechania
###############################
############################### INCLUDE
use Term::ANSIColor;
use Time::HiRes qw(time);
use File::Copy qw(copy);
use strict;
############################### GLOBALS
my $MLmodel = 'PROTGAMMAAUTO';
my $N = 0;
my $bestScores = {}; ### criteria -> min|max -> score
my $c = {};
### changes to @criteria must be reflected in calls to $bestScores, $c, $currentScores, and $currentTimes
my @criteria = (	'(Z) variance of ML support',
						'(Y) mean ML support',
						'(X) normalized sum of ML support',
						'(W) ML tree score',
						'(V) variance of parsimony jackknife',
						'(U) mean parsimony jackknife',
						'(T) normalized sum of parsimony jackknife',
						'(S) ensemble RI',
						'(R) ensemble CI',
						'(Q) parsimony tree length',
						'(P) parsimony information',
						'(O) minimum possible parsimony steps',
						'(N) maximum possible parsimony steps',
						'(M) TCS (transitive consistency score)',
						'(L) sum of pairs BLOSUM62',
						'(K) sum of pairs BLOSUM50',
						'(J) sum of pairs BLOSUM45',
						'(I) sum of pairs',
						'(H) variance of edit distance',
						'(G) mean edit distance',
						'(F) alignment length',
						'(E) align score',
						'(D) random',
						'(C) relative sequence length (s.d. above/below mean)',
						'(B) variance of sequence length',
						'(A) mean sequence length');
my $d = 1;
my $globalMean = 0;
my $globalSD = 0;
my $i = '';
my $m = 'lg';
my $n = ~0;
my $o = '';
my $r = 'raxml';
my $s = 0;
my $sequenceData = {}; ### species -> terminal -> sequence
my $t = 1;
my $terminalCount = {}; ### species -> number of sequences
my $w = 0;
my @criteriaLetters = ();
my @species = ();
############################### INITALIZE GLOBALS
for(my $k = $#criteria; $k >= 0; $k--){
	$bestScores->{$criteria[$k]}->{'max'} = -1*~0;
	$bestScores->{$criteria[$k]}->{'min'} = ~0;
	$c->{$criteria[$k]} = 1;
	$criteriaLetters[$k] = substr($criteria[$k], 1, 1);
}
############################### PARSE OPTIONS
print(STDERR color('red'));
for(my $k = $#ARGV; $k >= 0; $k--){
	if($ARGV[$k] eq '-c'){
		for(my $j = $#criteria; $j >= 0; $j--){
			if(index($ARGV[$k+1], $criteriaLetters[$j]) == -1){
				$c->{$criteria[$j]} = 0;
			} else {
				$c->{$criteria[$j]} = 1;
			}
		}
	} elsif($ARGV[$k] eq '-d'){
		$d = 0;
	} elsif($ARGV[$k] eq '-i'){
		if(-e $ARGV[$k+1]){
			$i = $ARGV[$k+1];
		} else {
			print(STDERR "\nInput file (-i) could not be read!\n");
		}
	} elsif($ARGV[$k] eq '-m'){
		if(lc($ARGV[$k+1]) eq 'wag'){
			$m = 'wag';
		}
	} elsif($ARGV[$k] eq '-n'){
		$ARGV[$k+1] =~ tr/0123456789//cd;
		if(($ARGV[$k+1] > 1) && ($ARGV[$k+1] < ~0)){
			$N = 1;
			$n = $ARGV[$k+1];
		}
	} elsif($ARGV[$k] eq '-o'){
		if(!-e $ARGV[$k+1]){
			$o = $ARGV[$k+1];
		} else {
			print(STDERR "\nOutput directory (-o) exists!\n");
		}
	} elsif($ARGV[$k] eq '-r'){
		$r = $ARGV[$k+1];
	} elsif($ARGV[$k] eq '-s'){
		$s = 1;
	} elsif($ARGV[$k] eq '-t'){
		if($ARGV[$k+1] =~ m/^[0-9]+$/){
			$t = $ARGV[$k+1];
		}
	} elsif($ARGV[$k] eq '-w'){
		$w = 1;
	}
}
print(STDERR color('reset'));
############################### START
if(length($i) && length($o)){
	############################### PARSE INPUT
	open(INFILE, $i) or die("Could not open $i!");
	my $M2 = 0;
	my $count = 0;
	my $currentSpecies = '';
	my $currentTerminal = '';
	my $mean = 0;
	while(my $line = <INFILE>){
		chomp($line);
		if(length($line)){
			if($line =~ m/^>(.+)#(.+)/){
				if(exists($sequenceData->{$currentSpecies}->{$currentTerminal})){
					($M2, $count, $mean) = mean($M2, $count, $mean, length($sequenceData->{$currentSpecies}->{$currentTerminal}));
				}
				$currentSpecies = $1;
				$currentTerminal = $2;
				$sequenceData->{$currentSpecies}->{$currentTerminal} = '';
				$terminalCount->{$currentSpecies}++;
			} else {
				$line = uc($line);
				$line =~ tr/ABCDEFGHIKLMNOPQRSTUVWXYZ//cd;
				$line =~ tr/BJOUZ/XXLCX/;
				$sequenceData->{$currentSpecies}->{$currentTerminal} .= $line;
			}
		}
	}
	($M2, $count, $mean) = mean($M2, $count, $mean, length($sequenceData->{$currentSpecies}->{$currentTerminal}));
	$globalMean = $mean;
	if($count > 1) {
		$globalSD = sqrt($M2/($count-1));
	} else {
		$globalSD = 0;
	}
	close(INFILE);
	############################### DEACTIVATE INAPPLICABLE CRITERIA
	my $sp = scalar(keys(%{$sequenceData}))-1;
	if(($sp < 4) && (($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1) || ($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1) || ($c->{'(W) ML tree score'} == 1) || ($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1))){
		$c->{'(N) maximum possible parsimony steps'} = 0;
		$c->{'(O) minimum possible parsimony steps'} = 0;
		$c->{'(P) parsimony information'} = 0;
		$c->{'(Q) parsimony tree length'} = 0;
		$c->{'(R) ensemble CI'} = 0;
		$c->{'(S) ensemble RI'} = 0;
		$c->{'(T) normalized sum of parsimony jackknife'} = 0;
		$c->{'(U) mean parsimony jackknife'} = 0;
		$c->{'(V) variance of parsimony jackknife'} = 0;
		$c->{'(W) ML tree score'} = 0;
		$c->{'(X) normalized sum of ML support'} = 0;
		$c->{'(Y) mean ML support'} = 0;
		$c->{'(Z) variance of ML support'} = 0;
		print(STDERR color('red'));
		print(STDERR "\nTree–based criteria are not valid for less than four species. The input file\ncontains only $sp, thus all tree based criteria have been disabled.\n");
		print(STDERR color('reset'));
	}
	############################### OUTPUT DIRECTORY
	qx`mkdir -p $o`;
	chdir($o) or die("Cannot change to $o!\n");
	############################### ML MODEL FOR RAXML
	if(($s == 1) && (($c->{'(W) ML tree score'} == 1) || ($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1))){
		open(OUTFILE, ">$i.fa") or die("Cannot open '$o/$i.fa'!\n");
		foreach my $species (sort(keys(%{$sequenceData}))){
			foreach my $terminal (sort(keys(%{$sequenceData->{$species}}))){
				print(OUTFILE formatFASTA($species . '#' . $terminal, $sequenceData->{$species}->{$terminal}));
			}
		}
		close(OUTFILE);
		qx`mafft --thread $t --auto --quiet $i.fa > $i.aln`;
		my $alignedSequenceData = {};
		my $id = '';
		open(INFILE, "$i.aln") or die("Cannot open '$o/$i.aln'!\n");
		while(my $line = <INFILE>){
			chomp($line);
			if($line =~ m/^>(.+)/){
				$id = $1;
			} else {
				$line =~ s/-/?/g;
				$alignedSequenceData->{$id} .= $line;
			}
		}
		my @terminals = sort({$b cmp $a} keys(%{$alignedSequenceData}));
		open(OUTFILE, ">$i.phy") or die("Cannot open '$o/$i.phy'!\n");
		print(OUTFILE scalar(@terminals) . ' ' . length($alignedSequenceData->{$id}) . "\n");
		for(my $k = $#terminals; $k >= 0; $k--){
			print(OUTFILE "$terminals[$k]\t$alignedSequenceData->{$terminals[$k]}\n");
		}
		print(OUTFILE "\n");
		close(OUTFILE);
		qx`$r -T $t -m $MLmodel -s $i.phy -n $i -p 5 -# 1 -x 5`;
		open(INFILE, "RAxML_info.$i") or die("Cannot open '$o/$i.aln'!\n");
		while(my $line = <INFILE>){
			chomp($line);
			if($line =~ m/best-scoring AA model: ([A-Z]+) likelihood -[0-9]+\.[0-9]+ with (empirical|fixed) base frequencies/){
				$MLmodel = 'PROTGAMMA' . $1;
				if($2 eq 'empirical'){
					$MLmodel .= 'F';
				} else {
					$MLmodel .= 'X';
				}
				last;
			}
		}
		close(INFILE);
		if($d == 1){
			unlink(glob($i . '.*'));
			unlink(glob('RAxML_*.' . $i . '*'));
		}
	}
	############################### GENERATE AND SCORE ALL COMBINATIONS
	my $combination = ''; ### csv of species#terminal
	foreach my $sp (sort({rand() <=> 0.5} keys(%{$terminalCount}))){
		if($terminalCount->{$sp} == 1){
			my @key = keys(%{$sequenceData->{$sp}});
			$combination .= $sp . '#' . $key[0] . ',';
		} else {
			push(@species, $sp);
		}
	}
	open(my $REPORT, '>report.tsv') or die("Cannot open '$o/report.tsv'!\n");
	print($REPORT "accessions\tfilename\t" . join("\t", reverse(@criteria)));
	if($w == 1){
		print($REPORT "\t" . join(" (s)\t", reverse(@criteria)) . ' (s)');
	}
	print($REPORT "\n");
	makeCombinations($#species, $combination, $REPORT);
	close($REPORT);
############################### OPTIONS
} else {
	print(STDERR "\noops.pl a script to select the 'best' set of paralogs from a gene\nfamily file.\n\n");
	print(STDERR "USAGE:\toops.pl [ -c A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z ]\n");
	print(STDERR "\t[ -d ] -i input.fasta [ -m lg|wag ] [ -n int ] -o output [ -r raxml ]\n");
	print(STDERR "\t[ -s ] [ -t number ] [ -w ]\n\n");
	print(STDERR "OPTIONS:\n");
	print(STDERR "-c\tOptimality criteria(on) for paralogs selection (default all):\n");
	print(STDERR "\t\t" . join("\n\t\t", reverse(@criteria)) . "\n");
	print(STDERR "-d\tDo not keep output directory clean (default = keep clean).\n");
	print(STDERR "-i\tInput sequences in FASTA format (names should be '>genus_species#sequenceID').\n");
	print(STDERR "-m\tUse the LG or WAG model for FastTree calculations (default = $m; RAxML\n\tconducts model testing).\n");
	print(STDERR "-n\tLimit analysis to n combinations (default = all).\n");
	print(STDERR "-o\tOutput directory.\n");
	print(STDERR "-r\tRAxML executable name (default = $r).\n");
	print(STDERR "-s\tUse RAxML in place of the default FastTree (slower);\n\tbootstrap is used in place of FastTree’s SH–like local supports.\n");
	print(STDERR "-t\tNumber of threads (default = $t) for MAFFT (and RAxML if -s is used).\n");
	print(STDERR "-w\tReport analysis times (default = do not).\n\n");
}
############################### FOWLER–NOLL–VO 32-BIT HASH
sub fnv32a {
	my @data = split(//, $_[0]);
	my $hash = 0x811c9dc5;
	for(my $k = 0; $k <= $#data; $k++){
		$hash = $hash^ord($data[$k]);
		$hash = ($hash*0x01000193)%4294967296; ### 2**32
	}
	return(sprintf("0x%X", $hash));
}
############################### FORMAT FASTA
sub formatFASTA {
	my $output = '>' . $_[0];
	my @sequence = split(//, $_[1]);
	for(my $k = 0; $k <= $#sequence; $k++){
		if(($k % 80) == 0){
			$output .= "\n" . $sequence[$k];
		} else {
			$output .= $sequence[$k];
		}
	}
	$output .= "\n";
	return($output);
}
############################### MAKE COMBINATIONS
sub makeCombinations {
	my $i = $_[0]; ### @species index
	my $combination = $_[1]; ### growing string (csv of species#terminal)
	my $REPORT = $_[2]; ### file handle
	if($i < 0){
		$n--;
		if(($N == 0) || (($N == 1) && ($n >= 0))){
			$combination =~ s/,$//;
			print($REPORT scoreCombination($combination));
		}
		return(0);
	} else {
		my $c = '';
		my @terminals = sort({rand() <=> 0.5} keys(%{$sequenceData->{$species[$i]}}));
		for(my $k = $#terminals; $k >= 0; $k--){
			$c .= makeCombinations($i-1, $combination . $species[$i] . '#' . $terminals[$k] . ',', $REPORT);
		}
		return($c);
	}
}
############################### ONLINE MEAN/VARIANCE (WELFORD 1962)
sub mean {
	my $M2 = $_[0];
	my $count = $_[1]+1;
	my $mean = $_[2];
	my $new = $_[3];
	my $nextMean = $mean+($new-$mean)/$count;
	$M2 = $M2+($new-$mean)*($new-$nextMean);
	return(($M2, $count, $nextMean));
}
############################### SCORE COMPINATION
sub scoreCombination {
	my @combination = sort({$b cmp $a} split(/,/, $_[0]));
	my $alignmentTime = 0;
	my $buffer = join('|', reverse(@combination));
	my $combinationHash = fnv32a($buffer);
	$buffer .= "\t$combinationHash";
	my $currentScores = {}; ### criteria -> score
	my $currentTimes = {}; ### criteria -> time
	for(my $k = $#criteria; $k >= 0; $k--){
		$currentScores->{$criteria[$k]} = 0;
		$currentTimes->{$criteria[$k]} = '—';
	}
	my $M2 = 0;
	my $count = 0;
	my $mean = 0;
	my $random = int(rand(100000));
	############################### SEQUENCE LENGTH
	my $start = time();
	open(OUTFILE, ">$combinationHash.fa") or die("Cannot open '$o/$combinationHash.fa'!\n");
	for(my $k = $#combination; $k >= 0; $k--){
		my ($species, $terminal) = split(/#/, $combination[$k]);
		($M2, $count, $mean) = mean($M2, $count, $mean, length($sequenceData->{$species}->{$terminal}));
		print(OUTFILE formatFASTA($combination[$k], $sequenceData->{$species}->{$terminal}));
	}
	close(OUTFILE);
	$currentScores->{'(A) mean sequence length'} = int($mean+0.5);
	$currentTimes->{'(A) mean sequence length'} = sprintf('%.2f', time()-$start);
	if($count > 1){
		$currentScores->{'(B) variance of sequence length'} = int(($M2/($count-1))+0.5);
		$currentTimes->{'(B) variance of sequence length'} = sprintf('%.2f', time()-$start);
	}
	if($c->{'(C) relative sequence length (s.d. above/below mean)'} == 1){
		if($globalSD != 0){
			if($mean > $globalMean){
				$currentScores->{'(C) relative sequence length (s.d. above/below mean)'} = '+' . sprintf('%.2f', (($mean-$globalMean)/$globalSD));
			} else {
				$currentScores->{'(C) relative sequence length (s.d. above/below mean)'} = '-' . sprintf('%.2f', (($globalMean-$mean)/$globalSD));
			}
			$currentTimes->{'(C) relative sequence length (s.d. above/below mean)'} = sprintf('%.2f', time()-$start);
		}
	}
	############################### RANDOM
	if($c->{'(D) random'} == 1){
		$start = time();
		$currentScores->{'(D) random'} = $random;
		$currentTimes->{'(D) random'} = sprintf('%.2f', time()-$start);
	}
	############################### RUN MAFFT
	my $alignedSequenceData = {};
	my $id = '';
	if(($c->{'(E) align score'} == 1) || ($c->{'(F) alignment length'} == 1) || ($c->{'(G) mean edit distance'} == 1) || ($c->{'(H) variance of edit distance'} == 1) || ($c->{'(I) sum of pairs'} == 1) || ($c->{'(J) sum of pairs BLOSUM45'} == 1) || ($c->{'(K) sum of pairs BLOSUM50'} == 1) || ($c->{'(L) sum of pairs BLOSUM62'} == 1) || ($c->{'(M) TCS (transitive consistency score)'} == 1) || ($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1) || ($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1) || ($c->{'(W) ML tree score'} == 1) || ($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1)){
		$start = time();
		qx`mafft --thread $t --auto --quiet $combinationHash.fa > $combinationHash.aln`;
		$alignmentTime = time()-$start;
		open(INFILE, "$combinationHash.aln") or die("Cannot open '$o/$combinationHash.aln'!\n");
		while(my $line = <INFILE>){
			chomp($line);
			if($line =~ m/^>(.+)/){
				$id = $1;
			} else {
				$line =~ s/-/?/g;
				$alignedSequenceData->{$id} .= $line;
			}
		}
		close(INFILE);
	}
	if(!length($id)){
		die('Alignment failed!');
	}
	my @terminals = sort({$b cmp $a} keys(%{$alignedSequenceData}));
	my $nodes = scalar(@terminals)-2;
	############################### ALIGN SCORE
	if($c->{'(E) align score'} == 1){
		$start = time();
		$currentScores->{'(E) align score'} = qx`alignAA.pl -i $combinationHash.aln`;
		chomp($currentScores->{'(E) align score'});
		$currentTimes->{'(E) align score'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	############################### ALIGNMENT LENGTH
	if(($c->{'(F) alignment length'} == 1) || ($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1) || ($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1) || ($c->{'(W) ML tree score'} == 1) || ($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1)){
		$start = time();
		$currentScores->{'(F) alignment length'} = length($alignedSequenceData->{$id});
		$currentTimes->{'(F) alignment length'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	############################### EDIT DISTANCE
	if(($c->{'(G) mean edit distance'} == 1) || ($c->{'(H) variance of edit distance'} == 1)){
		$start = time();
		($currentScores->{'(G) mean edit distance'}, $currentScores->{'(H) variance of edit distance'}) = split(/ /, qx`editAA.pl -i $combinationHash.aln`);
		chomp($currentScores->{'(H) variance of edit distance'});
		$currentTimes->{'(G) mean edit distance'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
		$currentTimes->{'(H) variance of edit distance'} = $currentTimes->{'(G) mean edit distance'};
	}
	############################### SUM OF PAIRS
	if($c->{'(I) sum of pairs'} == 1){
		$start = time();
		$currentScores->{'(I) sum of pairs'} = qx`pairsAA.pl -g 1 -G 10 -m 1 -i $combinationHash.aln`;
		chomp($currentScores->{'(I) sum of pairs'});
		$currentTimes->{'(I) sum of pairs'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	if($c->{'(J) sum of pairs BLOSUM45'} == 1){
		$start = time();
		$currentScores->{'(J) sum of pairs BLOSUM45'} = qx`pairsAA.pl -m -45 -i $combinationHash.aln`;
		chomp($currentScores->{'(J) sum of pairs BLOSUM45'});
		$currentTimes->{'(J) sum of pairs BLOSUM45'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	if($c->{'(K) sum of pairs BLOSUM50'} == 1){
		$start = time();
		$currentScores->{'(K) sum of pairs BLOSUM50'} = qx`pairsAA.pl -m -50 -i $combinationHash.aln`;
		chomp($currentScores->{'(K) sum of pairs BLOSUM50'});
		$currentTimes->{'(K) sum of pairs BLOSUM50'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	if($c->{'(L) sum of pairs BLOSUM62'} == 1){
		$start = time();
		$currentScores->{'(L) sum of pairs BLOSUM62'} = qx`pairsAA.pl -m -62 -i $combinationHash.aln`;
		chomp($currentScores->{'(L) sum of pairs BLOSUM62'});
		$currentTimes->{'(L) sum of pairs BLOSUM62'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	############################### TCS
	if($c->{'(M) TCS (transitive consistency score)'} == 1){
		$start = time();
		qx`t_coffee -infile=../$o/$combinationHash.aln -mode=evaluate -output=score_ascii -n_core=1 -multi_core=no 2>&1 /dev/null`; ##### Directory traversal used to circumvent a t_coffee bug. Forced to only use one processor due to a different bug.
		open(INFILE, "$combinationHash.score_ascii") or die("Cannot open '$o/$combinationHash.score_ascii'!\n");
		while(my $line = <INFILE>){
			chomp($line);
			if($line=~m/^SCORE=([0-9]+)$/){
				$currentScores->{'(M) TCS (transitive consistency score)'} = $1;
				last;
			}
		}
		close(INFILE);
		$currentTimes->{'(M) TCS (transitive consistency score)'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
	}
	############################### TNT
	if(($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1) || ($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1)){
		$start = time();
		open(OUTFILE, ">$combinationHash.tnt") or die("Cannot open '$o/$combinationHash.tnt'!\n");
		print(OUTFILE "nstates 32\n");
		print(OUTFILE "xread\n");
		print(OUTFILE $currentScores->{'(F) alignment length'} . ' ' . scalar(@terminals) . "\n");
		print(OUTFILE "& [prot]\n");
		for(my $k = $#terminals; $k >= 0; $k--){
			print(OUTFILE "$terminals[$k]\t$alignedSequenceData->{$terminals[$k]}\n");
		}
		print(OUTFILE ";\n");
		close(OUTFILE);
		my $tnt = "tnt bground log $combinationHash.tnt.log, watch=, mxram 10000 p $combinationHash.tnt, hold 1000, xi, col3, rs $random, ";
		if(($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1)){
			$tnt .= "log $combinationHash.minmax.log, minmax, ";
		}
		if(($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1)){
			$tnt .= "log $combinationHash.beststeps.log, mu=rep20ho1rat, bbreak=tbr, ";
		}
		if(($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1)){
			$tnt .= "log $combinationHash.jac.log, mu=rep1h2spr, k0, resample jak rep 100 frequency [mu=rep1h2spr], ";
		}
		$tnt .= 'log/, quit';
		my $tntTime = (time()-$start)+$alignmentTime;
		qx`$tnt`;
		open(INFILE, "$combinationHash.tnt.log") or die("Cannot open '$o/$combinationHash.tnt.log'!\n");
		while(my $line = <INFILE>){
			chomp($line);
			if($line =~ m/^xread ([0-9]+\.[0-9][0-9]) secs\. $/){
				$tntTime += $1;
			} elsif($line =~ m/^xinact ([0-9]+\.[0-9][0-9]) secs\. $/){
				$tntTime += $1;
				last;
			}
		}
		close(INFILE);
		if(($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1)){
			$start = time();
			open(INFILE, "$combinationHash.minmax.log") or die("Cannot open '$o/$combinationHash.minmax.log'!\n");
			while(my $line = <INFILE>){
				chomp($line);
				if($line =~ m/^Minimum possible steps \(total = ([0-9]+)\)/){
					$currentScores->{'(O) minimum possible parsimony steps'} = $1;
				} elsif($line =~ m/^Maximum possible steps \(total = ([0-9]+)\)/){
					$currentScores->{'(N) maximum possible parsimony steps'} = $1;
					last;
				}
			}
			close(INFILE);
			$currentScores->{'(P) parsimony information'} = $currentScores->{'(N) maximum possible parsimony steps'}-$currentScores->{'(O) minimum possible parsimony steps'};
			$currentTimes->{'(N) maximum possible parsimony steps'} = sprintf('%.2f', $tntTime);
			$currentTimes->{'(O) minimum possible parsimony steps'} = $currentTimes->{'(N) maximum possible parsimony steps'};
			$currentTimes->{'(P) parsimony information'} = $currentTimes->{'(N) maximum possible parsimony steps'};
			$currentTimes->{'(R) ensemble CI'} = $currentTimes->{'(N) maximum possible parsimony steps'};
			$currentTimes->{'(S) ensemble RI'} = $currentTimes->{'(N) maximum possible parsimony steps'};
		}
		if(($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1)){
			$start = time();
			my $searchTime = $tntTime;
			open(INFILE, "$combinationHash.beststeps.log") or die("Cannot open '$o/$combinationHash.beststeps.log'!\n");
			while(my $line = <INFILE>){
				chomp($line);
				if($line =~ m/^Best score \(TBR\): ([0-9]+)/){
					$currentScores->{'(Q) parsimony tree length'} = $1;
				} elsif($line =~ m/^mult ([0-9]+\.[0-9][0-9]) secs\. $/){
					$searchTime += $1;
				} elsif($line =~ m/^bbreak ([0-9]+\.[0-9][0-9]) secs\. $/){
					$searchTime += $1;
					last;
				}
			}
			close(INFILE);
			$currentTimes->{'(Q) parsimony tree length'} = sprintf('%.2f', (time()-$start)+$searchTime);
			if($c->{'(R) ensemble CI'} == 1){
				if($currentScores->{'(Q) parsimony tree length'} == 0){
					$currentScores->{'(R) ensemble CI'} = '1.00';
				} else {
					$currentScores->{'(R) ensemble CI'} = sprintf('%.2f', ($currentScores->{'(O) minimum possible parsimony steps'}/$currentScores->{'(Q) parsimony tree length'}));
				}
				$currentTimes->{'(R) ensemble CI'} = $currentTimes->{'(Q) parsimony tree length'};
			}
			if($c->{'(S) ensemble RI'} == 1){
				if($currentScores->{'(Q) parsimony tree length'} == 0){
					$currentScores->{'(S) ensemble RI'} = '1.00';
				} else {
					$currentScores->{'(S) ensemble RI'} = sprintf('%.2f', (($currentScores->{'(N) maximum possible parsimony steps'}-$currentScores->{'(Q) parsimony tree length'})/$currentScores->{'(P) parsimony information'}));
				}
				$currentTimes->{'(S) ensemble RI'} = $currentTimes->{'(Q) parsimony tree length'};
			}
		}
		if(($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1)){
			$M2 = 0;
			$count = 0;
			my $j = 1;
			my $jacTime = $tntTime;
			$mean = 0;
			$start = time();
			open(INFILE, "$combinationHash.jac.log") or die("Cannot open '$o/$combinationHash.jac.log'!\n");
			while(my $line = <INFILE>){
				chomp($line);
				if($line =~ m/^GC values/){
					$j = 0;
				} elsif($line =~ m/^mult ([0-9]+\.[0-9][0-9]) secs\. $/){
					$jacTime += $1;
				} elsif($line =~ m/^([0-9]+\.[0-9][0-9]) secs\. to complete resampling $/){
					$jacTime += $1;
				} elsif($j == 1){
					my @jac = ($line =~ /[`,]-{2,3}(100|[0-9][0-9])(?:[|+]|-{2,3} )/g);
					for(my $k = $#jac; $k >= 0; $k--){
						($M2, $count, $mean) = mean($M2, $count, $mean, $jac[$k]);
						$currentScores->{'(T) normalized sum of parsimony jackknife'} += $jac[$k];
					}
				}
			}
			close(INFILE);
			$currentScores->{'(T) normalized sum of parsimony jackknife'} = int(($currentScores->{'(T) normalized sum of parsimony jackknife'}/$nodes)+0.5);
			$currentTimes->{'(T) normalized sum of parsimony jackknife'} = sprintf('%.2f', (time()-$start)+$jacTime);
			$currentScores->{'(U) mean parsimony jackknife'} = int($mean+0.5);
			$currentTimes->{'(U) mean parsimony jackknife'} = $currentTimes->{'(T) normalized sum of parsimony jackknife'};
			if($count > 1){
				$currentScores->{'(V) variance of parsimony jackknife'} = int(($M2/($count-1))+0.5);
				$currentTimes->{'(V) variance of parsimony jackknife'} = $currentTimes->{'(T) normalized sum of parsimony jackknife'};
			}
		}
	}
	############################### ML
	if(($c->{'(W) ML tree score'} == 1) || ($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1)){
		if($s == 0){ ### FastTree
			$start = time();
			open(OUTFILE, ">$combinationHash.FastTree") or die("Cannot open '$o/$combinationHash.FastTree'!\n");
			for(my $k = $#terminals; $k >= 0; $k--){
				my $x = $alignedSequenceData->{$terminals[$k]};
				$x =~ tr/?X/--/;
				print(OUTFILE formatFASTA($terminals[$k], $x));
			}
			close(OUTFILE);
			qx`FastTree -quiet -boot 100 -$m -log $combinationHash.FastTree.log -gamma -seed $random -out $combinationHash.FastTree.newick $combinationHash.FastTree`;
			if($c->{'(W) ML tree score'} == 1){
				open(INFILE, "$combinationHash.FastTree.log") or die("Cannot open '$combinationHash.FastTree.log'!\n");
				while(my $line = <INFILE>){
					chomp($line);
					if($line =~ m/^Gamma20LogLk\t(-[0-9]{1,}\.[0-9]{1,})\t/){
						$currentScores->{'(W) ML tree score'} = $1;
						$currentTimes->{'(W) ML tree score'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
						last;
					}
				}
				close(INFILE);
			}
			if(($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1)){
				$M2 = 0;
				$count = 0;
				$mean = 0;
				open(INFILE, "$combinationHash.FastTree.newick") or die("Cannot open '$combinationHash.FastTree.newick'!\n");
				while(my $line = <INFILE>){
					chomp($line);
					my @sh = ($line =~ /(1\.000|0\.[0-9]{3,3}):/g);
					for(my $k = $#sh; $k >= 0; $k--){
						($M2, $count, $mean) = mean($M2, $count, $mean, 100*$sh[$k]);
						$currentScores->{'(X) normalized sum of ML support'} += 100*$sh[$k];
					}
				}
				close(INFILE);
				$currentScores->{'(X) normalized sum of ML support'} = int(($currentScores->{'(X) normalized sum of ML support'}/$nodes)+0.5);
				$currentTimes->{'(X) normalized sum of ML support'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
				$currentScores->{'(Y) mean ML support'} = int($mean+0.5);
				$currentTimes->{'(Y) mean ML support'} = $currentTimes->{'(X) normalized sum of ML support'};
				if($count > 1){
					$currentScores->{'(Z) variance of ML support'} = int(($M2/($count-1))+0.5);
				} else {
					$currentScores->{'(Z) variance of ML support'} = 0;
				}
				$currentTimes->{'(Z) variance of ML support'} = $currentTimes->{'(X) normalized sum of ML support'};
			}
		} else { ### RAxML
			$start = time();
			open(OUTFILE, ">$combinationHash.phy") or die("Cannot open '$o/$combinationHash.phy'!\n");
			print(OUTFILE scalar(@terminals) . ' ' . $currentScores->{'(F) alignment length'} . "\n");
			for(my $k = $#terminals; $k >= 0; $k--){
				print(OUTFILE "$terminals[$k]\t$alignedSequenceData->{$terminals[$k]}\n");
			}
			print(OUTFILE "\n");
			close(OUTFILE);
			my $raxml = "$r -T $t -m $MLmodel -s $combinationHash.phy -n $combinationHash -p $random";
			if(($c->{'(X) normalized sum of ML support'} == 0) && ($c->{'(Y) mean ML support'} == 0) && ($c->{'(Z) variance of ML support'} == 0)){
				qx`$raxml -f d`;
			} else {
				$raxml .= " -# 100 -x $random";
				if($c->{'(W) ML tree score'} == 1){
					$raxml .= ' -f a';
				}
				qx`$raxml`;
				qx`$r -T $t -m $MLmodel -J MR -z RAxML_bootstrap.$combinationHash -n $combinationHash.boo`;
			}
			if($c->{'(W) ML tree score'} == 1){
				open(INFILE, "RAxML_info.$combinationHash") or die("Cannot open '$o/RAxML_info.$combinationHash'!\n");
				while(my $line = <INFILE>){
					chomp($line);
					if($line =~ m/^Final ML Optimization Likelihood: (-[0-9]+\.[0-9]+)/){
						$currentScores->{'(W) ML tree score'} = $1;
						$currentTimes->{'(W) ML tree score'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
						last;
					}
				}
				close(INFILE);
			}
			if(($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1)){
				$M2 = 0;
				$count = 0;
				$mean = 0;
				open(INFILE, "RAxML_MajorityRuleConsensusTree.$combinationHash.boo") or die("Cannot open '$o/RAxML_MajorityRuleConsensusTree.$combinationHash.boo'!\n");
				while(my $line = <INFILE>){
					chomp($line);
					my @boo = ($line =~ /\[([0-9]{2,3})\]/g);
					for(my $k = $#boo; $k >= 0; $k--){
						($M2, $count, $mean) = mean($M2, $count, $mean, $boo[$k]);
						$currentScores->{'(X) normalized sum of ML support'} += $boo[$k];
					}
				}
				close(INFILE);
				$currentScores->{'(X) normalized sum of ML support'} = int(($currentScores->{'(X) normalized sum of ML support'}/$nodes)+0.5);
				$currentTimes->{'(X) normalized sum of ML support'} = sprintf('%.2f', (time()-$start)+$alignmentTime);
				$currentScores->{'(Y) mean ML support'} = int($mean+0.5);
				$currentTimes->{'(Y) mean ML support'} = $currentTimes->{'(X) normalized sum of ML support'};
				$currentScores->{'(Z) variance of ML support'} = int(($M2/($count-1))+0.5);
				$currentTimes->{'(Z) variance of ML support'} = $currentTimes->{'(X) normalized sum of ML support'};
			}
		}
	}
	############################### OUTPUT AND STORE SCORES
	my $extension = '.fa';
	if(($c->{'(E) align score'} == 1) || ($c->{'(F) alignment length'} == 1) || ($c->{'(G) mean edit distance'} == 1) || ($c->{'(H) variance of edit distance'} == 1) || ($c->{'(I) sum of pairs'} == 1) || ($c->{'(J) sum of pairs BLOSUM45'} == 1) || ($c->{'(K) sum of pairs BLOSUM50'} == 1) || ($c->{'(L) sum of pairs BLOSUM62'} == 1) || ($c->{'(M) TCS (transitive consistency score)'} == 1) || ($c->{'(N) maximum possible parsimony steps'} == 1) || ($c->{'(O) minimum possible parsimony steps'} == 1) || ($c->{'(P) parsimony information'} == 1) || ($c->{'(Q) parsimony tree length'} == 1) || ($c->{'(R) ensemble CI'} == 1) || ($c->{'(S) ensemble RI'} == 1) || ($c->{'(T) normalized sum of parsimony jackknife'} == 1) || ($c->{'(U) mean parsimony jackknife'} == 1) || ($c->{'(V) variance of parsimony jackknife'} == 1) || ($c->{'(W) ML tree score'} == 1) || ($c->{'(X) normalized sum of ML support'} == 1) || ($c->{'(Y) mean ML support'} == 1) || ($c->{'(Z) variance of ML support'} == 1)){
		$extension = '.aln'
	}
	my $times = '';
	for(my $k = $#criteria; $k >= 0; $k--){
		if($c->{$criteria[$k]} == 0){
			$buffer .= "\t—";
			$times .= "\t—";
		} else {
			$buffer .= "\t" . $currentScores->{$criteria[$k]};
			$times .= "\t" . $currentTimes->{$criteria[$k]};
			if($currentScores->{$criteria[$k]} < $bestScores->{$criteria[$k]}->{'min'}){
				$bestScores->{$criteria[$k]}->{'min'} = $currentScores->{$criteria[$k]};
				copy($combinationHash . $extension, "$criteriaLetters[$k]-min-final.fasta");
			}
			if($currentScores->{$criteria[$k]} > $bestScores->{$criteria[$k]}->{'max'}){
				$bestScores->{$criteria[$k]}->{'max'} = $currentScores->{$criteria[$k]};
				copy($combinationHash . $extension, "$criteriaLetters[$k]-max-final.fasta");
			}
		}
	}
	if($d == 1){
		unlink(glob($combinationHash . '.*'));
		if($s == 1){
			unlink(glob('RAxML_*.' . $combinationHash . '*'));
		}
	}
	if($w == 1){
		return($buffer . $times . "\n");
	} else {
		return($buffer . "\n");
	}
}
exit(0);
