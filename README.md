# OOPS (Optimal Orthlogy for Phylogenomic Studies)

OOPS selects one sequence per species from a set of input paralogs using one of the 26 supported criteria ([see below](#use)). Multiple simultaneous selections can be conducted if more than one criterion is desired. The default is to use all criteria.

### install
Copy [oops.pl](https://github.com/dpl10/oops/blob/master/oops.pl) to a convient location in your PATH and ensure that a [Perl](https://www.perl.org/) version 5 interpreter is available (tested with 5.18).
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) version 7 must be in your PATH if criteria E–Z are employed.
* The script [alignAA.pl](https://github.com/dpl10/oops/blob/master/alignAA.pl) must be in your PATH if criterion E is employed.
* The script [editAA.pl](https://github.com/dpl10/oops/blob/master/editAA.pl) must be in your PATH if criteria G–H are employed.
* The script [pairsAA.pl](https://github.com/dpl10/oops/blob/master/pairsAA.pl) must be in your PATH if criteria I–L are employed.
* [T-Coffee](http://tcoffee.org/) version 11 must be in your PATH if criterion M is employed.
* [TNT](http://www.lillo.org.ar/phylogeny/tnt/) version 1.5 must be in your PATH if criteria N–V are employed.
* [RAxML](https://github.com/stamatak/standard-RAxML) version 8 or [FastTree](http://meta.microbesonline.org/fasttree/#Install) version 2 must be in your PATH if criteria W–Z are employed (the -s option uses RAxML in place of the default FastTree).

### input files
The input FASTA format file should contain amino acid sequences from multiple species and sequence names should be in the form of ‘>genus_species#sequenceID’. The number (hash) symbol should not be used in the taxon name or the sequence identification number.
```plaintext
>Paradonea_presleyi#D123
ELVISISDEAD
>Preseucoela_imallshookupis#L456
ELVISLIVES
```
### use
Options -i (input file) and -o (output directory) are required. The other options control the analyses parameters and optimality criteria used.

```plaintext
-c	Optimality criteria(on) for paralogs selection (default all):
		(A) mean sequence length
		(B) variance of sequence length
		(C) relative sequence length (s.d. above/below mean)
		(D) random
		(E) align score
		(F) alignment length
		(G) mean edit distance
		(H) variance of edit distance
		(I) sum of pairs
		(J) sum of pairs BLOSUM45
		(K) sum of pairs BLOSUM50
		(L) sum of pairs BLOSUM62
		(M) TCS (transitive consistency score)
		(N) maximum possible parsimony steps
		(O) minimum possible parsimony steps
		(P) parsimony information
		(Q) parsimony tree length
		(R) ensemble CI
		(S) ensemble RI
		(T) normalized sum of parsimony jackknife
		(U) mean parsimony jackknife
		(V) variance of parsimony jackknife
		(W) ML tree score
		(X) normalized sum of ML support
		(Y) mean ML support
		(Z) variance of ML support
-d	Do not keep output directory clean (default = keep clean).
-i	Input sequences in FASTA format (names should be '>genus_species#sequenceID').
-m	Use the LG or WAG model for FastTree calculations (default = lg; RAxML
	conducts model testing).
-n	Limit analysis to n combinations (default = all).
-o	Output directory.
-r	RAxML executable name (default = raxml).
-s	Use RAxML in place of the default FastTree (slower);
	bootstrap is used in place of FastTree’s SH–like local supports.
-t	Number of threads (default = 1) for MAFFT (and RAxML if -s is used).
-w	Report analysis times (default = do not).
```

### citation
If you use this software, please cite: Little, D.P., G. Eshel, A. Narechania, and M. Tessler. Submitted. OOPS: Optimal Orthlogy for Phylogenomic Studies. [Molecular Phylogenetics and Evolution](https://doi.org/ADD_DOI)

### license
[GPL](https://github.com/dpl10/oops/blob/master/LICENSE)
