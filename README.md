# MutatorEvolution
This repository holds the following files used in the manuscript by Tracy et al., "Spontaneous polyploids and antimutators compete during the evolution of mutator cells".

1) eex_yeast_pileline.sh (A Unix shell pipeline used to guide whole genome sequencing alignments.)
2) S288C-masked-genome.fasta (the repeat-masked genome we used for the alignments.)
3) variant_deSNPer2013a.py (a program written in python 2, used in eex_yeast_pileline.sh to go through variant files and remove polymorphisms and mapping errors found in the AH0401.snp and AH0401indel.snp files).
4) AH0401.snp (a text file that contains SNPs segregating in AH0401 as well as recurrent mapping artifacts (designated MA).
5) AH0401indel.snp (a text file that contains indels segregating in AH0401 as well as recurrent mapping artifacts (designated MA).
6) EvolutionCompare4.py (a program written in python 2 that goes through single nucleotide variant lists from three independent subclones of an evolved culture and groups them into shared and unique variant lists.)
7) EvolutionCompare4Indel.py (a program written in python 2 that goes through indel variants lists from three independent subclones of an evolved culture and groups them into shared and unique variant lists.)
8) Sample input file for  




