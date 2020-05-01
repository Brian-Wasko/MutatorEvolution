# MutatorEvolution
This repository holds the following files used in the manuscript by Tracy et al., "Spontaneous polyploids and antimutators compete during the evolution of mutator cells".  All BAM files generated for this study can be downloaded from  https://www.ncbi.nlm.nih.gov/sra. (BioProject ID PRJNA629499).

1) eex_yeast_pileline.sh (A Unix shell pipeline used to guide whole genome sequencing alignments.)
2) S288C-masked-genome.fasta (the repeat-masked genome used for the alignments.)
3) variant_deSNPer2013a.py (a program written in python 2, used in eex_yeast_pileline.sh to go through variant files and remove polymorphisms and mapping errors found in the AH0401.snp and AH0401indel.snp files).
4) AH0401.snp (a text file that contains SNPs segregating in AH0401 as well as recurrent mapping artifacts (designated MA)).
5) AH0401indel.snp (a text file that contains indels segregating in AH0401 as well as recurrent mapping artifacts (designated MA)).
6) EvolutionCompare4.py (a program written in python 2 that goes through single nucleotide variant lists from three independent subclones of an evolved culture and groups them into shared and unique variant lists.)
7) EvolutionCompare4Indel.py (a program written in python 2 that goes through indel variants lists from three independent subclones of an evolved culture and groups them into shared and unique variant lists.)
8) StrainInput.csv (a Sample input file used for EvolutionCompare4.py and EvolutionCompare4Indel.py) 
9) annotatefs2.py (a modified version of Annotate 0.1 (Pashkova et al. 2013) that goes through bed files of indels and annotates the mutations, including the translated frameshift product.)
10) HypergeometicCosinSimTestPub.py (a program written in python 3 that utilizes the output of snv_spectrum https://github.com/aroth85/snv-spectrum/ as a starting point to test whether two spectra are statistically different using both the hypergeometric test and cosine similarity.)
11) Dataset 1 (Excel workbook that includes flow cytometry histograms reported in the paper, mutation rates, and simulation of the evolution of diploid mutator yeasts.  See Table of Contents for full details.).
12) Dataset 2 (Excel workbook that includes results from whole genome sequencing and hypergeometric and cosine similarity tests. See Table of Contents for full details.) 




