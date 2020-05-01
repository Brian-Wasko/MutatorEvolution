#!/bin/bash
set -e
#make subdirectories to keep things clean
mkdir tmp
mkdir discordant
mkdir split
REFPATH=~/Bioinformatics/genomes
SLICENUM=1
# export $REFPATH
strain_name='AH11304evoS-A5-2'
# 'AH11304evoS-A1-1
# AH11304evoS-A1-2
# AH11304evoS-A1-3
# AH11304evoS-A12-1
# AH11304evoS-A12-2
# AH11304evoS-A12-3
# AH11304evoS-A3-1
# AH11304evoS-A3-2
# AH11304evoS-A3-3
# AH11304evoS-A5-1
# AH11304evoS-A5-2
# AH11304evoS-A5-3
# AH11304evoS-A7-1
# AH11304evoS-A7-2
# AH11304evoS-A7-3
# AH11304evoS-A8-1
# AH11304evoS-A8-2
# AH11304evoS-A8-3
# AH11304evoS-A9-1
# AH11304evoS-A9-2
# AH11304evoS-A9-3
# AH11304evoS-B1-1
# AH11304evoS-B1-2
# AH11304evoS-B1-3
# AH11304evoS-B10-1
# AH11304evoS-B10-2
# AH11304evoS-B10-3
# AH11304evoS-B11-1
# AH11304evoS-B11-2
# AH11304evoS-B11-3
# AH11304evoS-B2-1
# AH11304evoS-B2-2
# AH11304evoS-B2-3
# AH11304evoS-B3-1
# AH11304evoS-B3-2
# AH11304evoS-B3-3
# AH11304evoS-B8-1
# AH11304evoS-B8-2
# AH11304evoS-B8-3
# AH11304evoS-B9-1
# AH11304evoS-B9-2
# AH11304evoS-B9-3
# AH11304evoS-C11-1
# AH11304evoS-C11-2
# AH11304evoS-C11-3
# AH11304evoS-C2-1
# AH11304evoS-C2-2
# AH11304evoS-C2-3
# AH11304evoS-C5-1
# AH11304evoS-C5-2
# AH11304evoS-C5-3
# AH11304evoS-C9-1
# AH11304evoS-C9-2
# AH11304evoS-C9-3
# AH11304evoS-D10-1
# AH11304evoS-D10-2
# AH11304evoS-D10-3
# AH11304evoS-D12-1
# AH11304evoS-D12-2
# AH11304evoS-D12-3
# AH11304evoS-D5-1
# AH11304evoS-D5-2
# AH11304evoS-D5-3
# AH11304evoS-D9-1
# AH11304evoS-D9-2
# AH11304evoS-D9-3
# AH11304evoS-E1-1
# AH11304evoS-E1-2
# AH11304evoS-E1-3
# AH11304evoS-E10-1
# AH11304evoS-E10-2
# AH11304evoS-E10-3
# AH11304evoS-E7-1
# AH11304evoS-E7-2
# AH11304evoS-E7-3
# AH11304evoS-F1-1
# AH11304evoS-F1-2
# AH11304evoS-F1-3
# AH11304evoS-F10-1
# AH11304evoS-F10-2
# AH11304evoS-F10-3
# AH11304evoS-F11-1
# AH11304evoS-F11-2
# AH11304evoS-F11-3
# AH11304evoS-F12-1
# AH11304evoS-F12-2
# AH11304evoS-F12-3
# AH11304evoS-F6-1
# AH11304evoS-F6-2
# AH11304evoS-F6-3
# AH11304evoS-F7-1
# AH11304evoS-F7-2
# AH11304evoS-F7-3
# AH11304evoS-G6-1
# AH11304evoS-G6-2
# AH11304evoS-G6-3
# AH11304evoS-G7-1
# AH11304evoS-G7-2
# AH11304evoS-G7-3
# AH11304evoS-G8-1
# AH11304evoS-G8-2
# AH11304evoS-G8-3
# AH11304evoS-H10-1
# AH11304evoS-H10-2
# AH11304evoS-H10-3
# AH11304evoS-H11-1
# AH11304evoS-H11-2
# AH11304evoS-H11-3
# AH11304evoS-H2-1
# AH11304evoS-H2-2
# AH11304evoS-H2-3
# AH11304evoS-H3-1
# AH11304evoS-H3-2
# AH11304evoS-H3-3
# AH11304evoS-H4-1
# AH11304evoS-H4-2
# AH11304evoS-H4-3
# AH11304evoS-H6-1
# AH11304evoS-H6-2
# AH11304evoS-H6-3
# AH11304evoS-H9-1
# AH11304evoS-H9-2
# AH11304evoS-H9-3'

#Alignment of Alan's strains against entire yeast genome  count%${slienum} is count modulo.  It refers to the remainder of division.

count=0
for strain in $strain_name; do
bwa mem -M $REFPATH/yeast/S288C-masked-genome.fasta AH11304fasta/${strain}_R1_001.fastq.gz AH11304fasta/${strain}_R2_001.fastq.gz |samblaster -M -d discordant/${strain}.disc.sam -s split/${strain}.split.sam | grep -v SA:Z |grep -v XA: | samtools view -Sb -q 10 - |samtools sort -o tmp/${strain}_pe.unique_sorted.bam - &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/picard.jar AddOrReplaceReadGroups INPUT=tmp/${strain}_pe.unique_sorted.bam OUTPUT=tmp/${strain}_pe.fixedhdr.bam RGID=1 RGLB=1 RGPL=illumina   RGPU=TTAATA RGSM=${strain} VALIDATION_STRINGENCY=LENIENT &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    samtools index tmp/${strain}_pe.fixedhdr.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T RealignerTargetCreator -I tmp/${strain}_pe.fixedhdr.bam -R ${REFPATH}/yeast/S288C-masked-genome.fasta -o ${strain}.intervals &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T IndelRealigner -R ${REFPATH}/yeast/S288C-masked-genome.fasta -I tmp/${strain}_pe.fixedhdr.bam -targetIntervals ${strain}.intervals -o tmp/${strain}_pe-realigned.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T LeftAlignIndels -R ${REFPATH}/yeast/S288C-masked-genome.fasta -I tmp/${strain}_pe-realigned.bam -o tmp/${strain}_pe-leftaligned.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T BaseRecalibrator -I tmp/${strain}_pe-leftaligned.bam -R ${REFPATH}/yeast/S288C-masked-genome.fasta -knownSites ${REFPATH}/yeast/BY4741-diploid_snp_sorted_final.vcf -o tmp/${strain}_pe-recalibrated.grp &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T PrintReads -R ${REFPATH}/yeast/S288C-masked-genome.fasta -I tmp/${strain}_pe-leftaligned.bam -BQSR tmp/${strain}_pe-recalibrated.grp -o ${strain}_pe-recalibrated.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
	samtools mpileup -Bf $REFPATH/yeast/S288C-masked-genome.fasta ${strain}_pe-recalibrated.bam > ${strain}_pe.mpileup &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -jar ~/Bioinformatics/programs/VarScan.v2.3.9.jar mpileup2snp ${strain}_pe.mpileup --min-coverage 18 --min-var-freq .22 --strand-filter 1 |python ~/Bioinformatics/programs/variant_deSNPer2013a.py -s ~/Bioinformatics/genomes/yeast/AH0401.snp > ${strain}.noSNP.var &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
	java -jar ~/Bioinformatics/programs/VarScan.v2.3.9.jar mpileup2indel ${strain}_pe.mpileup --min-coverage 18 --min-var-freq .3 --strand-filter 1 |python ~/Bioinformatics/programs/variant_deSNPer2013a.py -s ~/Bioinformatics/genomes/yeast/AH0401indel.snp > ${strain}.no-indelSNP.var &
	        let count+=1
	        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
  samtools view -Sb discordant/${strain}.disc.sam|samtools sort -o ${strain}.disc.sort.bam -
  samtools view -Sb split/${strain}.split.sam|samtools sort -o ${strain}.split.sort.bam -
  samtools index ${strain}.disc.sort.bam
  samtools index ${strain}.split.sort.bam
    let count+=1
    [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done
rm *.intervals
rm -r tmp/
rm -r discordant
rm -r split
