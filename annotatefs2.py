"""
annotate.py

takes input of BED formatted file with mutation locations (from mutation_caller.py) and annotates SNPs

python annotate.py -f <BED file containing mutations> -s <ORF sequences> -n <non-coding GFF file for annotations> -g <Genome sequence as FASTA>

"""

def main(BED, orfs, noncoding_file, genome_file):

	"""GET DATA TOGETHER"""

	#gather data about genes in a dictionary:
	#orf id, start, stop (both 5'-3'), exons, introns, protein length, chr
	genes = {}
	for record in SeqIO.parse( open(orfs, 'r'), 'fasta' ):
		start = 0
		stop = 0
		exons = map(lambda x: x.split('-'), record.description.split(', ')[1].split(' ')[-1].split(','))
		ch = chromosome_conversion(record.description.split('Chr ')[1].split(' ')[0])
		genename = record.description.split(' ')[1]
		category = record.description.split("\"")[0].split(', ')[-2]
		try:
			description = record.description.split("\"")[1]
		except IndexError:
			description = 'no record'

		#find introns in genes that have them
		introns = []
		if record.id.split('-')[0][-1] == 'C':
			exons = exons[::-1]
			start = int(exons[0][0])
			stop = int(exons[-1][1])
			if len(exons) > 1:
				for e in range(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		else:
			start = int(exons[0][0])
			stop = int(exons[-1][1])
			if len(exons) > 1:
				for e in range(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])

		genes[record.id] = [record.id] + [start, stop] + [exons] + [introns] + [len(record.seq)/3.0] + [ch] + [genename] + [category] + [description]


	#populate dictionary of chromosomes in genome
	genome = {}
	for record in SeqIO.parse( open(genome_file, 'r'), 'fasta' ):
		genome[chromosome_conversion(record.description.split(' [')[4].split('=')[1][0:-1])] = str(record.seq)

	#populate second dictionary of non-coding annotations
	#noncoding[ID] = [ ID, chrom, regiontype, start, stop ]
	noncoding = {}
	for line in open(noncoding_file, 'r'):
		if line.startswith('#') != True:
			l = line.strip().split('\t')
			noncoding[randint(1, 9999999)] = [ l[8].split(';')[0].split('=')[1], chromosome_conversion(l[0][3:]), l[2], int(l[3]), int(l[4]) ]

	"""ANNOTATE SNPS"""

	#open output file
	f_out = open(BED+'.annotated', 'w')

	#start reading BED file
    #chr, start, stop, ref, obs
	for line in open(BED, 'r'):
		#if it's a header, then print out more header to the next line
		if 'start' in line:
			print >> f_out, line.strip() + '\t' + '\t'.join(['annotation', 'region', 'protein'])
			continue
		l = line.strip().split('\t')

		fs_site = int(l[1])

        #make copy of chromosome with mutation: e.g. agctact becomes [a,g,c,t,a,c,t]
		mut_chr = list(copy(genome[chromosome_conversion(l[0])]))

        #Insert nucleotides
		if l[4][0] == '+':
			insertion = list(l[4][1:])
			indel = len(l[4][1:])
			mut_chr = mut_chr[0:fs_site] + insertion + mut_chr[fs_site:]
        #Delete nucleotides
		elif l[4][0] == '-':
			indel = -len(l[4][1:])
			mut_chr = mut_chr[0:fs_site] + mut_chr[fs_site - indel:]
		annotation = False
		fs_seq = ""
		wt_seq = ""
		fs_annotation = False
		#loop through genes, trying to find one containing snp
		for g in genes:
			#if gene on correct chromosome
			if genes[g][6] == chromosome_conversion(l[0]):
				#CRICK STRAND
				if genes[g][0].split('-')[0][-1] == 'C':
					#found gene containing snp
					if fs_site <= genes[g][1] and fs_site >= genes[g][2]:
						#check if snp is in an intron:
						for intron in genes[g][4]:
							if fs_site <= intron[0] and fs_site >= intron[1]:
								#found an intronic snp - check if it's near a splice site
								if fs_site in range(intron[0]-2, intron[0]+1) or fs_site in range(intron[1],intron[1]+3):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
								#otherwise it's just intronic (boring!)
								else:
									print >> f_out, '\t'.join(l + ['intron', genes[g][0], 'NA'])

						#if the snp is within gene start-stop, but not in an intron, then it's in an exon
						#remove intronic sequences from gene
						mut_gene = []
						wt_gene = []

						for exon in genes[g][3]:
							#check first if the snp is near a splice site
							if fs_site in range(int(exon[0])-2,int(exon[0])+1) or fs_site in range(int(exon[1]),int(exon[1])+3):
								#if it's not in the start/stop regions (these aren't splice-sites)
								if fs_site not in range(genes[g][1]-2,genes[g][1]+1) and fs_site not in range(genes[g][2],genes[g][2]+3):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])

							mut_gene += list    (               complement(        mut_chr[     (int(exon[1]) + indel)-1 : (int(exon[0]) + indel)     ] [::-1]          )                      )
							wt_gene += list     (               complement(        genome[genes[g][6]] [ int(exon[1])-1 :   int(exon[0])                 ] [::-1]         )                       )
							mutgenestr = "".join(str(x) for x in mut_gene)
							mutgenestr = Seq(mutgenestr, IUPAC.unambiguous_dna)
							mutpro = str(mutgenestr.translate())
							wtgenestr = "".join(str(x) for x in wt_gene)
							wtgenestr = Seq(wtgenestr, IUPAC.unambiguous_dna)
							wtpro = str(wtgenestr.translate())

						#loop through codons, find mismatch
							count = 0
							for aa in range(0, len(mutpro), 1):
								if fs_annotation == False:
									if mutpro[aa] == wtpro[aa]:
										count += 1
									elif mutpro[aa] != wtpro[aa]:
										wt_seq = wtpro[aa]
										if mutpro[aa] == '*':
											fs_seq = mutpro[aa]
											print >> f_out, '\t'.join(l + ['frameshift', genes[g][0], wt_seq[0] + str(count+1) + ":" + fs_seq, genes[g][7], genes[g][8], genes[g][9]])
											fs_annotation = True
										else:
											fs_seq = mutpro.split('*')[0][aa:]
											print >> f_out, '\t'.join(l + ['frameshift', genes[g][0], wt_seq[0] + str(count+1) + ":" + fs_seq + "*", genes[g][7], genes[g][8], genes[g][9]])
											fs_annotation = True
									annotation = True
							annotation = True

					#if snp_pos isn't in a gene, check if it's upstream of a gene
					elif fs_site in range(genes[g][1]+1,genes[g][1]+201):
						print >> f_out, '\t'.join(l + ["5'-upstream", genes[g][0], 'NA'])
						annotation = True

				#WATSON STRAND
				else:
					#found gene containing snp
					if fs_site >= genes[g][1] and fs_site <= genes[g][2]:

						#check if snp is in an intron:
						for intron in genes[g][4]:
							if fs_site >= intron[0] and fs_site <= intron[1]:
								#found an intronic snp - check if it's splice site
								if fs_site in range(intron[0],intron[0]+3) or fs_site in range(intron[1]-2,intron[1]+1):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
								#otherwise, it's just intronic
								else:
									print >> f_out, '\t'.join(l + ['intron', genes[g][0], 'NA'])
								annotation = True

						#if the snp is within the gene start-stop, but not in an intron, then it's in an exon
						#remove intronic sequences
						mut_gene = []
						wt_gene = []

						for exon in genes[g][3]:
							#first, check if in a splice-site
							if fs_site in range(int(exon[0]),int(exon[0])+3) or fs_site in range(int(exon[1])-2,int(exon[1])+1):
								#if it's not the start/stop
								if fs_site not in range(genes[g][1],genes[g][1]+3) and fs_site not in range(genes[g][2]-2,genes[g][2]+1):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])

							mut_gene += list(mut_chr[ int(exon[0])-1:int(exon[1]) ])
							wt_gene += list(genome[genes[g][6]][int(exon[0])-1:int(exon[1])])
							mutgenestr = "".join(str(x) for x in mut_gene)
							mutgenestr = Seq(mutgenestr, IUPAC.unambiguous_dna)
							mutpro = str(mutgenestr.translate())
							wtgenestr = "".join(str(x) for x in wt_gene)
							wtgenestr = Seq(wtgenestr, IUPAC.unambiguous_dna)
							wtpro = str(wtgenestr.translate())
							# print(mutgenestr)

						#loop through codons, find mismatch
						count = 0
						for aa in range(0, len(mutpro), 1):
							if fs_annotation == False:
								if mutpro[aa] == wtpro[aa]:
									count += 1
								elif mutpro[aa] != wtpro[aa]:
									wt_seq = wtpro[aa]
									if mutpro[aa] == '*':
										fs_seq = mutpro[aa]
										print >> f_out, '\t'.join(l + ['frameshift', genes[g][0], wt_seq[0] + str(count+1) + ":" + fs_seq, genes[g][7], genes[g][8], genes[g][9]])
										fs_annotation = True
									else:
										fs_seq = mutpro.split('*')[0][aa:]
										print >> f_out, '\t'.join(l + ['frameshift', genes[g][0], wt_seq[0] + str(count+1) + ":" + fs_seq + "*", genes[g][7], genes[g][8], genes[g][9]])
										fs_annotation = True

								annotation = True

						annotation = True

					#if snp isn't in a gene, check if it's in the upstream region of a gene
					elif fs_site in range(genes[g][1]-200,genes[g][1]):
						print >> f_out, '\t'.join(l + ["5'-upstream", genes[g][0], 'NA'])
						annotation = True

		#after checking all genes for genic or non-coding snps, check non-coding elements
		if annotation == False:
			#if snp isn't upstream of a gene, check if it's in a non-coding region of the genome
			for n in noncoding:
				#if correct chromosome
				if noncoding[n][1] == chromosome_conversion(l[0]):
					if fs_site >= noncoding[n][3] and fs_site <= noncoding[n][4]:
						print >> f_out, '\t'.join(l + [noncoding[n][2], noncoding[n][0], 'NA'])
						annotation = True
			#if annotation is still false after check non-coding elements, then the snp is just intergenic
			if annotation == False:
				print >> f_out, '\t'.join(l + ['intergenic', 'NA', 'NA'])

	"""OTHER USEFUL FUNCTIONS"""

def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10,
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':17, 'mitochondrion':17}
	if chrom_number.startswith('chr'):
		chrom_number = chrom_number[3:]

	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:
		return chrom_conv[chrom_number]

def complement(base):
	compbase = []
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	for i in range(len(base)):
		compbase.append(comp[base[i].upper()])
	return ''.join(compbase)

if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	from copy import copy
	from random import randint

	parser = OptionParser()
	parser.add_option('-f', '--input', action = 'store', type = 'string', dest = 'inputfile', help = 'file with mutations')
	parser.add_option('-s', '--sequences', action = 'store', type = 'string', dest = 'sequences', help = 'fasta file of coding sequences')
	parser.add_option('-n', '--non-coding', action = 'store', type = 'string', dest = 'noncoding', help = 'gff file containing non-coding regions')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'genome', help = 'fasta file containing genome sequence')
	(option, args) = parser.parse_args()

	main(option.inputfile, option.sequences, option.noncoding, option.genome)
