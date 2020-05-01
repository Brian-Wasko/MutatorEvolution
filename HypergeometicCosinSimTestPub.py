"""Implementation of cosine similarity and hypergeometric (see Adams and Skopek) tests on mutation spectra
generated using snv_spectrum (https://github.com/aroth85/snv-spectrum/).  First written by AHerr & BHearn 1/16/19.
Modified in April 2020 by AHerr to add cosine similarity."""

import glob
import pyfaidx
import snv_spectrum as Snv
import sys
import pandas as pd
import numpy as np
import random
from collections import defaultdict
from matplotlib import pyplot
import math
from itertools import combinations
from scipy.spatial.distance import cosine

# This is the reference file used to determine the context of each mutation in variant_context.
ref_file='/Users/Mutator/Bioinformatics/genomes/yeast/S288C-masked-genome.fasta'

def variant_context(ref_file, chrom, pos, k=3):
    """Return the genomic context surrounding a position from an indexed FASTA
    file. Logically, the length of the context must be odd such that a center
    can be considered.
    Parameters
    ----------
    ref_file : str
        Filepath location to a FASTA file reference (indexed).
    chrom : str
        The chromosome for the given position.
    pos : int
        The position in the chromosome that the context will be centered on.
    """
    reference = pyfaidx.Fasta(ref_file)

    if k % 2 != 1:
        raise ValueError('DNA sequence must be of odd length')

    flank = (k - 1) / 2
    context = str(reference[chrom][int(pos - flank - 1):int(pos + flank)])

    return context
    reference.close()

"""This function calculates factorials and converts them to logs so that they are managable."""
def logfact(x):
    return (math.log(math.factorial(x)))

"""This calculates the hypergeometric p-value of a dataframe under the assumption that the data
from the columns was drawn from the same process."""
def pcalc(df):
    df_ct = df.sum(axis=0)
    df_rt = df.sum(axis=1)
    total = sum(df_ct)
    # Calculates the log factorial of each cell in the dataframe
    df_f = df.applymap(logfact)
    # Calcuates the log factorial of the marginal totals.
    df_rtf = df_rt.apply(logfact)
    df_ctf = df_ct.apply(logfact)
    # Calculates the hypergeometric probability of the table.
    denominator = sum(df_rtf) + sum(df_ctf)
    numerator = sum(df_f.sum()) + logfact(total)
    logp = denominator - numerator
    return logp
#     p = math.exp(logp)

def get_table_hypergeo(perm, row_totals_cpy):
    random.shuffle(perm)
    table = [[0] * 96, [0] * 96]
    j = 0
    for i in perm:
        while row_totals_cpy[j] == 0:
            j += 1
        table[i][j] += 1
        row_totals_cpy[j] -= 1
    return table

"""This generates random tables with the same marginal totals as the observed table and then calculates the hypergeometric
probability and cosine similarity of each one, storing them in a table that will then be used to infer the significance
of the original table."""
def random_df(df, n): # df = incoming dataframe; n = number of iterations
    df_ctob = df.sum(axis=0)
    df_rtob = df.sum(axis=1)
    stat_list=[]
    list_rt = df_rtob.tolist()
#     print(list_rt)
    list_ct = df_ctob.tolist()
#     print(list_ct)
    count = 0
    permutation = [0] * list_ct[0] + [1] * list_ct[1]  # To implement Agresti (1979) method
    while count < n:
        rand_table = get_table_hypergeo(permutation, list_rt.copy())
        df_rand = pd.DataFrame()
        df_rand[a] = rand_table[0]
        df_rand[b] = rand_table[1]
#         print(df_rand)
        cosin_sim = (1 - cosine(df_rand[a], df_rand[b]))
        logpt = pcalc(df_rand)
        stat_list.append((cosin_sim,logpt))
        count += 1
    return stat_list

"""This goes through a directory of var files, opens each var file that meet certain criteria, and determines the mutation spectra and context.
For our studies, var files were generated using Varscan 2.3.9. The essential elements of the format are as follows:
ref|NC_001134|	563832	C	T (Chr, position, Ref, and Var nucleotides in tab delimited format).  Chr names must match the names in the ref seq file (e.g. S288C-masked-genome.fasta)."""

df_ob = pd.DataFrame()
for varfile in glob.iglob('/Volumes/Jan2020seq/30-337971423/00_fastq/MutEvo/AH11304/AH11304evoS*/*1_trunkSNP.var'):
    subdir_split = varfile.split("/")[8].split('.var')[0] #pulls out relevant name of file
    spectrum = Snv.Spectrum(k=3, reference_notation="pyrimidine")
    open_file = open(varfile, 'r')
    count = 0
    for line in open_file:
        try:
            splitline = line.split("\t")
            spectrum[Snv.Snv(reference=splitline[2], alternate=splitline[3], context=variant_context(ref_file, splitline[0], int(splitline[1]), k=3)).with_pyrimidine_reference]+=1
            count+=1
        except KeyError:
            continue
#         print (subdir_split, count)
        sl = list(spectrum)
        vector = [sublist[1] for sublist in sl]
# This generates figures of each spectrum.
    fig, (ax_main, ax_cbar) = Snv.plot_spectrum(spectrum, kind='density')
    fig.savefig(f"/Volumes/External4/170429_Yeast/Lane2/AH164_evolution/MutationSpectra/{subdir_split}.pdf")
    open_file.close()
# This generates a large table composed of columns of the above spectra.
    df_ob[subdir_split] = vector

# Number of desired iterations
n = 10000
# This makes a list of all possible pairwise combinations in the dataframe df_ob.
column_list=("cosine_sim","cos_p_est","logp_cutoff","Hypergeo_p_est")
table_dict={}
comb = list(combinations(df_ob.columns,2))
for (a,b) in comb:
    counter = defaultdict(list)
    df_test = pd.DataFrame()
    df_test[a] = df_ob[a]
    df_test[b] = df_ob[b]
    logpdf = pcalc(df_test)
    cosine_sim = (1 - cosine(df_ob[a], df_ob[b]))
    print(f"Cosine similarity between {a} v {b} is {cosine_sim}")
    print(f"for {a} v {b} logp cutoff is {logpdf}")
    stat_list = random_df(df_test, n)
#     print(stat_list)
    lower = 0
    for v in stat_list:
        if v[0] < cosine_sim:
            lower += 1
    cos_p_est = lower/n
    if cos_p_est == 0:
        print("cos_p_est < ", 1/n)
    else:
        print("cos_p_est = ", cos_p_est)
    lower = 0
    for v in stat_list:
        if v[1] < logpdf:
            lower += 1
    Hypergeo_p_est = lower/n
    if Hypergeo_p_est == 0:
        print("Hypergeo_p_est < ", 1/n)
    else:
        print("Hypergeo_p_est = ", Hypergeo_p_est)
    table_dict[a,b]=(cosine_sim,cos_p_est,logpdf,Hypergeo_p_est)
report_df = pd.DataFrame.from_dict(table_dict, orient='index',columns = column_list)
report_df.to_csv(f'/Volumes/GoogleDrive/My Drive/Spontaneous Polyploidization/AH11304_evoS_compare/trunks.csv') #creates csv of comparisons
