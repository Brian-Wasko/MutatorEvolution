#!/usr/bin/env python

#written in python 2
#Modified from JLS lineage caller to compare 3 isolates from evolved cultures
from collections import defaultdict
from argparse import ArgumentParser
import pandas as pd


def main():

    parser = ArgumentParser()
    parser.add_argument("--strain_info", action='store', dest="strain_info", help='csv file consisting of culture number followed by isolate names to be compared', required=True)
    parser.add_argument("--min_depth", action='store', dest="min_depth", type=int, help="Minimum depth to evaluate site", default=18)
    o = parser.parse_args()

    strain_dict = defaultdict(lambda: 0)
    scored_sites = defaultdict(lambda: 0)
    var_dict = defaultdict(lambda: [])

    strain_info_file = open(o.strain_info, 'r')

    for line in strain_info_file:
# creates a function called splitline that takes each line, removes the carriage returns and new line character and splits on ","
        splitline = line.strip('\r\n').split(',')
# "int" ensures the 0 index is read as a number that serves as a key with the rest of the line serving as the value
        strain_dict.update({int(splitline[0]): splitline[1:]})

        isolate_1 = splitline[1]
        isolate_2 = splitline[2]
        isolate_3 = splitline[3]

    strain_info_file.close()
# creates a writeable file site_counts.txt that will be written to
    site_count_file = open("site_counts.txt", 'w')
# creates an empty list called isolate_list
    isolate_list = []
#This names the keys in each key-value pair as cultures as it iterates through the dictionary
    for culture in strain_dict.keys():
# prints the culture ID number to stout so that we see it the shell window
        print culture
# This defines the term "isolate" as the second index in the key - value pair. This portion of the script is scoring all of the nucleotides in isolate that have min_depth, which we define as 18-fold depth.
# It counts them with nuc_ctr and writes the position to the dictionary of scored_sites.
        for isolate in strain_dict[culture][:3]:
            isolate_list.append(isolate)
            nuc_ctr = 0
            pileup_file = open(isolate + '_pe.mpileup', 'r')

            for line in pileup_file:

                splitline = line.strip('\n').split('\t')
# int turns a text string into a numerical value.  If the value of index 3 of the line in pileup file is greater than the specified min depth (18) it adds 1 to nuc_ctr.
                if int(splitline[3]) > o.min_depth:
                    nuc_ctr += 1
# This adds the chromosome and position number as a key in scored_sites dictionary if it's not there already and adds one to the value (see lambda: 0).  This is how
# the dictionary keeps track of how many times it has observed a site.
                    scored_sites[splitline[0] + ":" + str(splitline[1])] += 1
# You need to close the pile_up file because a new pile_up file will be opened as the program iterates through the strain list.
            pileup_file.close()
# This writes the tally of the number of scorable bases from the isolate to the site_count file
            site_count_file.write(isolate + '\t' + str(nuc_ctr) + '\n')

    scored_site_file = open("scored_sites.txt", 'w')
    scored_site_ctr = 0

#From the scored_sites dictionary that kept track of sites with at least 18-fold coverage in one strain, we now assess how many times we counted it and compare it to the
#the number of strains in our sublineage list.  For example, if we have 12 strains in our isolate_list and we observed it 12 times than we can add the site to the shared site set.
    shared_set=set()
    for site in scored_sites.keys():
        if scored_sites[site] == len(isolate_list):
 	    shared_set.add(site)
            scored_site_ctr += 1
#Documents all sites that are shared by writing it to the scored_site_file.
            scored_site_file.write(site + '\n')

    site_count_file.write("# Scored Sites:\t" + str(scored_site_ctr))

    scored_site_file.close()
    site_count_file.close()

    total_var_set = set()
# # The commented out line below creates a dictionary where the key = chr (ref|NC_001133|) and the value is a set() of mutations in the form of an integer "117324"
#     total_var_dict = defaultdict(lambda:set())  This was used to create a position file (commented out below).
    for culture in strain_dict.keys():
# We create a temporary list to write all of the vars that are within the shared set of positions
        temp_list = []
        for isolate in strain_dict[culture]:

            var_file = open(isolate + ".noSNP.var", 'r')
            isolate_var_set = set()

            for line in var_file:
                splitline = line.strip('\n').split('\t')

                if splitline[0] + ":" + splitline[1] in shared_set and splitline[0] != "Chrom":
                    isolate_var_set.add(splitline[0] + ":" + splitline[1])

# Could add "  + ":" + splitline[2] + ">" + splitline[3] " above and below if we want to include mutational information.  Not sure if that would muck things up.
                    total_var_set.add(splitline[0] + ":" + splitline[1])
# # This creates a key-value pair in total_var_dict in which the key is the chromosome and the value is a set of mutated positions
#                     if splitline[1] not in total_var_dict[splitline[0]]:
#                         total_var_dict[splitline[0]].add(int(splitline[1]))
# This appends temp list with the var_set.  So temp_list turns out to be a list of sets from each strain in a culture.  There will be
            temp_list.append(isolate_var_set)
# This associates each culture with a list of sets in the dictionary.
        strain_dict[culture] = temp_list
## This section uses the tot_var_dict to create a position file which gives the chromosome and nucleotide position of each variant.  Not necessary for our purposes.
    # pos_file = open("position_file.txt", 'w')
    # for chrom in sorted(total_var_dict.keys()):
    #
    #     pos_file.write(chrom + '\n')
    #     x = sorted(list(total_var_dict[chrom]))
    #
    #     for pos in x:
    #         pos_file.write(str(pos) + '\n')
    #
    # pos_file.close()

# This makes a sorted list of all of the cultures
    culture_keys = sorted(strain_dict.keys())

# # This creates an empty dataframe
    df_master = pd.DataFrame()

# Turns the total var set into a sorted list
    total_var_list = sorted(list(total_var_set))

    for x in culture_keys:
# This names all columns in our comparison with the culture ID number
        shared_set_vars = strain_dict[x][0].intersection((strain_dict[x][1]),(strain_dict[x][2]))
        shared_vars = sorted(list(shared_set_vars))
        isolate_1_branch_vars = sorted(list(strain_dict[x][0] - shared_set_vars))
        isolate_1_all_vars = sorted(list(strain_dict[x][0]))
        isolate_1_unique_vars = sorted(list(strain_dict[x][0] - strain_dict[x][1] - strain_dict[x][2]))
        isolate_2_branch_vars = sorted(list(strain_dict[x][1] - shared_set_vars))
        isolate_2_all_vars = sorted(list(strain_dict[x][1]))
        isolate_2_unique_vars = sorted(list(strain_dict[x][1] - strain_dict[x][0] - strain_dict[x][2]))
        isolate_3_branch_vars = sorted(list(strain_dict[x][2] - shared_set_vars))
        isolate_3_all_vars = sorted(list(strain_dict[x][2]))
        isolate_3_unique_vars = sorted(list(strain_dict[x][2] - strain_dict[x][1] - strain_dict[x][0]))
# Unlike in EvolutionCompare3, these next three represent the variants only shared between the two strains
        isolate_1v2_vars = sorted(list((strain_dict[x][0].intersection(strain_dict[x][1])) - shared_set_vars))
        isolate_1v3_vars = sorted(list((strain_dict[x][0].intersection(strain_dict[x][2])) - shared_set_vars))
        isolate_2v3_vars = sorted(list((strain_dict[x][1].intersection(strain_dict[x][2])) - shared_set_vars))


# This highly redundant bit of code is because I couldn't think of a way of specifying which unique_vars file should be used when sorting a noSNP.var.  I don't think
# that python would recognize isolate_1, for instance, as a variable.  I'm sure there is a way to do this.
    # for isolate_1:

        var_file = open(isolate_1 + ".noSNP.var", 'r')
        unique_var_file = open(isolate_1 + "_uniqueSNP.var", "w")
        branch_var_file = open(isolate_1 + "_branchSNP.var", "w")
        shared_var_file = open(isolate_1 + "_trunkSNP.var", "w")
        oneVtwo_file = open(isolate_1 + "_1v2_SNP.var", "w")
        oneVthree_file = open(isolate_1 + "_1v3_SNP.var", "w")

        for line in var_file:
            splitline = line.strip('\n').split('\t')

            if splitline[0] + ":" + splitline[1] in isolate_1_unique_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                unique_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_1_branch_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                branch_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in shared_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                shared_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_1v2_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                oneVtwo_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_1v3_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                oneVthree_file.write('\t'.join(newline) + '\n')

        var_file.close()
        branch_var_file.close()
        unique_var_file.close()
        oneVtwo_file.close()
        oneVthree_file.close()

    # for isolate_2:

        var_file = open(isolate_2 + ".noSNP.var", 'r')
        unique_var_file = open(isolate_2 + "_uniqueSNP.var", "w")
        branch_var_file = open(isolate_2 + "_branchSNP.var", "w")
        shared_var_file = open(isolate_2 + "_trunkSNP.var", "w")
        oneVtwo_file = open(isolate_2 + "_1v2_SNP.var", "w")
        twoVthree_file = open(isolate_2 + "_2v3_SNP.var", "w")

        for line in var_file:
            splitline = line.strip('\n').split('\t')

            if splitline[0] + ":" + splitline[1] in isolate_2_unique_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                unique_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_2_branch_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                branch_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in shared_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                shared_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_1v2_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                oneVtwo_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_2v3_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                twoVthree_file.write('\t'.join(newline) + '\n')

        var_file.close()
        branch_var_file.close()
        unique_var_file.close()
        oneVtwo_file.close()
        twoVthree_file.close()

    # for isolate_3:

        var_file = open(isolate_3 + ".noSNP.var", 'r')
        unique_var_file = open(isolate_3 + "_uniqueSNP.var", "w")
        branch_var_file = open(isolate_3 + "_branchSNP.var", "w")
        shared_var_file = open(isolate_3 + "_trunkSNP.var", "w")
        oneVthree_file = open(isolate_3 + "_1v3_SNP.var", "w")
        twoVthree_file = open(isolate_3 + "_2v3_SNP.var", "w")

        for line in var_file:
            splitline = line.strip('\n').split('\t')

            if splitline[0] + ":" + splitline[1] in isolate_3_unique_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                unique_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_3_branch_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                branch_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in shared_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                shared_var_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_1v3_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                oneVthree_file.write('\t'.join(newline) + '\n')

            if splitline[0] + ":" + splitline[1] in isolate_2v3_vars and splitline[0] != "Chrom":
                newline = splitline[0:4]
                newline.append(str(float(splitline[4].split(":")[4].split("%")[0])/100))
                twoVthree_file.write('\t'.join(newline) + '\n')

        var_file.close()
        branch_var_file.close()
        unique_var_file.close()
        oneVthree_file.close()
        twoVthree_file.close()

        column_titles = ['1_allSNPs','1_branchSNPs','1_uniqueSNPs','2_allSNPs','2_branchSNPs','2_uniqueSNPs','3_allSNPs','3_branchSNPs','3_uniqueSNPs','1 v 2 only','1 v 3 only','2 v 3 only','trunkSNPs']

# This creates a dataframe to report all mutations from each strain found in the shared genome at different branch points.
        df = pd.DataFrame({'1_allSNPs' : pd.Series(isolate_1_all_vars,index=isolate_1_all_vars),
                            '1_branchSNPs' : pd.Series(isolate_1_branch_vars,index=isolate_1_branch_vars),
                            '1_uniqueSNPs' : pd.Series(isolate_1_unique_vars,index=isolate_1_unique_vars),
                            '2_allSNPs' : pd.Series(isolate_2_all_vars,index=isolate_2_all_vars),
                            '2_branchSNPs' : pd.Series(isolate_2_branch_vars,index=isolate_2_branch_vars),
                            '2_uniqueSNPs' : pd.Series(isolate_2_unique_vars,index=isolate_2_unique_vars),
                            '3_allSNPs' : pd.Series(isolate_3_all_vars,index=isolate_3_all_vars),
                            '3_branchSNPs' : pd.Series(isolate_3_branch_vars,index=isolate_3_branch_vars),
                            '3_uniqueSNPs' : pd.Series(isolate_3_unique_vars,index=isolate_3_unique_vars),
                            '1 v 2 only' : pd.Series(isolate_1v2_vars,index=isolate_1v2_vars),
                            '1 v 3 only' : pd.Series(isolate_1v3_vars,index=isolate_1v3_vars),
                            '2 v 3 only' : pd.Series(isolate_2v3_vars,index=isolate_2v3_vars),
                            'trunkSNPs' : pd.Series(shared_vars,index=shared_vars)}).reindex(total_var_list)

# This step is necessary because the columns aren't indexed in df in the order they are entered.
        df = df.reindex(columns = column_titles)
        df_master = pd.concat([df_master,df], axis=1)

    df_master.to_csv("test.csv")

if __name__ == "__main__":
    main()
