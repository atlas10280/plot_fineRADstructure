# Author: Matthew Bootsma
# Last edit: 2/01/2019
# this script will take a haps.vcf and print each individual's haplotype call on a line for each locus.
# should also function with a SNPs.vcf but not tested

#set working directory
import os
os.getcwd()
os.chdir("I:\\WAE_RAD_Data\\novoseq2\\SNPs\\Genomics\\FIS_filtered_work")
# open file you are going to write new results to
out_file = open("fineRAD_haps.txt", "w")
# open vcf file read all lines into an array, close file
raw_vcf_file = open("v2_genomic_haps.vcf", "r")
raw_vcf_array = raw_vcf_file.readlines()
raw_vcf_file.close()
# for each line in vcf file
for i in raw_vcf_array:
    # this is trying to recognize the header line
    if i.startswith("#CHROM"):
        # First, write column header for the locus ID
        out_file.write("locus")
        # then individual names as their own column headers
        header_line = i.rstrip().split("\t")
        # isolate the names
        names = header_line[9:(len(header_line)+1)]
        # write the names iteratively with a "\t" in between
        for j in names:
            out_file.write("\t" + j)
        # add eol after last name
        out_file.write("\n")
        # at this point we should have header line

    # use the if not # to go to the data
    elif "#" not in i:
        # split_ind_line is going to hold our genotype calls e.g., 0/2
        split_ind_line = i.rstrip().split("\t")
        # Store the locus ID in variable locus_ID
        locus_ID = split_ind_line[0]
        # Write that in the first column
        out_file.write(locus_ID)
                
        #Store the haplotype alleles in an array, they will be delimited by ","
        haplotypes = split_ind_line[3]+","+split_ind_line[4]
        #immediately break this haplotype allele array into it's components for indexing down the line
        haplotypes = haplotypes.split(",")
        
        #iterate through each genotype and convert it into the desired format
        # data starts in col 10 so -1 for base 0 index
        # len calculates the true number of entries from a base 1 index so -1 for base 0 index
        # k is the genotype reference call at that index, e.g. 0/1
        for k in split_ind_line[9:(len(header_line)+1)]:

            
            # handle missing data first
            if k.startswith("."):
                out_file.write("\t" + "")
            else:
            #break the allele calls down, encapsulate both alleles
                genotype_k = k.split("/")
                allele_k_index1 = int(genotype_k[0])
                allele_call_1 = haplotypes[allele_k_index1]
                allele_k_index2 = int(genotype_k[1])
                allele_call_2 = haplotypes[allele_k_index2]
                # if the allele calls are identical only write 1
                if allele_call_1 == allele_call_2:
                    out_file.write("\t" + allele_call_1)
                # if the allele calls are NOT identical, write both, delimit with a "/"
                else:
                    out_file.write("\t" + allele_call_1 + "/" + allele_call_2)
        # add eol before moving to next row (haplotype)
        out_file.write("\n")
        

out_file.close()  # script parse vcf file get alleles per individual
