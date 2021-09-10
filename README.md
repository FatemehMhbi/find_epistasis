# find_epistasis
 
This script finds epistasis mutations and returns their fitness value as a 95 confidence interval. \
First it creates a mutation matrix by comparing each sequence to the referece file. To help the code run faster, the columns with less than 200 mutations (1's in the mutation matrix), all 1's columns, and the sequences with more than a threshold number of mutations (10 + average mutation number) are removed from the matrix. \
Then, it extracts all the mutation pairs with 00, 01, 10, and 11 haplotypes. The fitness value, F_H is calculated for each haplotype H for each pair and delta = F_11 - F_01 - F_10 - F_00 is returned for each mutation pair. \
Fitness calculation uses the frequencey of a haplotype over time. Using the daily frequency as lambda, random poisson samples are drawn 2000 time. The mean value of the random samples is used for fitness calculation, in order to reduce the sampling error. This is repeated 100 times and the 95 confidence interval of this 100 fitness values are returned. 


# How to run:
Type in the command line: \
python3.7 calculate_delta.py sequence_file.fasta referece_file.fasta metadata_file.csv\
Example: \
python3.7 calculate_delta.py korea.fasta s_gene_ref.fasta metadata_2021_07_29.csv\
Note metadata_2021_07_29 is in csv format \
metadata_file includes 'Virus name', 'Accession ID', and 'Collection date' which are extracted from metadata.tsv file downloaded from GISAID.

