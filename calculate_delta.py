# -*- coding: utf-8 -*-
"""
Created on Tue May 25 15:43:40 2021

@author: Fatemeh
"""


import pandas as pd
import sys, os, re, glob
from matplotlib.dates import date2num
from find_mutations import get_mutations_mat
from fitness_coefficient_calcu import fitness_CI_calculation, clusters_counts, calculte_frequency
import config


def check_pairs_CI(fitness_df):
    """calculates delta for each pairs fitness values delta = F_11 - F_01 - F_10 - F_00"""
    delta = []
    for col in fitness_df.columns:
        fitness_values = fitness_df[col]
        delta_1 = fitness_values[3][1] - fitness_values[0][0] - fitness_values[1][0] - fitness_values[2][0]
        delta_0 = fitness_values[3][0] - fitness_values[0][1] - fitness_values[1][1] - fitness_values[2][1]
        delta.append((delta_0, delta_1))
    fitness_df.loc['delta'] = delta
    return fitness_df
    
    
def multiple_labeling(labels, id_name, time_period, delay):
    """calculates delta for each pair (each column) in file and save it in a csv file"""
    fitnesses = pd.DataFrame()
    for col in labels.columns:
        clusters_name, [dates, counts] = clusters_counts(pd.DataFrame(labels[col]), metadata_info)
        frequency = calculte_frequency(counts, dates, clusters_name, time_period, delay)
        frequency.index = date2num(pd.to_datetime(frequency.index))
        # print(frequency)
        # frequency.to_csv(file.split('.csv')[0] + col + '_freq.csv')
        allele_fitness = []
        for i in ['-1', '0', '1', '2']:
            allele_pairs = ['-1', '0', '1', '2']
            allele_pairs.remove(i)
            sum_frequency = frequency[allele_pairs].sum(axis=1)
            Y_N_allele_freq = pd.DataFrame({i: frequency[i],'N': sum_frequency})
            fitness_CI = fitness_CI_calculation(Y_N_allele_freq, iterations_no, time_period)
            # fitness_CI = (np.power(1 + (fitness_CI_[0] / b), a), np.power(1 + (fitness_CI_[1] / b), a))
            print("Pair ", col, "fitness values: ", fitness_CI, "for haplotype ", i)
            allele_fitness.append(fitness_CI)
        fitnesses[col] = allele_fitness
    fitnesses.index = ['-1', '0', '1', '2']
    check_pairs_CI(fitnesses).to_csv(file.split('.csv')[0] + '_delta.csv')
    return check_pairs_CI(fitnesses)
            

if __name__ == '__main__':
    """file: the fasta file that includes sequences for a country, ref_file: reference fasta file, and metadata_file: metadata 
    as a csv file that includes Collection date, Virus name and Accession ID"""
    file = sys.argv[1]  
    ref_file = sys.argv[2]
    metadata_file = sys.argv[3]
    iterations_no = 100
    time_period = 1
    delay = 0
    #here id_name can be 'Virus name' or 'Accession ID'
    id_name = 'Accession ID' #'Virus name'
    labels = get_mutations_mat(file, ref_file, id_name)
    metadata_info = pd.read_csv(metadata_file)
    metadata_info = metadata_info.set_index(id_name)
    multiple_labeling(labels, id_name, time_period, delay)
