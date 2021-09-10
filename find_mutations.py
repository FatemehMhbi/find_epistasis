# -*- coding: utf-8 -*-
"""
Created on Tue May  4 19:47:31 2021

@author: Fatemeh
"""


from Bio import SeqIO
import numpy as np
import pandas as pd
from statistics import mean
from utils import extract_strain_name, extract_EPI_id
import re, os, sys
import multiprocessing as mp
import config


def remove_noisy_rows_cols(df):
    """removes all the columns with less than or more than a threshold ones,
    also removes all the rows that have more than 10 + average mutations"""
    print('Removing noisy sequences and columns...')
    sum_rowwise = df.sum(axis=1)
    avg_mutation_no = mean(sum_rowwise)
    print('average number of mutations in each seq: ', avg_mutation_no)
    rows_to_rm = df.index[np.where(sum_rowwise > avg_mutation_no + 10)]
    df = df.drop(rows_to_rm, axis=0)
    seq_no = df.shape[0]
    sum_colwise = df.sum(axis=0)
    small_cols_to_rm = set(map(str, np.where(sum_colwise < 200)[0].tolist()))
    big_cols_to_rm = set(map(str, np.where(sum_colwise > seq_no - 200)[0].tolist()))
    cols_to_rm = list(set(small_cols_to_rm).union(big_cols_to_rm))
    cols_to_rm = list(map(int, cols_to_rm))
    df = df.drop(cols_to_rm, axis=1)
    return df


def find_potential_ep_pairs(i):
    print(i)
    epistatis_pairs = []
    epistatis_labelings = []
    subtract_matrix = config.mutation_matrix - config.mutation_matrix.shift(periods=-i, axis="columns", fill_value=0)
    subtract_matrix = subtract_matrix.iloc[:,:-i]
    counts_10 = subtract_matrix[subtract_matrix > 0].sum(0)
    counts_01 = abs(subtract_matrix[subtract_matrix < 0].sum(0))
    logical_and = (config.mutation_matrix & config.mutation_matrix.shift(periods=-i, axis="columns", fill_value=0)).iloc[:,:-i]
    counts_11 = logical_and.sum()
    counts_df = pd.DataFrame({'10': counts_10,'01': counts_01, '11': counts_11})
    counts_df.index = zip(list(config.mutation_matrix.columns[:-i]), list(config.mutation_matrix.columns[i:]))
    pairs_with_larg_counts = counts_df[(counts_df > 100).all(1)]
    epistatis_pairs.append(pairs_with_larg_counts.index)
    labeling = subtract_matrix + 2 * logical_and
    labeling.columns = counts_df.index
    if len(pairs_with_larg_counts.index) != 0:
        epistatis_labelings.append(labeling[pairs_with_larg_counts.index])
    return epistatis_pairs, epistatis_labelings


def epist_mutations_from_mu_mat_parallel(mat_file):
    """checks all the pairs for 01, 10, 11 pairs"""
    a_pool = mp.Pool(processes = 10)
    epistatis_pairs_ = []
    epistatis_labelings_ = []
    for epistatis_pair, epistatis_labeling in a_pool.map(find_potential_ep_pairs, \
                                                          range(1,config.mutation_matrix.shape[1])):
        epistatis_pairs_.append(epistatis_pair)
        epistatis_labelings_.append(epistatis_labeling)
    epistatis_pairs_ = [item for elem in epistatis_pairs_ for item in elem]
    epistatis_labelings_ = [item for elem in epistatis_labelings_ for item in elem]
    epistatis_pairs_ = np.concatenate(epistatis_pairs_).ravel()
    epistatis_pairs_ = epistatis_pairs_[~pd.isnull(epistatis_pairs_)]
    epistatis_pairs = pd.DataFrame(epistatis_pairs_)
    epistatis_pairs.to_csv(str(mat_file).split(".fasta")[0] + '_potential_epistatis_pairs.csv')
    pairs_labels = pd.concat(epistatis_labelings_, axis=1)
    pairs_labels.to_csv(str(mat_file).split(".fasta")[0] + '_epistatis_pairs_labels.csv')
    return pairs_labels, epistatis_pairs



def epist_mutations_from_mu_mat():
    """checks all the pairs for 01, 10, 11 haplotypes"""
    print('Finding potential epistasis pairs...')
    mutation_matrix = config.mutation_matrix.fillna(0).astype(int)
    epistatis_pairs = []
    epistatis_labelings = []
    for i in range(1,mutation_matrix.shape[1]):
        subtract_matrix = mutation_matrix - mutation_matrix.shift(periods=-i, axis="columns", fill_value=0)
        subtract_matrix = subtract_matrix.iloc[:,:-i]
        counts_10 = subtract_matrix[subtract_matrix > 0].sum(0)
        counts_01 = abs(subtract_matrix[subtract_matrix < 0].sum(0))
        logical_and = (mutation_matrix & mutation_matrix.shift(periods=-i, axis="columns", fill_value=0)).iloc[:,:-i]
        counts_11 = logical_and.sum()
        counts_df = pd.DataFrame({'10': counts_10,'01': counts_01, '11': counts_11})
        counts_df.index = zip(list(mutation_matrix.columns[:-i]), list(mutation_matrix.columns[i:]))
        pairs_with_larg_counts = counts_df[(counts_df > 100).all(1)]
        epistatis_pairs.append(pairs_with_larg_counts.index)
        labeling = subtract_matrix + 2 * logical_and
        labeling.columns = counts_df.index
        if len(pairs_with_larg_counts.index) != 0:
            epistatis_labelings.append(labeling[pairs_with_larg_counts.index])
    epistatis_pairs = np.concatenate(epistatis_pairs).ravel()
    epistatis_pairs = epistatis_pairs[~pd.isnull(epistatis_pairs)]
    if epistatis_labelings:
        pair_labels = pd.concat(epistatis_labelings, axis=1)
    else:
        raise Exception('No epistasis mutations!')
    return pair_labels, pd.DataFrame(epistatis_pairs)


def get_mutations_mat(seq_file, ref_file, id_name):
    """makes a zero one mutation matrix and then checks all the pairs for 01, 10, 11 pairs, then save pairs 
    as a labeling file where -1 represent 01, 0 -> 00, 1 -> 10 and 2 -> 11"""
    print('Mutation matrix calculation...')
    ref = list(SeqIO.parse(ref_file, "fasta"))[0]
    sequences = list(SeqIO.parse(seq_file, "fasta"))
    mutation_list = []
    ids = []
    for seq in sequences:
        comparision = [int(i!=j) for i,j in zip(seq.seq, ref.seq)]
        mark_blanks = [int(i.isalpha()) for i in seq.seq]
        mutation_list.append([i*j for i,j in zip(comparision, mark_blanks)])
        ids.append(seq.id + seq.description)
        # print(seq.id + seq.description)
    if id_name == 'Virus name':
        ids = extract_strain_name(ids)
    elif id_name == 'Accession ID':
        ids = extract_EPI_id(ids)
    config.mutation_matrix = pd.DataFrame(mutation_list, index = ids).astype(int)
    # print(mutation_matrix)
    config.mutation_matrix.to_csv(str(seq_file).split(".fasta")[0] + '_mutation_mat.csv')
    config.mutation_matrix = remove_noisy_rows_cols(config.mutation_matrix)
    config.mutation_matrix.to_csv(seq_file.split(".fasta")[0] + '_cleaned.csv')
    if config.mutation_matrix.shape[1] > 50:
        [pair_labels, pairs] = epist_mutations_from_mu_mat_parallel(seq_file)
    else:
        [pair_labels, pairs] = epist_mutations_from_mu_mat()
    pairs.to_csv(seq_file.split(".fasta")[0] + '_potential_epistatis_pairs.csv')
    pair_labels.to_csv(seq_file.split(".fasta")[0] + '_epistatis_pairs_labels.csv')
    return pair_labels
    
