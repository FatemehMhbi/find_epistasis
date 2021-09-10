# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:34:46 2021

@author: Fatemeh
"""
import numpy as np
import pandas as pd
import os, re, sys, glob
from matplotlib.dates import date2num


def poisson_dist_sampling_mean(df, sample_num):
    """for each entry (value) in df, the mean value of random values drown from
    poisson distribution is replaced"""
    df_to_array = df.to_numpy()
    # print(df_to_array.shape)
    row = []
    for i in range(df_to_array.shape[0]):
        column = []
        for j in range(df_to_array.shape[1]):
            column.append(np.mean(np.random.poisson(df_to_array[i][j], sample_num)))
        row.append(column)
    row = pd.DataFrame(row, index = df.index)
    # row.to_csv("test/poisson_distr_"+str(sample_num)+".csv")
    return row


def extract_strain_name(labels):
    ids = []
    for label in labels:
        try:
            ids.append(label[label.find("/") + 1:].split('|')[0])
        except:
            ids.append(label)
    ids = np.asarray(ids)
    return ids


def extract_EPI_id(labels):
    ids = []
    for label in labels:
        try:
            ids.append(label[label.find("EPI"):].split('|')[0])
        except:
            ids.append(label)
    ids = np.asarray(ids)
    return ids


def remove_nans(list_1, list_2):
    """remove every corresponding element in list_2 that has nan value in list_1"""
    nans = np.isnan(list_1)
    nan_indecis = np.where(nans)[0]
    return np.delete(list_2, nan_indecis).tolist()


def read_meta(metadata, id_name):
    print("reading metadata...")
    tsv_read = pd.read_csv(metadata, sep='\t')
    metadata_info = tsv_read.set_index(id_name)
    metadata_info = metadata_info[~metadata_info.index.duplicated(keep='first')]
    return metadata_info[['Accession ID', 'Collection date', 'Virus name']]


def remove_non_val_dates(dates):
    """removes no valid dates from the list of string dates with format 0000-00-00"""
    valid_dates = []
    count = 0
    for date in dates:
        if len(date) == 10:
            try:
                valid = pd.to_datetime(date)
                valid_dates.append(date)
            except:
                count = count + 1
                continue
    # print(count, " non valid dates")
    return valid_dates


def generate_datapoints(dates, mode, time_period, start_delay):
    """from the start date (+ start_delay), it calculates all the dates 'time_period' number of days apart from 
    each other, and converts it to a list of numbers if it is datetime"""
    dates = np.sort(dates)
    start = dates[0]
    end = dates[len(dates)-1]
    timepoints = []
    if (mode == 'date'):
        start = pd.to_datetime(start)
        date = start + np.timedelta64(start_delay,'D')
    else:
        date = start + start_delay
    timepoints.append(start)
    while True:
        if (mode == 'date'):
            date =  date + np.timedelta64(time_period,'D')
        else:
            date = date + time_period
        timepoints.append(date)
        if(date >= end):
            break
    dates = timepoints
    if (mode == 'date'):
        timepoints = date2num(timepoints)
    return timepoints, dates
    