# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 18:56:08 2020

@author: Fatemeh
"""

import pandas as pd
import numpy as np, scipy.stats as st
from matplotlib.dates import date2num
from scipy.interpolate import CubicSpline
from scipy import integrate
from utils import generate_datapoints, remove_nans, read_meta,\
    remove_non_val_dates, extract_strain_name, poisson_dist_sampling_mean

    
def calculate_sum_per_timepoints(list_, dates, timepoints):
    """calculates the sum of elements in the list_ for all dates between timepoints, where dates and timepoints are
    numbers (date_to_num() results)"""
    last_timepoint = timepoints[0] - 1 
    summation_list = []
    # for timepoint in timepoints[1:]:
    for timepoint in timepoints:
        indices = np.where((last_timepoint < dates) & (dates <= timepoint))[0]
        last_timepoint = timepoint
        if indices.size == 0:
            summation_list.append(np.nan)
        else:
            summation_list.append(np.sum(np.array(list_)[indices]))
    return summation_list


def get_counts_for_cluster(cluster_counts, cluster):
    """extracts the count for a specific cluster for all tuples (each day has a tuple) in cluster_counts and
    returns it as a list"""
    counts = []
    for i in range(len(cluster_counts)):
        try:
            counts.append([item[1] for item in cluster_counts[i] if str(item[0]) == cluster][0])
        except:
            counts.append(0) 
    return counts


def calculte_frequency(cluster_counts, dates, clusters, time_period, start_delay):
    """reads the data out of the csv file and returns the frequencies as dataframes, rows are weeks and columns are 
    clusters' frequency"""
    dates = pd.to_datetime(dates)
    time_p_indices, dates_generated = generate_datapoints(dates, 'date', time_period, start_delay)
    dates_indices = date2num(dates)

    cluster_freq = pd.DataFrame()
    
    for cluster in clusters:
        frequencies = get_counts_for_cluster(cluster_counts, str(cluster))
        weekly_cluster_freq = calculate_sum_per_timepoints(frequencies, dates_indices, time_p_indices)
        cluster_freq[str(cluster)] = weekly_cluster_freq
    cluster_freq.index = dates_generated  # time_p_indices[1:] 
    return cluster_freq


def find_daily_counts(metadata_info, all_submission_dates, label_name):
    """find counts of members of each cluster for each date"""
    df = pd.DataFrame()
    dates = []
    daily_cluster_counts = []
    for date in all_submission_dates:
        df = metadata_info.groupby('Collection date').get_group(date)
        count = df.groupby(label_name).count()
        dates.append(date)
        daily_cluster_counts.append(list(zip(list(count.index), list(count['Collection date'].values))))
    dates, daily_cluster_counts = zip(*list(sorted(zip(dates, daily_cluster_counts))))
    return dates, daily_cluster_counts



def clusters_counts(labels, metadata_info):
    """matches ids in labeling file and metadata file to find the submission date for each id, 
    then returns the daily counts of each cluster"""
    labels = labels[~labels.index.duplicated(keep='first')]
    label_name = labels.columns[0]
    clusters_name = labels[label_name].unique() #sorted(labels[label_name].unique())
    ixs = metadata_info.index.intersection(labels.index)
    metadata_info = metadata_info.loc[ixs]
    metadata_info = metadata_info[~metadata_info.index.duplicated(keep='first')]
    labels = labels.loc[ixs]
    # all_submission_dates = metadata_info['Collection date'].unique() 
    all_submission_dates = remove_non_val_dates(metadata_info['Collection date'].unique())
    # print(all_submission_dates)
    metadata_info[label_name] = labels.values
    return clusters_name, find_daily_counts(metadata_info, all_submission_dates, label_name)



def spline_interpolate(x,y,xs):
    """x is datapoints, y is values and xs is the datapoints we need interpolation for"""
    cspl =  CubicSpline(x, y)
    cspl_dr = cspl.derivative() 
    return cspl_dr(xs), cspl(xs)


def calculate_cubic_derivatives(df, time_period):
    """caluclates derivative of cubic spline interpolation over all timepoints for each column"""
    weeks_indices, blank = generate_datapoints(df.index, 'index', time_period, 0)
    # print(weeks_indices)
    df = df.replace(0, np.nan)
    derivation = pd.DataFrame(index = weeks_indices)
    interp = pd.DataFrame(index = weeks_indices)
    for col in df.columns:
        no_nan_indices = remove_nans(df[col], df.index)
        no_nan_values = df[col].loc[no_nan_indices]
        try:
            derivation[col], interp[col] = spline_interpolate(no_nan_indices, no_nan_values, weeks_indices)
        except:
            continue
    return derivation.loc[df.index], interp.loc[df.index]
    
    
def calculate_normal_freq(frequencies_df):
    """calculates the normalized frequencies for each row"""
    sum_df = frequencies_df.sum(axis=1)
    normal_freq = frequencies_df.div(sum_df.values, axis='index')
    return normal_freq


def integral(x, y):
    """return the integral for a set of points"""
    return integrate.cumtrapz(y, x, initial=0)


def integral_over_columns(column):
    """calculates the average fitness for the input column"""
    no_nan_indices = remove_nans(column, column.index)
    # print(no_nan_indices)
    time = range(1, len(no_nan_indices) + 1)
    coefficients = integral(time, column.loc[no_nan_indices])[-1]
    return coefficients / len(no_nan_indices)


def return_resistance_coef(cluster_freq_df, time_period, num_of_rows_to_ignore):
    """calculates fitness coefficient for each columns of df as a seperate cluster
    after spline interpolation, time_period is 1, 7 ,.. for daily, weekly,.. frequencies,
    num_of_rows_to_ignore is number of rows to ignore at the start of a cluster since start 
    day may introduce some bias"""
    fitness_function = pd.DataFrame(index = cluster_freq_df.index)
    resistance_coefficient = []
    h = pd.DataFrame(cluster_freq_df.sum(axis = 1), index = cluster_freq_df.index).cumsum()
    h_derivation, h_interp = calculate_cubic_derivatives(h, time_period)
    u = calculate_normal_freq(cluster_freq_df.cumsum())
    u_appr, u_interp =  calculate_cubic_derivatives(u, time_period)
    u_division = u_appr.div(u_interp)
    for col in u_division.columns:
        start = list(cluster_freq_df.index).index(cluster_freq_df[col].ne(0).idxmax())
        fitness_function[col] = u_division[col] + h_derivation.div(h)[0]
        # fitness_function[col].iloc[:start] = 0
        # resistance_coefficient.append(integral_over_columns(fitness_function[col]))
        resistance_coefficient.append(integral_over_columns(fitness_function[col].iloc[start + 1:]))
    return resistance_coefficient, fitness_function



def fitness_CI_calculation(frequency, iterations, time_period):
    """calculates fitness values and returns 95% confidence interval"""
    frequency = frequency.replace(np.nan, 0)
    # m = 5.2
    # s = 1.72
    # a = (m * m) / (s * s)
    # b = m / (s * s)
    fitness_for_clusters = pd.DataFrame(columns = frequency.columns)
    for i in range(iterations):
        distributed_freq = poisson_dist_sampling_mean(frequency, 2000)
        # print(distributed_freq)
        resistance_coef, fitness_function = return_resistance_coef(distributed_freq, time_period, 1)
        for column in fitness_for_clusters.columns:
            try:
                fitness_for_clusters.at[i, column] = resistance_coef[list(fitness_for_clusters.columns).index(column)]
            except:
                fitness_for_clusters.at[i, column] = 0
    # fitness_for_clusters = (1 + (fitness_for_clusters / b)).pow(a)
    relative_fitness = fitness_for_clusters[fitness_for_clusters.columns[0]] / \
              fitness_for_clusters[fitness_for_clusters.columns[1]]
    fitness_for_cluster_CI = st.t.interval(0.95, len(relative_fitness) - 1, loc = np.mean(relative_fitness), \
                                            scale = st.sem(relative_fitness))
    return fitness_for_cluster_CI

# if __name__ == '__main__':