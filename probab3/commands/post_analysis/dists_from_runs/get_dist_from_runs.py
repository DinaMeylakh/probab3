from probab3.commands.post_analysis.get_sample_stats import get_all_properties_in_file, get_all_props_split_ims_fs, get_props_from_files
from probab3.commands.common.constants import *

import numpy as np
import matplotlib.pyplot as plt

def get_x_y_of_dist(samples, bin_num=10, log=False):
    bins = bin_num
    min_sample = np.min(np.array(samples))
    max_sample = np.max(np.array(samples))
    if log:
        bins = np.logspace(np.log10(min_sample), np.log10(max_sample), bin_num)
    counts, bins = np.histogram(samples, bins=bins)
    x_values = [0.5*(bins[i]+bins[i+1]) for i in range(len(bins)-1)]
    bin_widths = np.array([bins[i+1]-bins[i] for i in range(len(bins)-1)])
    y_values = counts/(len(samples)*bin_widths[0])

    return x_values, y_values

def plot_dist_of_file(file_path, property_name, property_label, y_label, bin_num=10, legend_label=''):
    all_properties = get_all_properties_in_file(file_path)
    x_values, y_values = get_x_y_of_dist(all_properties[property_name], bin_num=bin_num)
    plt.plot(x_values, y_values, label=legend_label)
    plt.xlabel(property_label)
    plt.ylabel(y_label)
    plt.legend(loc='lower left', bbox_to_anchor=(1.05, 0.0))

def plot_dist_of_files(file_paths, property_name, property_label, y_label, bin_num=10, legend_label='', endings=ALL_FINAL_STATES):

    all_relevant_properties = get_props_from_files(file_paths, property_name, endings=endings)
    
    if property_name == 'all_E0s':
        x_values, y_values = get_x_y_of_dist((-1)*np.array(all_relevant_properties), bin_num=bin_num, log=True)
        x_values = np.array(x_values)
    else:
        x_values, y_values = get_x_y_of_dist(all_relevant_properties, bin_num=bin_num)

    plt.plot(x_values, y_values, label=legend_label)
    plt.xlabel(property_label)
    plt.ylabel(y_label)
    plt.legend()
    return x_values, y_values

def plot_fs_ims_dists_of_file(file_path, property_name, property_label, y_label, bin_num=10, legend_label='', color='blue'):
    all_props_fss, all_props_imss = get_all_props_split_ims_fs(file_path)

    x_values_fs, y_values_fs = get_x_y_of_dist(all_props_fss[property_name], bin_num=bin_num)
    x_values_ims, y_values_ims = get_x_y_of_dist(all_props_imss[property_name], bin_num=bin_num) 
    plt.plot(x_values_fs, y_values_fs, label=f"{legend_label} FS", color=color)
    plt.plot(x_values_ims, y_values_ims, label=f"{legend_label} IMS", linestyle='dashed', color=color)
    plt.xlabel(property_label)
    plt.ylabel(y_label)

def plot_2d_scatter_of_files(file_paths, property_name_x, property_name_y, x_label, y_label, legend_label='', endings=ALL_FINAL_STATES, ev_stat_index=-1):

    all_x_properties = []
    all_y_properties = [] 
    for file_path in file_paths:
        file_properties = get_all_properties_in_file(file_path, endings=endings, ev_stat_index=ev_stat_index)

        all_x_properties = all_x_properties + file_properties.get(property_name_x, []) 
        all_y_properties = all_y_properties + file_properties.get(property_name_y, []) 
    
    plt.scatter(all_x_properties, (-1)*np.array(all_y_properties), label=legend_label, marker='x')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.yscale('log')
    plt.legend()
    return all_x_properties, all_y_properties