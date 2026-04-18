from probab3.commands.common.formulas.general import *
from probab3.commands.common.data_classes import *
from probab3.commands.common.general_code import *
from probab3.commands.common.formulas.star_cluster import StarCluster, TwoComponentCluster, UnequalMassCluster 
from probab3.commands.post_analysis.get_sample_stats import get_evolution_stats, get_run_details, get_file_stats
from probab3.commands.plot.plotter import plot_cluster_params_stats

import copy
import click
import itertools
import pprint

pp = pprint.PrettyPrinter(indent=2)


def get_triple_td_stats(file_path, verbose=False, cluster_class=UnequalMassCluster):
    if verbose:
        click.echo(f"Getting manual exchange stats from file {file_path}.")
    
    evolutions_df = read_dict_from_file(file_path)
    run_details = get_run_details(file_path)
    cluster = cluster_class(Mc_stars=run_details.get('Mc'), rh_stars=run_details.get('rh'), rc_stars=run_details.get('rc'))

    triple_td_stats_ids = {}
    counted_evolutions = len(evolutions_df['evolution'])
    
    for ev_index, evolution in enumerate(evolutions_df['evolution']):
        ev_stats = get_evolution_stats(evolution, cluster)
        if len(ev_stats.get('ms_td_Ns', [])) > 0:
            if evolutions_df['end'][ev_index] in triple_td_stats_ids.keys():
                triple_td_stats_ids[evolutions_df['end'][ev_index]].append(ev_index)
            else:
                triple_td_stats_ids[evolutions_df['end'][ev_index]] = [ev_index]

    if verbose:
        click.echo(f"triple td stats:\n{triple_td_stats_ids}")
    
    return run_details, counted_evolutions, triple_td_stats_ids


def get_few_file_triple_td_stats(file_paths, verbose=False, allowed_final_states=ALL_FINAL_STATES, run_param='Mc', cluster_class=UnequalMassCluster):
    triple_td_stats = {}
    param_values = []
    evolutions_num = []

    for file_path in file_paths:
        click.echo(f"Getting stats from file {file_path} with cluster class {cluster_class}.")
        run_details, counted_evolutions, triple_td_stats_ids = get_triple_td_stats(file_path, verbose=verbose, cluster_class=cluster_class)

        if run_param == 'c':
            param_values.append(run_details.get('rh')/run_details.get('rc'))
        else:
            param_values.append(run_details.get(run_param))

        evolutions_num.append(counted_evolutions) 
        for final_state in allowed_final_states:
            if triple_td_stats.get(final_state):
                triple_td_stats[final_state].append(len(triple_td_stats_ids.get(final_state, []))/counted_evolutions)
            else:
                triple_td_stats[final_state] = [len(triple_td_stats_ids.get(final_state, []))/counted_evolutions]

    return triple_td_stats, param_values, evolutions_num


def save_triple_td_stats(file_paths, run_param='Mc', verbose=False, save_path=None, cluster_class=UnequalMassCluster):
    
    triple_td_stats, param_values, evolutions_num = get_few_file_triple_td_stats(file_paths=file_paths, 
                                                                         verbose=verbose, 
                                                                         run_param=run_param,
                                                                         cluster_class=cluster_class)
    plot_stats_kwargs = {
        'param_values': param_values,
        'triple_td_stats': triple_td_stats,
        'evolutions_per_file': evolutions_num
    }

    if save_path:
        click.secho(f"Saving stats {plot_stats_kwargs} \n to file {save_path}", fg="green")
        save_numpy_dict_to_json_file(save_path, plot_stats_kwargs)


def plot_triple_td_stats_from_file(file_path, run_param='Mc', plot_key='man_ex_stats'):
    stats = read_dict_from_json_file(file_path)
    param_values = stats.get('param_values')
    if run_param == 'Mc':
        param_values = [x/MSun for x in param_values] 
    plot_cluster_params_stats(param_values, stats[plot_key], run_param)

def plot_triple_td_stats_vs_all_from_file(triple_td_file_path, total_file_path, run_param='Mc'):
    triple_td_dict = read_dict_from_json_file(triple_td_file_path)
    total_dict = read_dict_from_json_file(total_file_path)
    
    triple_td_stats = triple_td_dict.get('man_ex_stats')
    total_stats = total_dict.get('params_stats')

    triple_td_param_values = triple_td_dict.get('param_values')
    total_param_values = total_dict.get('param_values')

    # Step 1: Extract the list using the key and get the sorted indices
    triple_td_sorted_indices = sorted(range(len(triple_td_param_values)), key=lambda i: triple_td_param_values[i])
    total_stats_sorted_indices = sorted(range(len(total_param_values)), key=lambda i: total_param_values[i])

    # Step 2: Create a new dictionary with sorted lists
    triple_td_sorted_stats = {k: [v[i] for i in triple_td_sorted_indices] for k, v in triple_td_stats.items()}
    total_sorted_stats = {k: [v[i] for i in total_stats_sorted_indices] for k, v in total_stats.items()}

    ratios_stats = {}
    for man_ex_key in triple_td_sorted_stats.keys():
        ratios_stats[man_ex_key] = [(man_ex_val / total_val if total_val else 0) for man_ex_val, total_val in zip(triple_td_sorted_stats[man_ex_key], total_sorted_stats[man_ex_key])]

    if run_param == 'Mc':
        triple_td_param_values = [x/MSun for x in triple_td_param_values]
    
    plot_cluster_params_stats(triple_td_param_values, ratios_stats, run_param)
