from probab3.commands.common.formulas.general import *
from probab3.commands.common.data_classes import *
from probab3.commands.common.general_code import *
from probab3.commands.post_analysis.get_sample_stats import get_run_details, get_stats_errors 

from itertools import combinations_with_replacement
import click

def get_identity_file_stats(file_path, verbose=False, allowed_final_states=ALL_FINAL_STATES):
    evolutions_list = read_dict_from_file(file_path)
    total_evolutions = len(evolutions_list)
    if verbose:
        click.echo(f"total evolutions in file {total_evolutions}")

    stats = {}

    for binary_identity in list(combinations_with_replacement(OBJECT_IDENTITIES, 2)):
        stats[frozenset(binary_identity)] = {}

    total_counted_evolutions = 0
    for index, evolution in enumerate(evolutions_list):
        evolution_end_binary = evolution["evolution"][-2]
        end_ev_dbs_state = MetastableState(*evolution_end_binary) 
        binary_identity = frozenset({end_ev_dbs_state.binary_state.get('ma').get('identity'),
                                     end_ev_dbs_state.binary_state.get('mb').get('identity')})
        if (allowed_final_states and evolution["end"] in allowed_final_states):
            if stats[binary_identity].get(evolution["end"]):
                stats[binary_identity][evolution["end"]] = stats[binary_identity][evolution["end"]] + 1 
            else:
                stats[binary_identity][evolution["end"]] = 1

            total_counted_evolutions += 1
            if verbose:
                click.echo(f"evolution end: {evolution['end']} in index {index}")

    if total_counted_evolutions > 0:
        for binary_identity in stats.keys():
            for final_state in stats[binary_identity].keys():
                stats[binary_identity][final_state] = stats[binary_identity][final_state]/total_counted_evolutions


    for binary_identity in stats.keys():     
        for final_state in allowed_final_states:
            if not stats[binary_identity].get(final_state):
                stats[binary_identity][final_state] = 0

    return stats, total_counted_evolutions


def get_few_identity_file_stats(file_paths, verbose=False, allowed_final_states=ALL_FINAL_STATES, 
                       run_param='Mc'):
    all_stats = {}
    all_stats_errors = {}
    param_values = []

    for binary_identity in list(combinations_with_replacement(OBJECT_IDENTITIES, 2)):
        all_stats[frozenset(binary_identity)] = {}
        all_stats_errors[frozenset(binary_identity)] = {}

    for file_path in file_paths:
        run_details = get_run_details(file_path)
        if run_param == 'c':
            param_values.append(run_details.get('rc')/run_details.get('rh'))
        else:
            param_values.append(run_details.get(run_param))
        file_stats, counted_evolutions = get_identity_file_stats(file_path, verbose=verbose, allowed_final_states=allowed_final_states)
        for binary_identity in file_stats.keys():
            file_stats_errors = get_stats_errors(file_stats[binary_identity], counted_evolutions)
            for final_state in file_stats[binary_identity].keys():
                if all_stats[binary_identity].get(final_state):
                    all_stats[binary_identity][final_state].append(file_stats[binary_identity][final_state])
                    all_stats_errors[binary_identity][final_state].append(file_stats_errors[final_state])
                else:
                    all_stats[binary_identity][final_state] = [file_stats[binary_identity][final_state]]
                    all_stats_errors[binary_identity][final_state] = [file_stats_errors[final_state]]

    return all_stats, all_stats_errors, param_values