from probab3.commands.common.formulas.general import *
from probab3.commands.common.data_classes import *
from probab3.commands.common.general_code import *
from probab3.commands.common.formulas.star_cluster import StarCluster, TwoComponentCluster 
from probab3.commands.common.formulas.inspiral import Bounds

import copy
import click
import itertools
import pprint

pp = pprint.PrettyPrinter(indent=2)

def get_pre_evolution_stats(pre_evolution: list):

    all_aBs = []
    all_qBs = []
    all_aBs_tidal_limit = []
    all_qBs_for_merger = []
    all_qBs_for_EM = []
    all_CBs = []
    all_eBs = []
    crude_times = []

    pre_evolution_data = list(pre_evolution.pop("MSs2BHs_evolution"))
    
    inital_binary_mt = pre_evolution_data[0]
    inital_binary_dbs_state = MetastableState(*inital_binary_mt)

    initial_ma = inital_binary_dbs_state.binary_state.get('ma').get('mass')
    initial_mb = inital_binary_dbs_state.binary_state.get('mb').get('mass')
    all_aBs.append(aB(initial_ma, initial_mb, inital_binary_dbs_state.binary_state.get('EB')))
    all_qBs.append(qB_from_EBLB(initial_ma, initial_mb, inital_binary_dbs_state.binary_state.get('EB'),
                                inital_binary_dbs_state.binary_state.get('LB')))
    all_eBs.append(eB_from_EBLB(initial_ma, initial_mb, inital_binary_dbs_state.binary_state.get('EB'),
                                inital_binary_dbs_state.binary_state.get('LB')))
    all_aBs_tidal_limit.append(aB(ma=initial_ma,
                                  mb=initial_mb,
                                  EB=inital_binary_dbs_state.EB_cluster_tides))
    all_qBs_for_merger.append(inital_binary_dbs_state.qB_merger)
    all_qBs_for_EM.append(inital_binary_dbs_state.qB_EM)

    for ev_dbs_state_json in pre_evolution_data[1:]:
        if type(ev_dbs_state_json) != dict:
            ev_dbs_state = MetastableState(*ev_dbs_state_json)
        else:
            ev_dbs_state = MetastableState(**ev_dbs_state_json)
        
        if ev_dbs_state.state != DBSEvolutionState.NEW_MS.value:
            ev_dbs_ma = ev_dbs_state.binary_state.get('ma').get('mass')
            ev_dbs_mb = ev_dbs_state.binary_state.get('mb').get('mass')
            all_aBs.append(aB(ev_dbs_ma, ev_dbs_mb, ev_dbs_state.binary_state.get('EB')))
            all_qBs.append(qB_from_EBLB(ev_dbs_ma, ev_dbs_mb, ev_dbs_state.binary_state.get('EB'),
                                                            ev_dbs_state.binary_state.get('LB')))
            all_eBs.append(eB_from_EBLB(ev_dbs_ma, ev_dbs_mb, ev_dbs_state.binary_state.get('EB'), 
                                        ev_dbs_state.binary_state.get('LB')))
            crude_times.append(ev_dbs_state.crude_time)
            all_CBs.append(ev_dbs_state.binary_state.get('CB'))
            all_aBs_tidal_limit.append(aB(ma=ev_dbs_state.binary_state.get('ma').get('mass'),
                                          mb=ev_dbs_state.binary_state.get('mb').get('mass'),
                                          EB=ev_dbs_state.EB_cluster_tides))
            all_qBs_for_merger.append(ev_dbs_state.qB_merger)
            all_qBs_for_EM.append(ev_dbs_state.qB_EM)

    return {"all_aBs": all_aBs,
            "all_qBs": all_qBs,
            "all_qBs_for_merger": all_qBs_for_merger,
            "all_qBs_for_EM": all_qBs_for_EM,
            "crude_times": crude_times,
            "all_CBs": all_CBs,
            "all_eBs": all_eBs, 
            "all_aBs_tidal_limit": all_aBs_tidal_limit}


def get_evolution_stats(evolution: list, cluster=None):
    all_aBs = []
    all_qBs = []
    fs_Ns = []
    ims_Ns = []
    ejected_Ns = []
    crude_times = []
    dbs_kicks = []
    all_qBs_for_merger = []
    all_qBs_for_EM = []
    all_E0s = []
    all_L0s = []
    all_sampled_CSs = []
    all_CBs = []
    CBs_from_teritiary = []
    all_eBs = []
    all_L0s_norm = []
    all_aBs_tidal_limit = []
    all_as_cluster_tides = []
    ms_td_Ns = []
    manual_exchange_Ns = []
    manual_exchange_merger_Ns = []
    manual_exchange_em_merger_Ns = []
    previous_state = None

    N = 0
    ending = evolution[-1]
    for ev_dbs_state_json in evolution[:-1]:
        if type(ev_dbs_state_json) != dict:
            ev_dbs_state = MetastableState(*ev_dbs_state_json)
        else:
            ev_dbs_state = MetastableState(**ev_dbs_state_json)

        if ev_dbs_state.state != DBSEvolutionState.NEW_MS.value:
            all_aBs.append(ev_dbs_state.aB)
            all_qBs.append(ev_dbs_state.qB)
            all_eBs.append(eB_from_aBqB(ev_dbs_state.aB, ev_dbs_state.qB))
            crude_times.append(ev_dbs_state.crude_time)
            all_CBs.append(ev_dbs_state.binary_state.get('CB'))
            if ev_dbs_state.state != DBSEvolutionState.TRIPLE_TD.value:
                all_aBs_tidal_limit.append(aB(ma=ev_dbs_state.binary_state.get('ma').get('mass'),
                                            mb=ev_dbs_state.binary_state.get('mb').get('mass'),
                                            EB=ev_dbs_state.EB_cluster_tides) if ev_dbs_state.EB_cluster_tides else None)
                all_as_cluster_tides.append(ev_dbs_state.a_s_cluster_tides)
                if ev_dbs_state.state != DBSEvolutionState.MANUAL_EXCHANGE.value:
                    all_qBs_for_merger.append(ev_dbs_state.qB_merger)
                    all_qBs_for_EM.append(ev_dbs_state.qB_EM)

            if ev_dbs_state.state == DBSEvolutionState.IMS.value:
                ims_Ns.append(N)
            
            if ev_dbs_state.state == DBSEvolutionState.FS.value:
                fs_Ns.append(N)
                dbs_kicks.append(ev_dbs_state.dbs_kick)
            
            if ev_dbs_state.state == DBSEvolutionState.EJECTED_FS.value:
                ejected_Ns.append(N)

            if ev_dbs_state.state == DBSEvolutionState.TRIPLE_TD.value:
                all_aBs_tidal_limit.append(aB(ma=previous_state.binary_state.get('ma').get('mass'),
                                            mb=previous_state.binary_state.get('mb').get('mass'),
                                            EB=previous_state.EB_cluster_tides))
                all_as_cluster_tides.append(previous_state.a_s_cluster_tides)
                all_qBs_for_merger.append(previous_state.qB_merger)
                all_qBs_for_EM.append(previous_state.qB_EM)
                ms_td_Ns.append(N)
            
            if ev_dbs_state.state == DBSEvolutionState.MANUAL_EXCHANGE.value:
                manual_exchange_Ns.append(N)
                if cluster:
                    T_ref = cluster.binary_single_scattering_time(
                        ma=ev_dbs_state.binary_state.get('ma').get('mass'), mb=ev_dbs_state.binary_state.get('mb').get('mass'), R=ev_dbs_state.aB)
                
                    EB_merger = Bounds.EB_GW_analytic(ma=ev_dbs_state.binary_state.get('ma').get('mass'), 
                                                  mb=ev_dbs_state.binary_state.get('mb').get('mass'), 
                                                  eB=eB_from_aBqB(ev_dbs_state.aB, ev_dbs_state.qB),
                                                  T_ref=T_ref)
                
                    EB_EM = Bounds.EB_EM(ma=ev_dbs_state.binary_state.get('ma').get('mass'), 
                                         mb=ev_dbs_state.binary_state.get('mb').get('mass'), 
                                         eB=eB_from_aBqB(ev_dbs_state.aB, ev_dbs_state.qB))

                    if ev_dbs_state.binary_state.get('EB') < EB_merger:
                        if ev_dbs_state.binary_state.get('EB') < EB_EM:
                            manual_exchange_em_merger_Ns.append(N)
                        else: 
                            manual_exchange_merger_Ns.append(N)
                    
                    qB_merge = qB_from_EBeB(ma=ev_dbs_state.binary_state.get('ma').get('mass'), mb=ev_dbs_state.binary_state.get('mb').get('mass'), 
                                            EB=EB_merger, eB=eB_from_aBqB(ev_dbs_state.aB, ev_dbs_state.qB))
                    qB_EM = qB_from_EBeB(ma=ev_dbs_state.binary_state.get('ma').get('mass'), mb=ev_dbs_state.binary_state.get('mb').get('mass'), 
                                         EB=EB_EM, eB=eB_from_aBqB(ev_dbs_state.aB, ev_dbs_state.qB))
                    all_qBs_for_merger.append(qB_merge)
                    all_qBs_for_EM.append(qB_EM)
            
            N = N + 1
        
        else:
            all_E0s.append(ev_dbs_state.dbs_state.get('E0'))
            sample_L0 = ev_dbs_state.dbs_state.get('L0')
            all_L0s.append(sample_L0)
            all_sampled_CSs.append(ev_dbs_state.teritiary_state.get('Cs') if ev_dbs_state.teritiary_state else None) 
            all_CBs.append(ev_dbs_state.teritiary_state.get('CBs') if ev_dbs_state.teritiary_state else None)
            CBs_from_teritiary.append(ev_dbs_state.teritiary_state.get('CBs') if ev_dbs_state.teritiary_state else None)
            if ev_dbs_state.binary_state:
                all_L0s_norm.append(sample_L0/L0_max(ma=ev_dbs_state.binary_state.get('ma').get('mass'),
                                                 mb=ev_dbs_state.binary_state.get('mb').get('mass'),
                                                 ms=ev_dbs_state.dbs_state.get('ms').get('mass'),
                                                 E0=ev_dbs_state.dbs_state.get('E0'),
                                                 alpha=2.5))
                all_as_cluster_tides.append(ev_dbs_state.a_s_cluster_tides)
        
        previous_state = ev_dbs_state
        
    
    return {"all_aBs": all_aBs,
            "all_qBs": all_qBs,
            "all_qBs_for_merger": all_qBs_for_merger,
            "all_qBs_for_EM": all_qBs_for_EM,
            "fs_Ns": fs_Ns,
            "ims_Ns": ims_Ns,
            "ejected_Ns": ejected_Ns,
            "crude_times": crude_times,
            "dbs_kicks": dbs_kicks,
            "all_E0s": all_E0s,
            "all_L0s": all_L0s,
            "all_sampled_CSs": all_sampled_CSs,
            "all_CBs": all_CBs,
            "CBs_from_teritiary": CBs_from_teritiary,
            "all_eBs": all_eBs, 
            "all_L0s_norm": all_L0s_norm,
            "all_aBs_tidal_limit": all_aBs_tidal_limit,
            "all_as_cluster_tides": all_as_cluster_tides,
            "ms_td_Ns": ms_td_Ns,
            "manual_exchange_Ns": manual_exchange_Ns, 
            "manual_exchange_merger_Ns": manual_exchange_merger_Ns,
            "manual_exchange_em_merger_Ns": manual_exchange_em_merger_Ns}

def get_all_properties_in_file(file_path, endings=ALL_FINAL_STATES, ev_stat_index=-1):
    evolutions_df = read_dict_from_file(file_path)
    all_properties = {}
    filtered_evolutions = evolutions_df[evolutions_df['end'].isin(endings)]

    for evolution in filtered_evolutions['evolution']:
        ev_stats = get_evolution_stats(evolution)
        for key in ev_stats.keys():
            if all_properties.get(key):
                all_properties[key] += ev_stats[key]
            else:
                all_properties[key] = ev_stats[key]
    
    return all_properties

def get_certain_property_by_indexes(ev_stats_in_key, all_properties, key, Ns = []):
    if key not in {"all_aBs", "all_qBs", "all_qBs_for_merger", "all_qBs_for_EM", "all_eBs"}:
        return all_properties
    relevant_ev_stats = ev_stats_in_key[Ns] if any(Ns) else ev_stats_in_key
    if any(all_properties.get(key, [])):
        all_properties[key] = np.concatenate((all_properties[key], relevant_ev_stats))
        return all_properties
    all_properties[key] = relevant_ev_stats
    return all_properties

def get_all_props_split_ims_fs(file_path):
    evolutions_df = read_dict_from_file(file_path)
    all_props_fss = {}
    all_props_imss = {}
    for evolution in evolutions_df['evolution']:
        ev_stats = get_evolution_stats(evolution)
        fs_Ns = np.array(ev_stats.pop('fs_Ns'))
        ims_Ns = np.array(ev_stats.pop('ims_Ns'))
        
        for key in ev_stats.keys():
            all_props_imss = get_certain_property_by_indexes(np.array(ev_stats[key]), all_props_imss, key, ims_Ns)
            all_props_fss = get_certain_property_by_indexes(np.array(ev_stats[key]), all_props_fss, key, fs_Ns)

    return all_props_fss, all_props_imss 

def get_props_from_files(file_paths, property_name, endings=ALL_FINAL_STATES, ev_stat_index=-1):
    all_relevant_properties = []
    for file_path in file_paths:
        file_properties = get_all_properties_in_file(file_path, endings=endings, ev_stat_index=ev_stat_index)
        all_relevant_properties = all_relevant_properties + file_properties[property_name]

    return all_relevant_properties

def get_run_details(file_path, evolution_num=0):
    evolutions_df = read_dict_from_file(file_path)
    run_details = dict(evolutions_df.iloc[evolution_num])
    run_details.pop('evolution', None)
    run_details.pop('pre_evolution', None)
    run_details['evolutions_in_file'] = evolutions_df.shape[0]
    return run_details

def get_evolution_from_file(file_path, evolution_num=0):
    evolutions_df = read_dict_from_file(file_path)
    return dict(evolutions_df.iloc[evolution_num])

def get_stats_errors(file_stats, num_of_runs):
    stats_errors = {}
    for stat_key in file_stats.keys():
        p = file_stats.get(stat_key)
        stats_errors[stat_key] = (p*(1-p)/num_of_runs)**(1/2)

    return stats_errors

def get_file_stats(file_path, verbose=False, allowed_final_states=ALL_FINAL_STATES):
    evolutions_df = read_dict_from_file(file_path)

    run_details = copy.deepcopy(dict(evolutions_df.iloc[0]))
    run_details.pop('evolution', None)
    run_details.pop('pre_evolution', None)

    if verbose:
        total_evolutions = evolutions_df.shape[0]
        click.echo(f"total evolutions in file {total_evolutions}")

    stats = {}
    total_counted_evolutions = 0

    for allowed_final_state in allowed_final_states:
        stats[allowed_final_state] = (evolutions_df['end'] == allowed_final_state).sum()
        total_counted_evolutions += stats[allowed_final_state]

        if verbose:
            final_state_indexes = list(evolutions_df[evolutions_df['end'] == allowed_final_state].index) 
            click.echo(f"found evolution end: {allowed_final_state} in indexes: {final_state_indexes}")

    if total_counted_evolutions > 0:
        for final_state in stats.keys():
            stats[final_state] = stats[final_state]/total_counted_evolutions
    
    for final_state in allowed_final_states:
        if not stats.get(final_state):
            stats[final_state] = 0

    return stats, total_counted_evolutions, run_details

def get_merger_only_file_stats(file_path, verbose=False):

    allowed_final_states = MERGER_ONLY_FINAL_STATES  

    return get_file_stats(file_path, verbose=verbose, allowed_final_states=allowed_final_states)


def get_few_file_stats(file_paths, verbose=False, allowed_final_states=ALL_FINAL_STATES, 
                       run_param='Mc'):
    all_stats = {}
    all_stats_errors = {}
    param_values = []
    evolutions_num = []

    for file_path in file_paths:
        click.echo(f"Getting stats from file {file_path}.")
        file_stats, counted_evolutions, run_details = get_file_stats(file_path, verbose=verbose, allowed_final_states=allowed_final_states)

        if run_param == 'c':
            param_values.append(run_details.get('rh')/run_details.get('rc'))
        else:
            param_values.append(run_details.get(run_param))

        file_stats_errors = get_stats_errors(file_stats, counted_evolutions)
        evolutions_num.append(counted_evolutions) 
        for final_state in file_stats.keys():
            if all_stats.get(final_state):
                all_stats[final_state].append(file_stats[final_state])
                all_stats_errors[final_state].append(file_stats_errors[final_state])
            else:
                all_stats[final_state] = [file_stats[final_state]]
                all_stats_errors[final_state] = [file_stats_errors[final_state]]

    return all_stats, all_stats_errors, param_values, evolutions_num

def renormalize_stats_data(params_stats, allowed_final_states, evolution_nums):
    
    relevant_stats = {}
    relevant_evolution_nums = np.array([]) 
    relevant_stats_errors = {}

    for final_state, final_state_stats in params_stats.items():
        if final_state in allowed_final_states:
            relevant_stats[final_state] = np.array(final_state_stats)*np.array(evolution_nums)
            if len(relevant_evolution_nums) == 0:
                relevant_evolution_nums = relevant_stats[final_state]
            else:
                relevant_evolution_nums = relevant_evolution_nums + relevant_stats[final_state]

    for final_state, final_state_stats in relevant_stats.items(): 
        relevant_stats[final_state] = relevant_stats[final_state]/relevant_evolution_nums 
        p = relevant_stats[final_state] 
        relevant_stats_errors[final_state] = np.sqrt((p*(1-p)/relevant_evolution_nums))

    return relevant_stats, relevant_stats_errors 
        
def get_cluster_param_from_files(file_paths, cluster_param):
    normalization = {
        'c': 1,
        'rc': parsec,
        'rh': parsec,
        'Mc': MSun,
        'h': 1
    }
    cluster_param_values = [] 
    dispersion_vels = []

    for file_path in file_paths:
        run_details = get_run_details(file_path)
        star_cluster = StarCluster(Mc=run_details.get('Mc'),
                                   rh=run_details.get('rh'),
                                   rc=run_details.get('rc'))
        dispersion_vels.append(star_cluster.vrms/1000)
        if cluster_param == 'c':
            cluster_param_value = run_details['rh']/run_details['rc']
        else:
            cluster_param_value = run_details[cluster_param]/normalization[cluster_param]
        cluster_param_values.append(cluster_param_value)

    return cluster_param_values, dispersion_vels 