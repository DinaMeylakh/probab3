from probab3.commands.common.constants import *
from probab3.commands.common.formulas.general import *
from probab3.commands.common.formulas.inspiral import Bounds 
from probab3.commands.common.formulas import triple_phase_space
from probab3.commands.plot.plotter import read_dict_from_file, read_dict_from_json_file
from probab3.commands.pre_calculation.detectors_sensitivity.detectors_sensitivity import detectors_sensitivity
from probab3.commands.common.data_classes import *
from probab3.commands.common.general_code import EnhancedJSONEncoder
import json

import numpy as np
from matplotlib import pyplot as plt

EM_SENSITIVITY_RANGE = {
    EMSensitivityParams.e_min.value: list(np.logspace(-2, -0.5, 10)),
    EMSensitivityParams.f_min.value: list(4*np.logspace(0, 1, 13))
}

MULT = 1

def is_merger_eccentric_ims(ma, mb, qB, e_min, f_min):
    if qB < Bounds.qB_EM(ma, mb, ef=e_min, f=f_min):
        return 1
    return 0

def is_merger_eccentric_fs(ma, mb, eB, EB, e_min, f_min):
    if EB < Bounds.EB_EM(ma=ma, mb=mb, eB=eB, ef=e_min, f=f_min):
        return 1
    return 0


def calc_e_min_sensitivty_ims(ma, mb, qB):
    counts_as_eccentric = []
    for e_min in EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]:
        if qB < Bounds.qB_EM(ma, mb, ef=e_min):
            counts_as_eccentric.append(1)
        else:
            counts_as_eccentric.append(0)
    
    return np.array(counts_as_eccentric) * MULT

def calc_f_min_sensitivty_ims(ma, mb, qB):
    counts_as_eccentric = []
    for f_min in EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]:
        if qB < Bounds.qB_EM(ma, mb, f=f_min):
            counts_as_eccentric.append(1)
        else:
            counts_as_eccentric.append(0)
    return np.array(counts_as_eccentric) * MULT

def calc_e_min_sensitivty_fs(ma, mb, eB, EB):
    counts_as_eccentric = []
    for e_min in EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]:
        if EB < Bounds.EB_EM(ma=ma, mb=mb, eB=eB, ef=e_min):
            counts_as_eccentric.append(1)
        else:
            counts_as_eccentric.append(0)
    
    return np.array(counts_as_eccentric) * MULT

def calc_f_min_sensitivty_fs(ma, mb, eB, EB):
    counts_as_eccentric = []
    for f_min in EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]:
        if EB < Bounds.EB_EM(ma=ma, mb=mb, eB=eB, f=f_min):
            counts_as_eccentric.append(1)
        else:
            counts_as_eccentric.append(0)
    
    return np.array(counts_as_eccentric) * MULT


def get_scramble(scramble_data):
    scramble = triple_phase_space.StateSample(ma=scramble_data["ma"], 
                                              mb=scramble_data["mb"], 
                                              EB=scramble_data["EB"], 
                                              LB=scramble_data["LB"], 
                                              CB=scramble_data["CB"])
    return scramble

def get_merger_stats_errors(stats, total_mergers):

    logger.info(f"Total mergers: {total_mergers}")
    errors = {}
    for percision_dependency_param in stats.keys():
        errors[percision_dependency_param] = []
        for p in list(stats[percision_dependency_param]):
            errors[percision_dependency_param].append((p*(1-p)/total_mergers)**(1/2))

    return errors


def get_sensitivity_stats(file_path, verbose=False):
    
    evolutions_df = read_dict_from_file(file_path)
    merged_evolutions = evolutions_df[evolutions_df["end"].isin(MERGER_ONLY_FINAL_STATES)]
    e_min_list_len = len(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value])
    f_min_list_len = len(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value])

    stats = {
        MetastableFinalState.IMS_EM_MERGER.value: {EMSensitivityParams.e_min.value: np.array([0]*e_min_list_len),
                                                   EMSensitivityParams.f_min.value: np.array([0]*f_min_list_len)},
        MetastableFinalState.IMS_MERGER.value: 0,
        MetastableFinalState.FS_EM_MERGER.value: {EMSensitivityParams.e_min.value: np.array([0]*e_min_list_len),
                                                   EMSensitivityParams.f_min.value: np.array([0]*f_min_list_len)},
        MetastableFinalState.FS_MERGER.value: 0,
        MetastableFinalState.EJECTED_FS_EM.value: {EMSensitivityParams.e_min.value: np.array([0]*e_min_list_len),
                                                   EMSensitivityParams.f_min.value: np.array([0]*f_min_list_len)},
        MetastableFinalState.EJECTED_FS_MERGER.value: 0
    }

    total_mergers = 0
    for index, evolution in merged_evolutions.iterrows():
        if evolution["end"] in MERGER_ONLY_FINAL_STATES:
            total_mergers += 1
            ev_dbs_state_json = evolution.get("evolution")[-2]
            ev_dbs_state = MetastableState(*ev_dbs_state_json)
            binary_state = BinaryState(**ev_dbs_state.binary_state)
            binary_state.ma = MassObj(**binary_state.ma)
            binary_state.mb = MassObj(**binary_state.mb)

            if evolution["end"] in {MetastableFinalState.IMS_EM_MERGER.value,
                                    MetastableFinalState.IMS_MERGER.value}:
                stats[MetastableFinalState.IMS_MERGER.value] += 1 
                stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.e_min.value] += calc_e_min_sensitivty_ims(binary_state.ma.mass, binary_state.mb.mass, ev_dbs_state.qB)
                stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.f_min.value] += calc_f_min_sensitivty_ims(binary_state.ma.mass, binary_state.mb.mass, ev_dbs_state.qB) 
            if evolution["end"] in {MetastableFinalState.FS_EM_MERGER.value,
                                    MetastableFinalState.FS_MERGER.value}:
                stats[MetastableFinalState.FS_MERGER.value] += 1
                stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.e_min.value] += calc_e_min_sensitivty_fs(binary_state.ma.mass, binary_state.mb.mass, binary_state.eB(), binary_state.EB) 
                stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.f_min.value] += calc_f_min_sensitivty_fs(binary_state.ma.mass, binary_state.mb.mass, binary_state.eB(), binary_state.EB) 
            if evolution["end"] in {MetastableFinalState.EJECTED_FS_EM.value,
                                    MetastableFinalState.EJECTED_FS_MERGER.value}:
                stats[MetastableFinalState.FS_MERGER.value] += 1
                stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.e_min.value] += calc_e_min_sensitivty_fs(binary_state.ma.mass, binary_state.mb.mass, binary_state.eB(), binary_state.EB) 
                stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.f_min.value] += calc_f_min_sensitivty_fs(binary_state.ma.mass, binary_state.mb.mass, binary_state.eB(), binary_state.EB) 

    total_mergers = total_mergers * MULT # a check
    merger_stats = {EMSensitivityParams.e_min.value: 
                    (stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.e_min.value] 
                    + stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.e_min.value] 
                    + stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.e_min.value])/ total_mergers,
                    EMSensitivityParams.f_min.value: 
                    (stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.f_min.value] 
                    + stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.f_min.value] 
                    + stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.f_min.value])/ total_mergers} 
    
    merger_stats_errors = get_merger_stats_errors(merger_stats, total_mergers)
    
    return stats, merger_stats, merger_stats_errors


def get_sensitivity_stats_concrete(file_path, e_min, f_min, verbose=False):
    
    evolutions_df = read_dict_from_file(file_path)
    merged_evolutions = evolutions_df[evolutions_df["end"].isin(MERGER_ONLY_FINAL_STATES)]

    stats = {
        EMSensitivityParams.e_min.value: e_min,
        EMSensitivityParams.f_min.value: f_min,
        MetastableFinalState.IMS_EM_MERGER.value: 0,
        MetastableFinalState.IMS_MERGER.value: 0,
        MetastableFinalState.FS_EM_MERGER.value: 0,
        MetastableFinalState.FS_MERGER.value: 0,
        MetastableFinalState.EJECTED_FS_EM.value: 0,
        MetastableFinalState.EJECTED_FS_MERGER.value: 0
    }

    total_mergers = 0
    for index, evolution in merged_evolutions.iterrows():
        if evolution["end"] in MERGER_ONLY_FINAL_STATES:
            total_mergers += 1
            ev_dbs_state_json = evolution.get("evolution")[-2]
            ev_dbs_state = MetastableState(*ev_dbs_state_json)
            binary_state = BinaryState(**ev_dbs_state.binary_state)
            binary_state.ma = MassObj(**binary_state.ma)
            binary_state.mb = MassObj(**binary_state.mb)

            if evolution["end"] in {MetastableFinalState.IMS_EM_MERGER.value,
                                    MetastableFinalState.IMS_MERGER.value}:
                stats[MetastableFinalState.IMS_MERGER.value] += 1 
                stats[MetastableFinalState.IMS_EM_MERGER.value] += is_merger_eccentric_ims(binary_state.ma.mass, binary_state.mb.mass, ev_dbs_state.qB, e_min, f_min)
            if evolution["end"] in {MetastableFinalState.FS_EM_MERGER.value,
                                    MetastableFinalState.FS_MERGER.value}:
                stats[MetastableFinalState.FS_MERGER.value] += 1
                stats[MetastableFinalState.FS_EM_MERGER.value] += is_merger_eccentric_fs(binary_state.ma.mass, binary_state.mb.mass, binary_state.eB(), binary_state.EB, e_min, f_min) 
            if evolution["end"] in {MetastableFinalState.EJECTED_FS_EM.value,
                                    MetastableFinalState.EJECTED_FS_MERGER.value}:
                stats[MetastableFinalState.EJECTED_FS_MERGER.value] += 1
                stats[MetastableFinalState.EJECTED_FS_EM.value] += is_merger_eccentric_fs(binary_state.ma.mass, binary_state.mb.mass, binary_state.eB(), binary_state.EB, e_min, f_min) 

    total_mergers = total_mergers * MULT # a check
    merger_stats = ((stats[MetastableFinalState.IMS_EM_MERGER.value] 
                    + stats[MetastableFinalState.FS_EM_MERGER.value] 
                    + stats[MetastableFinalState.EJECTED_FS_EM.value])/ total_mergers)

    merger_stats_errors = (merger_stats*(1-merger_stats)/total_mergers)**(1/2)
    
    return stats, merger_stats, merger_stats_errors


def get_sensitivity_stats_old(file_path, verbose=False):
    
    evolutions_list = read_dict_from_file(file_path)
    e_min_list_len = len(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value])
    f_min_list_len = len(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value])

    stats = {
        MetastableFinalState.IMS_EM_MERGER.value: {EMSensitivityParams.e_min.value: np.array([0]*e_min_list_len),
                                                   EMSensitivityParams.f_min.value: np.array([0]*f_min_list_len)},
        MetastableFinalState.IMS_MERGER.value: 0,
        MetastableFinalState.FS_EM_MERGER.value: {EMSensitivityParams.e_min.value: np.array([0]*e_min_list_len),
                                                   EMSensitivityParams.f_min.value: np.array([0]*f_min_list_len)},
        MetastableFinalState.FS_MERGER.value: 0,
        MetastableFinalState.EJECTED_FS_EM.value: {EMSensitivityParams.e_min.value: np.array([0]*e_min_list_len),
                                                   EMSensitivityParams.f_min.value: np.array([0]*f_min_list_len)},
        MetastableFinalState.EJECTED_FS_MERGER.value: 0
    }

    total_mergers = 0
    for index, evolution in enumerate(evolutions_list):
        if evolution["end"] not in {MetastableFinalState.EJECTED_BINARY.value, 
                                    MetastableFinalState.IN_CLUSTER_BINARY.value, 
                                    MetastableFinalState.IONIZED_BINARY.value}:
            total_mergers += 1
            if evolution["end"] in {MetastableFinalState.IMS_EM_MERGER.value,
                                    MetastableFinalState.IMS_MERGER.value}:
                stats[MetastableFinalState.IMS_MERGER.value] += 1 

                ims_state = evolution.get("evolution")[-2].get("IMSes")[-1]
                scramble_data = ims_state.get("ims_sample")
                scramble = get_scramble(scramble_data)
                ims_qB = qB_from_EBeB(ma=scramble.ma, mb=scramble.mb, EB=scramble.EB, eB=scramble.eB())
                stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.e_min.value] += calc_e_min_sensitivty_ims(scramble.ma, scramble.mb, ims_qB)
                stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.f_min.value] += calc_f_min_sensitivty_ims(scramble.ma, scramble.mb, ims_qB) 
            if evolution["end"] in {MetastableFinalState.FS_EM_MERGER.value,
                                    MetastableFinalState.FS_MERGER.value}:
                stats[MetastableFinalState.FS_MERGER.value] += 1

                fs_sample_data = evolution.get("evolution")[-2].get("final_state").get("fs_sample")
                fs_sample = get_scramble(fs_sample_data)
                stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.e_min.value] += calc_e_min_sensitivty_fs(fs_sample.ma, fs_sample.mb, fs_sample.eB(), fs_sample.EB) 
                stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.f_min.value] += calc_f_min_sensitivty_fs(fs_sample.ma, fs_sample.mb, fs_sample.eB(), fs_sample.EB) 
            if evolution["end"] in {MetastableFinalState.EJECTED_FS_EM.value,
                                    MetastableFinalState.EJECTED_FS_MERGER.value}:
                stats[MetastableFinalState.FS_MERGER.value] += 1

                fs_sample_data = evolution.get("evolution")[-2].get("ejected_state").get("fs_sample")
                fs_sample = get_scramble(fs_sample_data)
                stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.e_min.value] += calc_e_min_sensitivty_fs(fs_sample.ma, fs_sample.mb, fs_sample.eB(), fs_sample.EB) 
                stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.f_min.value] += calc_f_min_sensitivty_fs(fs_sample.ma, fs_sample.mb, fs_sample.eB(), fs_sample.EB) 

    total_mergers = total_mergers * MULT # a check
    merger_stats = {EMSensitivityParams.e_min.value: 
                    (stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.e_min.value] 
                    + stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.e_min.value] 
                    + stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.e_min.value])/ total_mergers,
                    EMSensitivityParams.f_min.value: 
                    (stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.f_min.value] 
                    + stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.f_min.value] 
                    + stats[MetastableFinalState.EJECTED_FS_EM.value][EMSensitivityParams.f_min.value])/ total_mergers} 
    
    merger_stats_errors = get_merger_stats_errors(merger_stats, total_mergers)
    
    return stats, merger_stats, merger_stats_errors

def plot_grid(output_files1, output_files2, names, add_detector_lines=True):
    fig = plt.figure(figsize=[9.6, 7.2])
    gs = fig.add_gridspec(2,2, hspace=0, wspace=0)

    (e_plot1, f_plot1), (e_plot2, f_plot2) = gs.subplots(sharex='col', sharey='row')

    colors = plt.cm.tab20b(np.linspace(0, 1, len(output_files1)))
    color_index = 0
    for file_path1, file_path2, name in zip(output_files1, output_files2, names):
        try:
            stats1, merger_stats1, merger_stats_errors1 = get_sensitivity_stats(file_path=file_path1)
            stats2, merger_stats2, merger_stats_errors2 = get_sensitivity_stats(file_path=file_path2)
        except Exception as e:
            print(f"Error in file {file_path1} or {file_path2}: {e}")
            raise e 

        e_plot1.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], 
                        merger_stats1[EMSensitivityParams.e_min.value], 
                        yerr=merger_stats_errors1[EMSensitivityParams.e_min.value], 
                        color=colors[color_index], label=name, capsize=2)
        
        e_plot2.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], 
                        merger_stats2[EMSensitivityParams.e_min.value], 
                        yerr=merger_stats_errors2[EMSensitivityParams.e_min.value], 
                        color=colors[color_index], label=name, capsize=2)

        f_plot1.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], 
                        merger_stats1[EMSensitivityParams.f_min.value], 
                        yerr=merger_stats_errors1[EMSensitivityParams.f_min.value], 
                        color=colors[color_index], label=name, capsize=2)
        
        f_plot2.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], 
                        merger_stats2[EMSensitivityParams.f_min.value], 
                        yerr=merger_stats_errors2[EMSensitivityParams.f_min.value], 
                        color=colors[color_index], label='_nolegend_', capsize=2)
        
        color_index += 1

    colors = plt.cm.Pastel1(np.linspace(0, 1, len(detectors_sensitivity.keys()) + 1))
    color_index = 0
    if add_detector_lines:
        for detector_name in detectors_sensitivity.keys():
            detector = detectors_sensitivity[detector_name]
            f_plot1.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label='_nolegend_') 
            f_plot2.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label=f"{detector_name}")
            color_index += 1

    e_plot1.set_xscale('log')
    e_plot1.set_yscale('log')
    e_plot1.text(x=(10**(-2)), y=9*(10**(-3)), s="$Unequal\,mass\,BBH$")
    
    e_plot2.set_xscale('log')
    e_plot2.set_yscale('log')
    e_plot2.set_xlabel('$e_{min}$')
    e_plot2.text(x=(10**(-2)), y=3*(10**(-2)), s="$Equal\,mass\,BBH$") 

    f_plot1.set_xscale('log')
    f_plot1.set_yscale('log')
    
    f_plot1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    f_plot2.set_xscale('log')
    f_plot2.set_yscale('log')
    f_plot2.set_xlabel('$f_{min}$')
    f_plot2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale('log')
    plt.xscale('log')
    for ax in fig.get_axes():
        ax.label_outer()

    fig.supylabel('$N_{EM}/N_{Merger}$')

def plot_grid_from_files(file_path1, file_path2, names, names_labels=None, add_detector_lines=True, add_linear_fits=False):
    fig = plt.figure(figsize=[9.6, 7.2])
    gs = fig.add_gridspec(2,2, hspace=0, wspace=0)
    
    try:
        all_stats1 = read_dict_from_file(file_path1)
    except Exception as e:
        print(f"Error in file {file_path1}: {e}")
        raise e
    try:
        all_stats2 = read_dict_from_file(file_path2)
    except Exception as e:
        print(f"Error in file {file_path2}: {e}")
        raise e
    
    (e_plot1, f_plot1), (e_plot2, f_plot2) = gs.subplots(sharex='col', sharey='row')

    colors = plt.cm.tab20b(np.linspace(0, 1, len(names)))
    color_index = 0
    linear_fit_label_pos_e_min1 = [(0.017, 0.3), (0.045, 0.065), (0.05, 0.10), (0.05, 0.095)]
    linear_fit_label_pos_e_min2 = [(0.017, 0.035), (0.053, 0.022), (0.07, 0.065), (0.05, 0.11)]
    linear_fit_label_pos_f_min1 = [(4.5, 0.3), (7, 0.055), (12, 0.048), (5, 0.095)]
    linear_fit_label_pos_f_min2 = [(4.5, 0.028), (11, 0.02), (8, 0.065), (5, 0.11)]
    for name in names:
        stats1, merger_stats1, merger_stats_errors1 = all_stats1[name]
        stats2, merger_stats2, merger_stats_errors2 = all_stats2[name]
        name_label = name
        if names_labels and name in names_labels:
            name_label = names_labels[name]

        e_plot1.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], 
                        merger_stats1[EMSensitivityParams.e_min.value], 
                        yerr=merger_stats_errors1[EMSensitivityParams.e_min.value], 
                        color=colors[color_index], label=name_label, capsize=2)
        
        e_plot2.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], 
                        merger_stats2[EMSensitivityParams.e_min.value], 
                        yerr=merger_stats_errors2[EMSensitivityParams.e_min.value], 
                        color=colors[color_index], label=name_label, capsize=2)

        f_plot1.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], 
                        merger_stats1[EMSensitivityParams.f_min.value], 
                        yerr=merger_stats_errors1[EMSensitivityParams.f_min.value], 
                        color=colors[color_index], label=name_label, capsize=2)
        
        f_plot2.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], 
                        merger_stats2[EMSensitivityParams.f_min.value], 
                        yerr=merger_stats_errors2[EMSensitivityParams.f_min.value], 
                        color=colors[color_index], label='_nolegend_', capsize=2)

        
        if add_linear_fits:

            e_min_fit1 = np.polyfit(np.log10(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]), 
                                    np.log10(merger_stats1[EMSensitivityParams.e_min.value]), 1)
            e_min_fit2 = np.polyfit(np.log10(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]), 
                                    np.log10(merger_stats2[EMSensitivityParams.e_min.value]), 1)
            f_min_fit1 = np.polyfit(np.log10(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]), 
                                    np.log10(merger_stats1[EMSensitivityParams.f_min.value]), 1)
            f_min_fit2 = np.polyfit(np.log10(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]), 
                                    np.log10(merger_stats2[EMSensitivityParams.f_min.value]), 1)
            
            extended_e_range = np.linspace(min(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]) * 0.5, 
                             max(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]) * 1.5, 
                             len(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value]))

            extended_f_range = np.linspace(min(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]) * 0.5, 
                             max(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]) * 1.5, 
                             len(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value]))
            e_plot1.plot(extended_e_range, 
                        10**(e_min_fit1[1])*(extended_e_range**(e_min_fit1[0])), 
                        linestyle='--', color=colors[color_index])
            e_plot1.text(linear_fit_label_pos_e_min1[color_index][0], linear_fit_label_pos_e_min1[color_index][1], 
                         f"{e_min_fit1[0]:.2f}",
                        verticalalignment='bottom', horizontalalignment='left', color=colors[color_index])
            e_plot2.plot(extended_e_range, 
                        10**(e_min_fit2[1])*(extended_e_range**(e_min_fit2[0])), 
                        linestyle='--', color=colors[color_index])
            e_plot2.text(linear_fit_label_pos_e_min2[color_index][0], linear_fit_label_pos_e_min2[color_index][1], 
                         f"{e_min_fit2[0]:.2f}",
                        verticalalignment='bottom', horizontalalignment='left', color=colors[color_index])
            f_plot1.plot(extended_f_range, 
                        10**(f_min_fit1[1])*(extended_f_range**(f_min_fit1[0])), 
                        linestyle='--', color=colors[color_index])
            f_plot1.text(linear_fit_label_pos_f_min1[color_index][0], linear_fit_label_pos_f_min1[color_index][1], 
                         f"{f_min_fit1[0]:.2f}",
                        verticalalignment='bottom', horizontalalignment='left', color=colors[color_index])
            f_plot2.plot(extended_f_range, 
                        10**(f_min_fit2[1])*(extended_f_range**(f_min_fit2[0])), 
                        linestyle='--', color=colors[color_index])
            f_plot2.text(linear_fit_label_pos_f_min2[color_index][0], linear_fit_label_pos_f_min2[color_index][1], 
                         f"{f_min_fit2[0]:.2f}",
                        verticalalignment='bottom', horizontalalignment='left', color=colors[color_index])
        
        color_index += 1

    colors = plt.cm.Pastel1(np.linspace(0, 1, len(detectors_sensitivity.keys()) + 1))
    color_index = 0
    if add_detector_lines:
        for detector_name in detectors_sensitivity.keys():
            detector = detectors_sensitivity[detector_name]
            f_plot1.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label='_nolegend_') 
            f_plot2.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label=f"{detector_name}")
            color_index += 1

    e_plot1.set_xscale('log')
    e_plot1.set_yscale('log')
    e_plot1.text(x=(10**(-2)), y=5.5*(10**(-2)), s="$Equal\,mass\,BBH$")
    e_plot1.set_xlim(0.008, 0.4)

    e_plot2.set_xscale('log')
    e_plot2.set_yscale('log')
    e_plot2.set_xlabel('$e_{min}$')
    e_plot2.text(x=(10**(-2)), y=2*(10**(-2)), s="$Unequal\,mass\,BBH$") 
    e_plot2.set_xlim(0.008, 0.4)

    f_plot1.set_xscale('log')
    f_plot1.set_yscale('log')
    f_plot1.set_xlim(3.5, 50)

    f_plot1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    f_plot2.set_xscale('log')
    f_plot2.set_yscale('log')
    f_plot2.set_xlabel('$f_{min}$')
    f_plot2.set_xlim(3.5, 50)
    f_plot2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale('log')
    plt.xscale('log')
    for ax in fig.get_axes():
        ax.label_outer()

    fig.supylabel('$N_{OEM}/N_{Merger}$')




def plot(output_files, names, add_detector_lines=True, save_to_file=None, verbose=False):
    fig, ax = plt.subplots(nrows=2, ncols=1)
    e_plot = ax[0]
    f_plot = ax[1]

    all_stats = {}

    colors = plt.cm.tab20b(np.linspace(0, 1, len(output_files)))
    color_index = 0
    for file_path, name in zip(output_files, names):
        try:
            stats, merger_stats, merger_stats_errors = get_sensitivity_stats(file_path=file_path)
        except Exception as e:
            print(f"Error in file {file_path}: {e}")
            raise e

        merger_ratio_e = merger_stats[EMSensitivityParams.e_min.value]
        merger_ratio_e_errors = merger_stats_errors[EMSensitivityParams.e_min.value] 
        e_plot.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], merger_ratio_e, yerr=merger_ratio_e_errors, 
                        color=colors[color_index], label=name, capsize=2)

        merger_ratio_f = merger_stats[EMSensitivityParams.f_min.value]
        merger_ratio_f_errors = merger_stats_errors[EMSensitivityParams.f_min.value]
        f_plot.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], merger_ratio_f, yerr=merger_ratio_f_errors, 
                        color=colors[color_index], label='_nolegend_', capsize=2)

        all_stats[name] = [stats, merger_stats, merger_stats_errors]
        color_index += 1

    if save_to_file:
        with open(save_to_file, "w+") as outfile: 
            outfile.write(json.dumps(all_stats, cls=EnhancedJSONEncoder))


    colors = plt.cm.Pastel1(np.linspace(0, 1, len(detectors_sensitivity.keys()) + 1))
    color_index = 0
    if add_detector_lines:
        for detector_name in detectors_sensitivity.keys():
            detector = detectors_sensitivity[detector_name]
            f_plot.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label=f"{detector_name}") 
            color_index += 1

    e_plot.set_xscale('log')
    e_plot.set_yscale('log')
    e_plot.set_xlabel('$e_{min}$')
    e_plot.set_ylabel('$N_{EM}/N_{Merger}$')
    e_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    f_plot.set_xscale('log')
    f_plot.set_yscale('log')
    f_plot.set_xlabel('$f_{min}$')
    f_plot.set_ylabel('$N_{EM}/N_{Merger}$')
    f_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale('log')
    plt.xscale('log')
    fig.tight_layout()


def plot_from_file(output_file, add_detector_lines=True, names=None):
    fig, ax = plt.subplots(nrows=2, ncols=1)
    e_plot = ax[0]
    f_plot = ax[1]

    all_stats = read_dict_from_file(output_file)

    colors = plt.cm.tab20b(np.linspace(0, 1, len(all_stats.keys())))
    color_index = 0
    if not names:
        names = all_stats.keys()
    
    for name in names:
        if name not in all_stats:
            print(f"Name {name} not in stats")
            continue
        stats, merger_stats, merger_stats_errors = all_stats[name]
        
        merger_ratio_e = merger_stats[EMSensitivityParams.e_min.value]
        merger_ratio_e_errors = merger_stats_errors[EMSensitivityParams.e_min.value] 
        if any(merger_ratio_e):
            e_plot.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], merger_ratio_e, yerr=merger_ratio_e_errors, 
                            color=colors[color_index], label=name, capsize=2)

        merger_ratio_f = merger_stats[EMSensitivityParams.f_min.value]
        merger_ratio_f_errors = merger_stats_errors[EMSensitivityParams.f_min.value]
        if any(merger_ratio_f):
            f_plot.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], merger_ratio_f, yerr=merger_ratio_f_errors, 
                        color=colors[color_index], label='_nolegend_', capsize=2)

        color_index += 1


    colors = plt.cm.Pastel1(np.linspace(0, 1, len(detectors_sensitivity.keys()) + 1))
    color_index = 0
    if add_detector_lines:
        for detector_name in detectors_sensitivity.keys():
            detector = detectors_sensitivity[detector_name]
            f_plot.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label=f"{detector_name}") 
            color_index += 1

    e_plot.set_xscale('log')
    e_plot.set_yscale('log')
    e_plot.set_xlabel('$e_{min}$')
    e_plot.set_ylabel('$N_{EM}/N_{Merger}$')
    e_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    f_plot.set_xscale('log')
    f_plot.set_yscale('log')
    f_plot.set_xlabel('$f_{min}$')
    f_plot.set_ylabel('$N_{EM}/N_{Merger}$')
    f_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale('log')
    plt.xscale('log')
    fig.tight_layout()


def plot_ims_fs(output_files, names, add_detector_lines=True, ims=True):
   
    fig, ax = plt.subplots(nrows=2, ncols=1)
    e_plot = ax[0]
    f_plot = ax[1]

    colors = plt.cm.tab20b(np.linspace(0, 1, len(output_files)))
    color_index = 0
    for file_path, name in zip(output_files, names):
        stats, _, _ = get_sensitivity_stats(file_path=file_path)

        if ims:
            total_mergers = stats[MetastableFinalState.IMS_MERGER.value] 
            merger_stats = {EMSensitivityParams.e_min.value: 
                        (stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.e_min.value])/total_mergers,
                        EMSensitivityParams.f_min.value: 
                        (stats[MetastableFinalState.IMS_EM_MERGER.value][EMSensitivityParams.f_min.value])/total_mergers}
        else:
            # just FS not FS ejected
            total_mergers = stats[MetastableFinalState.FS_MERGER.value] 
            merger_stats = {EMSensitivityParams.e_min.value: 
                        (stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.e_min.value])/total_mergers,
                        EMSensitivityParams.f_min.value: 
                        (stats[MetastableFinalState.FS_EM_MERGER.value][EMSensitivityParams.f_min.value])/total_mergers}
        
        merger_errors = get_merger_stats_errors(merger_stats, total_mergers)

        merger_ratio_e = merger_stats[EMSensitivityParams.e_min.value]
        merger_ratio_e_errors = merger_errors[EMSensitivityParams.e_min.value] 
        e_plot.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.e_min.value], merger_ratio_e, yerr=merger_ratio_e_errors, 
                        color=colors[color_index], label=name, capsize=2)

        merger_ratio_f = merger_stats[EMSensitivityParams.f_min.value]
        merger_ratio_f_errors = merger_errors[EMSensitivityParams.f_min.value]
        f_plot.errorbar(EM_SENSITIVITY_RANGE[EMSensitivityParams.f_min.value], merger_ratio_f, yerr=merger_ratio_f_errors, 
                        color=colors[color_index], label='_nolegend_', capsize=2)
        
        color_index += 1

    colors = plt.cm.Pastel1(np.linspace(0, 1, len(detectors_sensitivity.keys()) + 1))
    color_index = 0
    if add_detector_lines:
        for detector_name in detectors_sensitivity.keys():
            detector = detectors_sensitivity[detector_name]
            f_plot.axvline(x=detector['freq_range'][0], color=colors[color_index], linestyle='-', label=f"{detector_name}") 
            color_index += 1

    y_label = r'$N_{EM, IMS}/N_{GW, IMS}$' if ims else r'$N_{EM, FS}/N_{GW, FS}$' 

    e_plot.set_xscale('log')
    e_plot.set_yscale('log')
    e_plot.set_xlabel('$e_{min}$')
    e_plot.set_ylabel(y_label)
    e_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    f_plot.set_xscale('log')
    f_plot.set_yscale('log')
    f_plot.set_xlabel('$f_{min}$')
    f_plot.set_ylabel(y_label)
    f_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.tight_layout()
