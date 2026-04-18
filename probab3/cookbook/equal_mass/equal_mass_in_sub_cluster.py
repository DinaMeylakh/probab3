"""A cookbook script to sample equal mass triple systems in a 2-component star cluster. 

The fist componet are the stars and the second component is the black holes sub-cluster.
"""
from probab3.commands.common.formulas.general import *
from probab3.commands.sample_dbs import dbs_sampler_equal_mass
from probab3.commands.plot import plotter
from probab3.commands.pre_calculation.cluster_params.cluster_params_relations import select_rh_for_Mc
from probab3.commands.post_analysis import sensitivity_dependence, manual_exchange_stats
from probab3.commands.common.general_code import merge_jsonl_files
from probab3.commands.common.formulas.star_cluster import TwoComponentCluster
import click
import math
import glob
import os

# Parameters
# number of samples
N = 4000
# The folder where the samples will be saved
OUTPUT_FOLDER = "sampling_logs/equal_mass_bbh/" 
# The suffix of the sample files
file_suffix = "hopefully_final_try2"

defualt_rc_index = 2
default_Mc_index = 4

Mcs = np.logspace(5, (8), 10)
rhs = select_rh_for_Mc(Mcs) 
rcs = np.logspace(math.log((1/20), 10), math.log(0.5, 10), 5)

def sample_Mc_rh_rc(Mc_value: int, Mi: int, rh_value: float, rhi: int, rc_value: float, rci: int, 
                    output_folder: str, option_suffix:str):
    """Sample the triple system N times and save the sample to a file.
    
    Args:
        Mc_value (int): The mass of the cluster in solar masses.
        Mi (int): The index of the mass of the cluster in the Mcs array.
        rh_value (float): The half-mass radius of the cluster in parsecs.
        rhi (int): The index of the half-mass radius of the cluster in the rhs array.
        rc_value (float): The core radius of the cluster in units of rh.
        rci (int): The index of the core radius of the cluster in the rcs array.
        output_folder (str): The folder where the samples will be saved.
    """
    sample_path = f"{output_folder}/sample_Mc{Mi}_rh{rhi}_rc{rci}_n{N}_{option_suffix}_{file_suffix}.jsonl"
    dbs_sampler_equal_mass.ms_to_BH_to_end_sample(alpha=2.5, 
                                                 Mc=Mc_value*MSun, 
                                                 rh=rh_value*parsec, 
                                                 rc=rc_value*rh_value*parsec,
                                                 size=N,
                                                 output_file_path=sample_path)

def run(option="A01"):
    """Run the cookbook script. 
    
    Initialize most probable cluster parameters and sample the triple system N times for each parameter set.
    Different options run different parameter sets. This is used to run the script in parallel with different processes.

    Args:
        option (str, optional): The option to run. Defaults to "1".
    """
    option_num = int(option[-2:])
    suffix_letter = option[0]
    if option_num < len(Mcs):
        click.secho(f"Got option num for Mcs: {option_num} with suffix {suffix_letter}", fg="green")
        click.secho(f"Starting sampling Mc={Mcs[option_num]} ({option_num}/{len(Mcs)-1}), "
                    f"rh={rhs[option_num]} ({option_num}/{len(rhs)-1})," 
                    f"rc={rcs[defualt_rc_index]} ({defualt_rc_index}/{len(rcs)-1})", fg="green")

        sample_Mc_rh_rc(Mc_value=Mcs[option_num], Mi=option_num, rh_value=rhs[option_num], 
                        rhi=option_num, rc_value=rcs[defualt_rc_index], rci=defualt_rc_index, 
                        output_folder=OUTPUT_FOLDER, option_suffix=suffix_letter)
    
    elif option_num < len(Mcs) + len(rcs):
        rcs_option_num = option_num - len(Mcs)
        click.secho(f"Got option num {option_num} for rcs: {rcs_option_num}", fg="green")
        click.secho(f"Starting sampling Mc={Mcs[default_Mc_index]} ({default_Mc_index}/{len(Mcs)-1}), "
                    f"rh={rhs[default_Mc_index]} ({default_Mc_index}/{len(rhs)-1})," 
                    f"rc={rcs[rcs_option_num]} ({rcs_option_num}/{len(rcs)-1})", fg="green")

        sample_Mc_rh_rc(Mc_value=Mcs[default_Mc_index], Mi=default_Mc_index, rh_value=rhs[default_Mc_index], 
                        rhi=default_Mc_index, rc_value=rcs[rcs_option_num], rci=rcs_option_num, 
                        output_folder=OUTPUT_FOLDER, option_suffix=suffix_letter)
    

def options():
    """Display the options of the cookbook script."""
    letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
    nums = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"]
    return [f"{letter}{num}" for letter in letters for num in nums] 

def plot2(folder_path: str = OUTPUT_FOLDER, include_binary: bool = True, merger_only: bool = False):
    """Plot the results."""
    
    if os.path.exists(f'{folder_path}Mcs_equal_bbh_stats_data_{file_suffix}.json'):
        plotter.plot_parameters_stats_from_file(file_path=f'{folder_path}Mcs_equal_bbh_stats_data_{file_suffix}.json',
                                                merger_only=merger_only)
    
    else:
        sample_paths_Mcs = glob.glob(f"{folder_path}*_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl", recursive=False)
        plotter.plot_few_cluster_files_stats(sample_paths_Mcs, 'Mc', include_binary=include_binary, merger_only=merger_only,
                                             save_plot_data_path=f'{folder_path}Mcs_equal_bbh_stats_data_{file_suffix}.json')

    if os.path.exists(f'{folder_path}cs_equal_bbh_stats_data_{file_suffix}.json'):
        plotter.plot_parameters_stats_from_file(f'{folder_path}cs_equal_bbh_stats_data_{file_suffix}.json',
                                                merger_only=merger_only)
    else:
        sample_paths_rcs = glob.glob(f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc*_n{N}_{file_suffix}.jsonl", recursive=False)

        plotter.plot_few_cluster_files_stats(sample_paths_rcs, 'c', include_binary=include_binary, merger_only=merger_only,
                                            save_plot_data_path=f'{folder_path}cs_equal_bbh_stats_data_{file_suffix}.json')

    plot_sensetivity_dep(folder_path)

def plot_from_saved(folder_path: str = OUTPUT_FOLDER):
    plotter.plot_parameters_stats_from_file(f'{folder_path}Mcs_equal_bbh_stats_data_{file_suffix}.json')
    plotter.plot_parameters_stats_from_file(f'{folder_path}cs_equal_bbh_stats_data_{file_suffix}.json')


def plot_sensetivity_dep(folder_path: str = OUTPUT_FOLDER):
    #sample_paths_Mcs = glob.glob(f"{folder_path}*_rc{defualt_rc_index}_n{N}_{file_suffix}.json", recursive=False)
    sample_paths_Mcs = [
        f"{folder_path}sample_Mc0_rh0_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc1_rh1_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc2_rh2_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc3_rh3_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc4_rh4_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc5_rh5_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc6_rh6_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc7_rh7_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc8_rh8_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc9_rh9_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl"
    ]

    sensitivity_dependence.plot(output_files=sample_paths_Mcs, 
                                names=["$M_{tot}=10^5 M_{\\odot}$", 
                                       "$M_{tot}=10^{5.33} M_{\\odot}$",
                                       "$M_{tot}=10^{5.67} M_{\\odot}$",
                                       "$M_{tot}=10^6 M_{\\odot}$",
                                       "$M_{tot}=10^{6.33} M_{\\odot}$",
                                       "$M_{tot}=10^{6.67} M_{\\odot}$", 
                                       "$M_{tot}=10^7 M_{\\odot}$",
                                       "$M_{tot}=10^{7.33} M_{\\odot}$",
                                       "$M_{tot}=10^{7.67} M_{\\odot}$",
                                       "$M_{tot}=10^8 M_{\\odot}$"],
                                save_to_file=f"{folder_path}Mcs_sensitivity_dependence_{file_suffix}.json")
    
    #sample_paths_rcs = glob.glob(f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc*_n{N}_{file_suffix}.json", recursive=False)
    sample_paths_rcs = [
        f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc0_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc1_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc2_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc3_n{N}_{file_suffix}.jsonl",
        f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc4_n{N}_{file_suffix}.jsonl"
    ]

    sensitivity_dependence.plot(output_files=sample_paths_rcs, 
                            names=["$r_c=0.05 r_h$", "$r_c=0.089 r_h$", 
                                   "$r_c=0.158 r_h$", "$r_c=0.281 r_h$", 
                                   "$r_c=0.5 r_h$"],
                            save_to_file=f"{folder_path}cs_sensitivity_dependence_{file_suffix}.json")


def plot_sensetivity_dep_from_saved(folder_path: str = OUTPUT_FOLDER):
    sensitivity_dependence.plot_from_file(f"{folder_path}Mcs_sensitivity_dependence_{file_suffix}.json")
    sensitivity_dependence.plot_from_file(f"{folder_path}cs_sensitivity_dependence_{file_suffix}.json")


def plot(folder_path: str = OUTPUT_FOLDER):
    sample_paths_Mcs = glob.glob(f"{folder_path}*_rc{defualt_rc_index}_n{N}_{file_suffix}.jsonl", recursive=False)
    manual_exchange_stats.save_manual_exchange_stats(file_paths=sample_paths_Mcs, save_path=f"{folder_path}Mcs_manual_exchange_stats_{file_suffix}.json", verbose=True, cluster_class=TwoComponentCluster)

    sample_paths_rcs = glob.glob(f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc*_n{N}_{file_suffix}.jsonl", recursive=False)
    manual_exchange_stats.save_manual_exchange_stats(file_paths=sample_paths_rcs, save_path=f"{folder_path}cs_manual_exchange_stats_{file_suffix}.json", verbose=True, cluster_class=TwoComponentCluster, run_param='c')


def merge_files(folder_path: str = OUTPUT_FOLDER):
    """Merge the sample files into one file."""
    for Mi, Mc_value in enumerate(Mcs):
        for rhi, rh_value in enumerate(rhs):
            for rci, rc_value in enumerate(rcs):
                input_pattern = f"{folder_path}/sample_Mc{Mi}_rh{rhi}_rc{rci}_n{N}_*_{file_suffix}.jsonl"
                output_path = f"{folder_path}/sample_Mc{Mi}_rh{rhi}_rc{rci}_n{N}_{file_suffix}.jsonl"
                click.echo(f"Merging {input_pattern} into {output_path}")
                merge_jsonl_files(input_pattern, output_path)
                click.echo(f"Merge complete.")