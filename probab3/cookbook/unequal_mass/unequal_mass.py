"""A cookbook script to sample unequal mass triple systems in a 2-component star cluster. 

The fist componet are the stars and the second component is the black holes sub-cluster.

All masses are sampled from a present day mass function (PDMF).
"""
from probab3.commands.common.formulas.general import *
from probab3.commands.sample_dbs import dbs_sampler_unequal_mass
from probab3.commands.plot import plotter
from probab3.commands.pre_calculation.cluster_params.cluster_params_relations import *
from probab3.commands.post_analysis import sensitivity_dependence
import click
import math
import glob

# Parameters
# number of samples
N = 2500 
# The folder where the samples will be saved
new_NC_output_folder = "sampling_logs/unequal_masses/" 
# The suffix of the sample files
file_suffix = "try16_mb_new"

defualt_rc_index = 2
default_Mc_index = 4


def sample_Mc_rh_rc(Mc_value: int, Mi: int, rh_value: float, rhi: int, rc_value: float, rci: int, 
                    output_folder: str):
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
    sample_path = f"{output_folder}/sample_Mc{Mi}_rh{rhi}_rc{rci}_n{N}_{file_suffix}.json"
    dbs_sampler_unequal_mass.unequal_mass_sample(alpha=2.5, 
                                                Mc=Mc_value*MSun, 
                                                rh=rh_value*parsec, 
                                                rc=rc_value*rh_value*parsec,
                                                size=N,
                                                output_file_path=sample_path)

def run_cluster(Mcs, rhs, rcs, output_folder):
    """Iterate over different cluster parameters and sample the triple system N times in each.

    Save the samples to files in the output folder.
    
    Args:
        Mcs (list): The mass of the cluster in solar masses.
        rhs (list): The half-mass radius of the cluster in parsecs.
        rcs (list): The core radius of the cluster in units of rh.
        output_folder (str): The folder where the samples will be saved.
    """
    for Mi, Mc_value in enumerate(Mcs):
        click.secho(f"Starting sampling Mc={Mc_value} ({Mi}/{len(Mcs)}), "
                    f"rh={rhs[0]} ({0}/{len(rhs)})," 
                    f"rc={rcs[0]} ({0}/{len(rcs)})", fg="green")
        sample_Mc_rh_rc(Mc_value, Mi, rh_value=rhs[0], rhi=0, rc_value=rcs[0], rci=0, output_folder=output_folder)

    for rhi, rh_value in enumerate(rhs[1:]):
        rhi = rhi + 1
        click.secho(f"Starting sampling Mc={Mcs[0]} ({0}/{len(Mcs)}), "
                    f"rh={rh_value} ({rhi}/{len(rhs)})," 
                    f"rc={rcs[0]} ({0}/{len(rcs)})", fg="green")
        sample_Mc_rh_rc(Mc_value=Mcs[0], Mi=0, rh_value=rh_value, rhi=rhi, rc_value=rcs[0], rci=0, output_folder=output_folder)

    for rci, rc_value in enumerate(rcs[1:]):
        rci = rci +1
        click.secho(f"Starting sampling Mc={Mcs[0]} ({0}/{len(Mcs)}), "
                    f"rh={rhs[0]} ({0}/{len(rhs)})," 
                    f"rc={rc_value} ({rci}/{len(rcs)})", fg="green")
        sample_Mc_rh_rc(Mc_value=Mcs[0], Mi=0, rh_value=rhs[0], rhi=0, rc_value=rc_value, rci=rci, output_folder=output_folder)


def run(option="1"):
    """Run the cookbook script. 
    
    Initialize most probable cluster parameters and sample the triple system N times for each parameter set.
    Different options run different parameter sets. This is used to run the script in parallel with different processes.

    Args:
        option (str, optional): The option to run. Defaults to "1".
    """
    Mcs = np.logspace(5, (8), 10)
    rhs = select_rh_for_Mc(Mcs) 
    rcs = np.logspace(math.log((1/20), 10), math.log(0.5, 10), 5)
 

    option_num = int(option)
    if option_num < len(Mcs):
        click.secho(f"Got option num for Mcs: {option_num}", fg="green")
        click.secho(f"Starting sampling Mc={Mcs[option_num]} ({option_num}/{len(Mcs)-1}), "
                    f"rh={rhs[option_num]} ({option_num}/{len(rhs)-1})," 
                    f"rc={rcs[defualt_rc_index]} ({defualt_rc_index}/{len(rcs)-1})", fg="green")

        sample_Mc_rh_rc(Mc_value=Mcs[option_num], Mi=option_num, rh_value=rhs[option_num], 
                        rhi=option_num, rc_value=rcs[defualt_rc_index], rci=defualt_rc_index, 
                        output_folder=new_NC_output_folder)
    
    elif option_num < len(Mcs) + len(rcs):
        rcs_option_num = option_num - len(Mcs)
        click.secho(f"Got option num {option_num} for rcs: {rcs_option_num}", fg="green")
        click.secho(f"Starting sampling Mc={Mcs[default_Mc_index]} ({default_Mc_index}/{len(Mcs)-1}), "
                    f"rh={rhs[default_Mc_index]} ({default_Mc_index}/{len(rhs)-1})," 
                    f"rc={rcs[rcs_option_num]} ({rcs_option_num}/{len(rcs)-1})", fg="green")

        sample_Mc_rh_rc(Mc_value=Mcs[default_Mc_index], Mi=default_Mc_index, rh_value=rhs[default_Mc_index], 
                        rhi=default_Mc_index, rc_value=rcs[rcs_option_num], rci=rcs_option_num, 
                        output_folder=new_NC_output_folder)
    

def options():
    """Display the options of the cookbook script."""
    return [f"{i}" for i in range(20)] 

def plot(folder_path: str = new_NC_output_folder, include_binary: bool = True, merger_only: bool = False):
    """Plot the results."""

    sample_paths_Mcs = glob.glob(f"{folder_path}*_rc{defualt_rc_index}_n{N}_{file_suffix}.json", recursive=False)

    plotter.plot_few_cluster_files_stats(sample_paths_Mcs, 'Mc', include_binary=include_binary, merger_only=merger_only)

    sample_paths_rcs = glob.glob(f"{folder_path}sample_Mc{default_Mc_index}_rh{default_Mc_index}_rc*_n{N}_{file_suffix}.json", recursive=False)

    plotter.plot_few_cluster_files_stats(sample_paths_rcs, 'c', include_binary=include_binary, merger_only=merger_only)

