from probab3.commands.common.formulas.general import *
from probab3.commands.common.general_code import *
from probab3.commands.post_analysis.get_sample_stats import * 

from matplotlib import pyplot as plt
import click

stats_colors = {}

def plot_evolution_from_data(data_dict: dict, pre_evolution_dict: dict = {}, verbose=False):

    if pre_evolution_dict: 
        pre_ev_aBs = np.array(pre_evolution_dict.pop('all_aBs'))/AU
        pre_ev_qBs = np.array(pre_evolution_dict.pop('all_qBs'))/AU
        pre_aBs_tidal_limit = np.array(pre_evolution_dict.pop('all_aBs_tidal_limit'))/AU
        pre_ev_Ns = np.array(list(range(-len(pre_ev_aBs)+1, 1)))
        if verbose:
            click.echo(f"pre_ev_Ns: {pre_ev_Ns}, pre_ev_aBs: {pre_ev_aBs}, pre_ev_qBs: {pre_ev_qBs}")

    all_aBs = np.array(data_dict.pop('all_aBs'))/AU
    all_qBs = np.array(data_dict.pop('all_qBs'))/AU

    all_qBs_for_merger = [] 
    if None not in data_dict.get('all_qBs_for_merger'):
        all_qBs_for_merger = np.array(data_dict.pop('all_qBs_for_merger'))/AU

    all_qBs_for_EM = []
    if None not in data_dict.get('all_qBs_for_EM'):
        all_qBs_for_EM = np.array(data_dict.pop('all_qBs_for_EM'))/AU

    fs_Ns = np.array(data_dict.pop('fs_Ns', []))
    ims_Ns = np.array(data_dict.pop('ims_Ns', []))
    ejected_Ns = np.array(data_dict.pop('ejected_Ns', []))

    all_aBs_tidal_limit = []
    if None not in data_dict.get('all_aBs_tidal_limit'):
        all_aBs_tidal_limit = np.array(data_dict.pop('all_aBs_tidal_limit'))/AU
    N = len(all_aBs)
    
    fs_aBs = all_aBs[fs_Ns] if len(fs_Ns) > 0 else []
    fs_qBs = all_qBs[fs_Ns] if len(fs_Ns) > 0 else [] 

    ims_aBs = all_aBs[ims_Ns] if len(ims_Ns) > 0 else [] 
    ims_qBs = all_qBs[ims_Ns] if len(ims_Ns) > 0 else []  
    if verbose:
        click.echo(f"fs_Ns: {fs_Ns}, ims_Ns: {ims_Ns}, ejected_Ns: {ejected_Ns}, ims_aBs: {ims_aBs}, ims_qBs: {ims_qBs}, all_aBs: {all_aBs}, all_qBs: {all_qBs}")
        click.echo(f"all_qBs_for_merger: {all_qBs_for_merger}, all_qBs_for_EM: {all_qBs_for_EM}")
    Ns = list(range(1, N + 1)) 
    plt.yscale('log')

    if pre_evolution_dict:
        united_Ns = np.concatenate((pre_ev_Ns,Ns)) 
        united_aBs = np.concatenate((pre_ev_aBs, all_aBs))
        united_qBs = np.concatenate((pre_ev_qBs, all_qBs))
        united_tidal_aBs = np.concatenate((pre_aBs_tidal_limit, all_aBs_tidal_limit))
        plt.plot(united_Ns, united_aBs)
        plt.plot(united_Ns, united_qBs)
        plt.plot(pre_ev_Ns, pre_ev_aBs, "s", label=r'PEV $a_B$', marker="h", markersize=4, markerfacecolor="blue", markeredgecolor="purple")
        plt.plot(pre_ev_Ns, pre_ev_qBs, "s", label=r'PEV $q_B$', marker='h', markersize=4, markerfacecolor="orange", markeredgecolor="red")
        plt.plot(united_Ns, united_tidal_aBs, label=r'$a_B$ Tidal Limit', linestyle="dashed", color="gray")

    else:
        plt.plot(Ns, all_aBs)
        plt.plot(Ns, all_qBs)

    if len(ims_Ns) > 0:
        plt.plot(ims_Ns + 1, ims_aBs, "s", label=r'IMS $a_B$', marker=".", markerfacecolor="blue")
        plt.plot(ims_Ns + 1, ims_qBs, "s", label=r'IMS $q_B$', marker='.', markerfacecolor="orange", markeredgecolor="none")
    if len(fs_Ns) > 0:
        plt.plot(fs_Ns + 1, fs_aBs, "s", label=r'FS $a_B$', marker="o", markersize=4, markerfacecolor="blue", markeredgecolor="purple")
        plt.plot(fs_Ns + 1, fs_qBs, "s", label=r'FS $q_B$', marker='o', markersize=4, markerfacecolor="orange", markeredgecolor="red")
    if len(ejected_Ns) > 0:
        plt.plot([Ns[-1]], [all_aBs[ejected_Ns]], "s", label='ejected aB', marker='s', markersize=4, markerfacecolor="purple", markeredgecolor="purple") 
    
    if ejected_Ns:
        plt.plot([Ns[-1]], [all_qBs[ejected_Ns]], "s", label=r'Ejected $q_B$', marker='s', markersize=4, markerfacecolor="orange", markeredgecolor="red") 

    if len(all_qBs_for_merger) == len(Ns):
        plt.plot(Ns, all_qBs_for_merger, "c", label=r'$q_B$ Merger', marker='.', markersize=4, markerfacecolor="none", color="turquoise")
   
    if len(all_qBs_for_EM) == len(Ns): 
        plt.plot(Ns, all_qBs_for_EM, label=r'$q_B$ OEM', markerfacecolor="none", color="hotpink")

    if not pre_evolution_dict and len(all_aBs_tidal_limit) == len(Ns):
        plt.plot(Ns, all_aBs_tidal_limit, label=r'$a_B$ Tidal Limit', linestyle="dashed", color="gray")

    plt.legend(loc='lower left', bbox_to_anchor=(-0.15, 1), ncol=6, columnspacing=0.8)
    plt.xlabel('N')
    plt.ylabel('separation [AU]')
    #return plt
    plt.show()

def plot_evolution_from_data_file(file_path, evolution_num=0, print_details=False, save_plot_data_path=None, pre_evolution=True, verbose=False):
    if print_details or save_plot_data_path:
        run_details = get_run_details(file_path=file_path, evolution_num=evolution_num)
        if print_details:
            click.echo(run_details)
    evolution = get_evolution_from_file(file_path=file_path, evolution_num=evolution_num)
    pre_evolution_dict = get_pre_evolution_stats(evolution.pop('pre_evolution')) if pre_evolution else {}

    evolution_stats = get_evolution_stats(evolution.pop('evolution'))
    if save_plot_data_path:
        evolution_stats['run_details'] = run_details
        save_numpy_dict_to_json_file(save_plot_data_path, evolution_stats)

    plot_evolution_from_data(evolution_stats, pre_evolution_dict, verbose=verbose)

def plot_evolution_from_stats_file(file_path, print_details=False, verbose=False):
    evolution_stats = read_dict_from_json_file(file_path)
    run_details = evolution_stats.pop('run_details')
    if print_details:
        click.echo(run_details) 
    plot_evolution_from_data(evolution_stats, verbose=verbose)


def plot_few_files_stats(file_paths, param_for_plot='x', normali_name="units", include_binary=True, save_plot_data_path=None):
    param_stats, param_stats_errors, param_values, ev_nums = get_few_file_stats(file_paths=file_paths, 
                                                                                verbose=False,
                                                                                run_param=param_for_plot)

    plot_stats_kwargs = {
        'param_values': param_values,
        'params_stats': param_stats,
        'param_name': f"{param_for_plot} [{normali_name}]",
        'include_binary': include_binary,
        'y_err': param_stats_errors,
        'evolutions_per_file': ev_nums
    }

    if save_plot_data_path and include_binary:
        click.echo(f"Saving stats: {plot_stats_kwargs}\n to file {save_plot_data_path}")
        save_numpy_dict_to_json_file(save_plot_data_path, plot_stats_kwargs)

    plot_stats_kwargs.pop('evolutions_per_file') 

    return plot_cluster_params_stats(**plot_stats_kwargs) 


def plot_few_cluster_files_stats(file_paths, cluster_param_for_plot='rc', include_binary=True, 
                                 merger_only=False, check_in_cluster=False, save_plot_data_path=None):

    allowed_final_states = ALL_FINAL_STATES 
    if merger_only:
        allowed_final_states = MERGER_ONLY_FINAL_STATES
    cluster_param_stats, cluster_param_stats_errors, _, ev_nums = get_few_file_stats(file_paths=file_paths, 
                                                                         verbose=False, 
                                                                         allowed_final_states=allowed_final_states,
                                                                         run_param=cluster_param_for_plot)
    normali_name = {
        'c': '',
        'rc': r'$[pc]$',
        'rh': r'$[pc]$',
        'Mc': r'$[M_{\odot}]$',
        'h': ''
    }
    param_label = {
        'c': r'$c_{\star}$',
        'rc': r'$r_c$',
        'rh': r'$r_h$',
        'Mc': r'$M_{\rm tot}$',
        'h': r'$h$'
    }

    click.secho(f"Getting {cluster_param_for_plot} values for plotting.", fg="green")
    cluster_param_values, dispersion_vels = get_cluster_param_from_files(file_paths, cluster_param_for_plot)
    
    plot_stats_kwargs = {
        'param_values': cluster_param_values,
        'params_stats': cluster_param_stats,
        'param_name': f"{param_label[cluster_param_for_plot]} {normali_name[cluster_param_for_plot]}",
        'include_binary': include_binary,
        'y_err': cluster_param_stats_errors,
        'top_x_values': dispersion_vels,
        'top_x_name': r"$\sigma_{c}$ [km/s]",
        'evolutions_per_file': ev_nums
    }

    if save_plot_data_path and include_binary and not merger_only:
        click.secho(f"Saving stats {plot_stats_kwargs} \n to file {save_plot_data_path}", fg="green")
        save_numpy_dict_to_json_file(save_plot_data_path, plot_stats_kwargs)

    plot_stats_kwargs.pop('evolutions_per_file') 

    click.secho(f"Plotting.", fg="green")
    plot_cluster_params_stats(**plot_stats_kwargs) 


def plot_cluster_params_in_ax(ax, param_values, params_stats, include_binary=True, y_err=None, in_cluster_mergers_stats=[], in_cluster_merger_param_values=[]):
    minimal_value = 1
    colors = plt.cm.Set3(np.linspace(0, 1, len(params_stats.keys())))
    for color_index, stat in enumerate(sorted(params_stats.keys())):
        if not include_binary and stat in [MetastableFinalState.IN_CLUSTER_BINARY.value, 
                                           MetastableFinalState.EJECTED_BINARY.value, 
                                           MetastableFinalState.IONIZED_BINARY.value]:
            continue

        if not any(params_stats[stat]):
            continue

        color = colors[color_index]
        if stat in FS_COLORS.keys():
            color = FS_COLORS[stat] 
        stat_label = stat
        if stat in FS_LABELS.keys():
            stat_label = FS_LABELS[stat]
        sorted_param_values = sorted(np.array(param_values))
        sorted_param_stats = np.array(params_stats[stat])[np.argsort(np.array(param_values))]
        if (y_err is not None) and (len(y_err.get(stat, [])) > 0):
            sorted_stat_errors = np.array(y_err.get(stat))[np.argsort(np.array(param_values))]
            ax.errorbar(sorted_param_values, sorted_param_stats, yerr=sorted_stat_errors, fmt=".-", label=stat_label, capsize=2, color=color)
        else:
            ax.plot(sorted_param_values, sorted_param_stats, ".-", label=stat_label, color=color)
    
        if in_cluster_mergers_stats and stat in in_cluster_mergers_stats.keys():
            ax.plot(in_cluster_merger_param_values, in_cluster_mergers_stats[stat], marker='o', color='black', markerfacecolor='none', linestyle='none')
        
        y_values = np.array(params_stats[stat]) 
        positive_ys = y_values[y_values > 0]
        if positive_ys.any() and min(positive_ys) < minimal_value:
            minimal_value = min(positive_ys)

    stat = MetastableFinalState.EJECTED_FS_EM.value 
    sorted_param_values = sorted(np.array(param_values))
    sorted_param_stats = np.array(params_stats[stat])[np.argsort(np.array(param_values))]
    if (y_err is not None) and (len(y_err.get(stat, [])) > 0):
        sorted_stat_errors = np.array(y_err.get(stat))[np.argsort(np.array(param_values))]
        ax.errorbar(sorted_param_values, sorted_param_stats, yerr=sorted_stat_errors, fmt=".-", label=FS_LABELS[stat], capsize=2, color=FS_COLORS[stat])

    ax.set_ylim(ymin=minimal_value/10, ymax=2)
    ax.set_xlim(xmin=min(sorted_param_values)/1.1, xmax=max(sorted_param_values)*1.1)


def plot_cluster_params_stats(param_values, params_stats, param_name, include_binary=True, y_err=None, 
                              top_x_values=[], top_x_name='', in_cluster_mergers_stats=[], 
                              in_cluster_merger_param_values=[]):


    fig, ax1 = plt.subplots(1,1)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    if top_x_values != []:
        ax12 = ax1.twiny()
        ax12.set_xscale('log')
        ax12.minorticks_off()

    plot_cluster_params_in_ax(ax1, param_values, params_stats, include_binary, y_err, in_cluster_mergers_stats, in_cluster_merger_param_values)

    if top_x_values != []:
        ax12.set_xlim(ax1.get_xlim())
        ax12.set_xticks(param_values)
        ax12.set_xticklabels(["{:.0f}".format(top_x_value) for top_x_value in top_x_values])

        if top_x_name != '':
            ax12.set_xlabel(top_x_name)

    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax1.set_xlabel(param_name)
    ax1.set_ylabel('Sampled Ratio')
    fig.tight_layout()
 
    return fig


def plot_two_panels_same_x_values_from_file_datas(top_data, bottom_data, x_label, y_label, top_label='', bottom_label=''):

    fig = plt.figure(figsize=(6.4, 9.6))
    gridspec = fig.add_gridspec(2, hspace=0)
    axs = gridspec.subplots(sharex=True)
    
    ax1 = axs[0]
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    # top x ticks for top panel 
    ax12 = ax1.twiny()
    ax12.set_xscale('log')
    ax12.minorticks_off()

    plot_cluster_params_in_ax(ax1, top_data.get('param_values'), top_data.get('params_stats'), 
                              top_data.get('include_binary'), top_data.get('y_err'), 
                              top_data.get('in_cluster_mergers_stats'), top_data.get('in_cluster_merger_param_values'))

    if top_data.get('top_x_values'):
        ax12.set_xlim(ax1.get_xlim())
        ax12.set_xticks(top_data.get('param_values'))
        ax12.set_xticklabels(["{:.0f}".format(top_x_value) for top_x_value in top_data.get('top_x_values')])

        if top_data.get('top_x_name'):
            ax12.set_xlabel(top_data.get('top_x_name'))

    

    ax2 = axs[1]

    plot_cluster_params_in_ax(ax2, bottom_data.get('param_values'), bottom_data.get('params_stats'), 
                              bottom_data.get('include_binary'), bottom_data.get('y_err'), 
                              bottom_data.get('in_cluster_mergers_stats'), bottom_data.get('in_cluster_merger_param_values')) 

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(x_label)
    fig.supylabel(y_label)
    ax1.set_yscale('log')

    ax2.legend(loc='lower left', bbox_to_anchor=(-0.15, -0.4), ncol=4, columnspacing=0.7)

    plt.show()

def plot_parameters_stats_from_file(file_path, merger_only=False):
    
    plot_data = get_plot_stats_from_file(file_path, merger_only) 

    plot_cluster_params_stats(**plot_data) 

def get_plot_stats_from_file(file_path, merger_only=False):

    allowed_final_states = ALL_FINAL_STATES 
    if merger_only:
        allowed_final_states = MERGER_ONLY_FINAL_STATES

    with open(file_path, 'r') as file_obj:
        file_data = json.loads(file_obj.read())
    
    evolution_nums = file_data.pop('evolutions_per_file')

    plot_data = copy.deepcopy(file_data)
    if merger_only:
        params_stats = plot_data.pop('params_stats')
        plot_data.pop('y_err') 
        new_stats, new_stats_errors = renormalize_stats_data(params_stats=params_stats, allowed_final_states=allowed_final_states, 
                                          evolution_nums=evolution_nums) 
        plot_data['params_stats'] = new_stats
        plot_data['y_err'] = new_stats_errors 
        plot_data['include_binary'] = False

    return plot_data 


def plot_two_panels_same_x_values_from_files(top_file_path, bottom_file_path, x_label, y_label, top_label='', bottom_label='', merger_only=False):

    top_plot_data = get_plot_stats_from_file(top_file_path, merger_only) 
    bottom_plot_data = get_plot_stats_from_file(bottom_file_path, merger_only)

    return plot_two_panels_same_x_values_from_file_datas(top_plot_data, bottom_plot_data, x_label, y_label, top_label, bottom_label)