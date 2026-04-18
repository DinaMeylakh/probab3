#!/usr/bin/python

import click
import pkg_resources
from probab3.commands.common.formulas.general import *
from probab3.commands.common.constants import *
from probab3.commands.common.general_code import *
from probab3.commands.plot import plotter
from probab3.commands.post_analysis import get_sample_stats, manual_exchange_stats
from probab3.cli_sample_dbs import sample_dbs
import time
import pprint
import logging
import click_log
from colorama import init as colorama_init
from colorama import Fore
from colorama import Style
import sys
import importlib.util
import importlib.metadata

logging.basicConfig(level=logging.INFO, filename=LOG_FILE_NAME, filemode='a', format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(PROJECT_NAME)
colorama_init()

version_num = importlib.metadata.version(PROJECT_NAME)

@click.group()
def cli():
    click.secho(f'You are using {PROJECT_NAME} version {version_num}', fg="yellow")

@cli.command(name='version', help='show version and exit.')
def version():
    pass

@cli.command(name='plot', help='plot a dynamical binary sequence evolution from a file.')
@click.option('--input-path', default='./current_sample.json', help='input file path of the sample.')
@click.option('--show-run-details', is_flag=True, default=False, help='Print run details along with the plot.')
@click.option('--N', default=0, help='The run number to plot.')
@click.option('--save-stats-path', default=None, help='Save the stats to a file.')
@click.option('--from-stats', is_flag=True, default=False, help='Plot from stats file instead of data file.')
@click.option('--pre-evolution', is_flag=True, default=False, help='Include the pre-evolution in the evolution plot. Only available with the --N option.')
@click.option('--verbose', is_flag=True, default=False, help='Add more details to the output.')
def plot_triple(input_path, show_run_details, n, save_stats_path, from_stats, pre_evolution, verbose):
    if not input_path:
        click.secho(f"Please provide input path", fg="red")
        raise click.Abort()
    if from_stats:
        plotter.plot_evolution_from_stats_file(input_path, show_run_details, verbose)
    else:
        plotter.plot_evolution_from_data_file(input_path, n, show_run_details, save_stats_path, pre_evolution, verbose)

@cli.command(name='show-stats', help='show statistics of triple evolutions from a file.')
@click.option('--input-path', default='./current_sample.json', help='input file path of the sample.')
@click.option('--N', is_flag=False, default=0, help='The evolution number in the file to show stats for.')
@click.option('--outcomes', is_flag=True, default=False, help='Show all outcome stats.')
@click.option('--verbose', is_flag=True, default=False, help='Add more details to the output.')
@click.option('--merger-only', is_flag=True, default=False, help='Show only the merger stats.')
@click.option('--manual-exchange', is_flag=True, default=False, help='Show only the manual exchange stats.')
def show_stats(input_path, n, outcomes, verbose, merger_only, manual_exchange):
    pp = pprint.PrettyPrinter()
    click.echo(pp.pprint(get_sample_stats.get_run_details(input_path, n)))
    if outcomes:
        if merger_only:
            click.echo(pp.pprint(get_sample_stats.get_merger_only_file_stats(input_path, verbose=verbose)))
        else:
            click.echo(pp.pprint(get_sample_stats.get_file_stats(input_path, verbose=verbose))) 
    if manual_exchange:
        click.echo(manual_exchange_stats.get_manual_exchange_stats(input_path, verbose=verbose))

@cli.command(name='run-cookbook', help='run a certain cookbook script. See () for more information about this command.')
@click.option('--input-path', default='probab3/cookbook/measure_aHB.py', help='input dir path of the file to run.')
@click.option('--module-name', default='measure_aHB', help='The python module name to run, omit .py suffix.')
@click.option('--plot', is_flag=True, default=False, help='Plot the results.')
@click.option('--option', default=None, help='Provide an option to either run or plot')
@click.option('--merge-files', is_flag=True, default=False, help='Merge the data files and exit.')
@click_log.simple_verbosity_option(logger)
def run_cookbook(input_path, module_name, plot, option, merge_files):
    spec = importlib.util.spec_from_file_location(module_name, input_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    if option:
        options = module.options()
        if option not in options:
            click.secho(f"Option {option} not found in options: {options} \ncontinuing without.", fg="red")
            option = None
    
    if plot:
        if option:
            module.plot(option)
        else:
            module.plot()
    elif merge_files:
        module.merge_files()
    else:
        if option:
            module.run(option)
        else:
            module.run() 

cli.add_command(sample_dbs)

if __name__ == '__main__':
    cli()
