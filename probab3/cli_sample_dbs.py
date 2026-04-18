import click
import time
import click_log
import logging

from probab3.commands.sample_dbs import dbs_sampler_unequal_mass, dbs_sampler_equal_mass 
from probab3.commands.common.constants import *
from probab3.commands.common.general_code import *
from probab3.commands.common.formulas.general import EB, LB_from_EBeB

from colorama import init as colorama_init
logging.basicConfig(level=logging.INFO, filename=LOG_FILE_NAME, filemode='a', format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(PROJECT_NAME)
colorama_init()

@click.group(name='sample-dbs', help='Sample a dynamical binary sequence system in a star cluster. Choose a sampling method. ' 
                                    'View the help for a list of methods and the help of each method for more details. '
                                    'The arguments for this command should be provided before those of the specific method, '
                                    'in the form of: sample-dbs --arg_name=arg_value <method> --method_arg_name=method_arg_value.')
@click.option('--output-path', default='./current_sample.json', help='output file path to store the sample.')
@click.option('--Mc', default=10**(5.5), help='The Cluster Mass in [MSun]')
@click.option('--rh', default=1.0, help='The rh of the cluster in [parsec]')
@click.option('--rc', default=0.1, help='The rc of the cluster in [parsec]')
@click.option('--alpha', default=2.5, help='The alpha parameter of the chaotic region.')
@click.option('--N', default=1, help='The number of dynamical binary sequences to sample.')
@click.option('--only-inverse-cdf', default=False, is_flag=True, help='Only use inverse cdf to sample.')
@click.pass_context
def sample_dbs(ctx, output_path: str, alpha:float, mc:int, rh:int, rc:int, n:int, only_inverse_cdf:bool):
    ctx.ensure_object(dict)
    ctx.obj['output_path'] = output_path
    ctx.obj['alpha'] = alpha
    ctx.obj['mc'] = mc
    ctx.obj['rh'] = rh
    ctx.obj['rc'] = rc
    ctx.obj['n'] = n
    ctx.obj['only_inverse_cdf'] = only_inverse_cdf


@sample_dbs.command(name='one-equal-mass-triple', help='One equal mass triple with mass provided.')
@click.option('--m', default=20.0, help='The ma mass in [MSun]')
@click.option('--E0', default=EB(20*MSun, 20*MSun, 10*AU), help='The initial energy of the triple in SI units.')
@click.option('--L0', default=LB_from_EBeB(20*MSun, 20*MSun, EB(20*MSun, 20*MSun, 10*AU), 0.5), help='The initial angular momentum of the triple in SI units.')
@click_log.simple_verbosity_option(logger)
@click.pass_context
def sample_one_equal_triple(ctx, m: int, e0: int, l0: int):
    start = time.time()
    dbs_sampler_equal_mass.one_triple_equal_mass_sample(m1=m*MSun,  
                                                  alpha=ctx.obj.get('alpha'), 
                                                  E0=e0, 
                                                  L0=l0,
                                                  output_file_path=ctx.obj.get('output_path'),
                                                  size=ctx.obj.get('n'),
                                                  rejection=not ctx.obj.get('only_inverse_cdf'))
    end = time.time()
    click.secho(f'time: {end - start}', fg='yellow')


@sample_dbs.command(name='one-unequal-mass-triple', help='One unequal mass triple with masses provided.')
@click.option('--m1', default=20.0, help='The first mass in [MSun]')
@click.option('--m2', default=10.0, help='The second mass in [MSun]')
@click.option('--m3', default=15.0, help='The third mass in [MSun]')
@click.option('--E0', default=EB(20*MSun, 20*MSun, 10*AU), help='The initial energy of the triple in SI units.')
@click.option('--L0', default=LB_from_EBeB(20*MSun, 20*MSun, EB(20*MSun, 20*MSun, 10*AU), 0.5), help='The initial angular momentum of the triple in SI units.')
@click_log.simple_verbosity_option(logger)
@click.pass_context
def sample_one_equal_triple(ctx, m1: float, m2: float, m3:float, e0: int, l0: int):
    start = time.time()
    dbs_sampler_unequal_mass.one_triple_unequal_mass_sample(m1=m1*MSun,  
                                                            m2=m2*MSun,
                                                            m3=m3*MSun,
                                                            alpha=ctx.obj.get('alpha'), 
                                                            E0=e0, 
                                                            L0=l0,
                                                            output_file_path=ctx.obj.get('output_path'),
                                                            size=ctx.obj.get('n'),
                                                            rejection=not ctx.obj.get('only_inverse_cdf'))
    end = time.time()
    click.secho(f'time: {end - start}', fg='yellow')


@sample_dbs.command(name='equal-mass-basic', help='Equal mass with mass provided in a cluster.')
@click.option('--m1', default=20.0, help='The ma mass in [MSun]')
@click_log.simple_verbosity_option(logger)
@click.pass_context
def sample_dbs_basic(ctx, m1: int):
    start = time.time()
    dbs_sampler_equal_mass.basic_equal_mass_sample(m1=m1*MSun, m2=m1*MSun, m3=m1*MSun, 
                                                  alpha=ctx.obj.get('alpha'), 
                                                  Mc=ctx.obj.get('mc')*MSun, 
                                                  rh=ctx.obj.get('rh')*parsec,
                                                  rc=ctx.obj.get('rc')*parsec,
                                                  output_file_path=ctx.obj.get('output_path'),
                                                  size=ctx.obj.get('n'),
                                                  rejection=not ctx.obj.get('only_inverse_cdf'))
    end = time.time()
    click.secho(f'time: {end - start}', fg='yellow')

@sample_dbs.command(name='unequal-mass-basic', help='Sample all masses from present day mass function.')
@click_log.simple_verbosity_option(logger)
@click.pass_context
def sample_dbs_unequal_mass(ctx):
    start = time.time()
    dbs_sampler_unequal_mass.unequal_mass_sample(alpha=ctx.obj.get('alpha'), 
                                                Mc=ctx.obj.get('mc')*MSun, 
                                                rh=ctx.obj.get('rh')*parsec,
                                                rc=ctx.obj.get('rc')*parsec,
                                                output_file_path=ctx.obj.get('output_path'),
                                                size=ctx.obj.get('n'),
                                                rejection=not ctx.obj.get('only_inverse_cdf'))
    end = time.time()
    click.secho(f'time: {end - start}', fg='yellow') 

@sample_dbs.command(name='unequal-mass-bbh', help='Evolve stars binary to BHs binary and then evolve it in the cluster. Sample all masses from present day mass function.')
@click_log.simple_verbosity_option(logger)
@click.pass_context
def sample_dbs_unequal_mass(ctx):
    start = time.time()
    dbs_sampler_unequal_mass.unequal_mass_bbh_sample(alpha=ctx.obj.get('alpha'), 
                                                    Mc=ctx.obj.get('mc')*MSun, 
                                                    rh=ctx.obj.get('rh')*parsec,
                                                    rc=ctx.obj.get('rc')*parsec,
                                                    output_file_path=ctx.obj.get('output_path'),
                                                    size=ctx.obj.get('n'),
                                                    rejection=not ctx.obj.get('only_inverse_cdf'))
    end = time.time()
    click.secho(f'time: {end - start}', fg='yellow') 

@sample_dbs.command(name='equal-mass-bbh', help='Evolve stars binary to equal mass BHs binary and evolve as equal mass triple.')
@click_log.simple_verbosity_option(logger)
@click.pass_context
def sample_dbs_to_end(ctx):
    start = time.time()
    dbs_sampler_equal_mass.ms_to_BH_to_end_sample(alpha=ctx.obj.get('alpha'), 
                                                 Mc=ctx.obj.get('mc')*MSun, 
                                                 rh=ctx.obj.get('rh')*parsec,
                                                 rc=ctx.obj.get('rc')*parsec,
                                                 output_file_path=ctx.obj.get('output_path'),
                                                 size=ctx.obj.get('n'),
                                                 rejection=not ctx.obj.get('only_inverse_cdf'))
    end = time.time()
    click.secho(f'time: {end - start}', fg='yellow')

if __name__ == '__main__':
    sample_dbs(obj={})