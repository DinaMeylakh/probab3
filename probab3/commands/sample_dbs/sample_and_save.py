import click
import logging
from colorama import init as colorama_init
from colorama import Fore, Style
import datetime
import time
import traceback
from probab3.commands.common.general_code import * 


logging.basicConfig(level=logging.INFO, filename=LOG_FILE_NAME, filemode='a', format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)
colorama_init()


def sample_and_save(func):

    def wrapper(*args, **kwargs):
        sample_dict = kwargs.copy()
        output_file_path = sample_dict.get("output_file_path")
        size = sample_dict.get("size")
        
        if not output_file_path:
            raise Exception("No output file specified. Aborting.")

        if os.path.exists(output_file_path):
            file_current_size = get_file_lines_len(output_file_path)
            if file_current_size < size:
                size = size - file_current_size
                click.secho(f"File {output_file_path} already exists with {file_current_size} lines. Continuing with {size} more.", fg="yellow")
            else:
                click.secho(f"File {output_file_path} already exists with {file_current_size} lines. Aborting.", fg="red")
                return

        batch_start_time = time.time()
        click.secho(f"Starting batch of size {size} with PID {os.getpid()} and log file {LOG_FILE_NAME}")
        logger.info(f"{batch_start_time}::{Fore.MAGENTA}Starting batch of size {size} with {sample_dict}{Style.RESET_ALL}")
        with click.progressbar(range(size)) as size_range:
            for i in size_range:
                click.secho(f"{str(datetime.datetime.now())}:: Starting evolution {i+1}/{size}", fg="cyan", bold=True)
                logger.info(f"{Fore.MAGENTA}{str(datetime.datetime.now())}::Starting evolution {i+1}/{size}{Style.RESET_ALL}")
                try:
                    start_time = time.time()
                    logger.info(f"{Fore.MAGENTA}{str(datetime.datetime.now())}::Starting {func.__name__}: {i+1}/{size}{Style.RESET_ALL}")
                    sample_updates = func(*args, **kwargs)
                    sample_dict.update(sample_updates)
                    end_time = time.time()
                    end = sample_dict.get("end", "unknown")
                    click.secho(f"{Fore.MAGENTA}{str(datetime.datetime.now())}::evolution ended with {end} after {end_time - start_time} seconds{Style.RESET_ALL}")
                    logger.info(f"{Fore.MAGENTA}{str(datetime.datetime.now())}::evolution ended with {end} after {end_time - start_time} seconds{Style.RESET_ALL}")
                    
                    fast_slim_append_dict_to_jsonl_file(output_file_path, sample_dict)
                    export_end_time = time.time()
                    click.secho(f"{Fore.MAGENTA}{str(datetime.datetime.now())}::Exporting evolution {i+1}/{size} to file {output_file_path} took {export_end_time - end_time} seconds{Style.RESET_ALL}", bold=True)
                    logger.info(f"{Fore.MAGENTA}{str(datetime.datetime.now())}::Exporting evolution {i+1}/{size} to file {output_file_path} took {export_end_time - end_time} seconds {Style.RESET_ALL}")

                except Exception as err:
                    error_time = time.time()
                    logger.error(f"{Fore.RED}{str(datetime.datetime.now())}::{error_time}::Got ERROR {err} when trying to evolve binary."
                                 f"Traceback:\n{traceback.format_exc()}\n"
                                 f"Continuing to next..{Style.RESET_ALL}")
                    click.echo(f"{Fore.RED}{str(datetime.datetime.now())}::{error_time}::Got ERROR {err} when trying to evolve binary. "
                               f"Continuing to next..{Style.RESET_ALL}")
                    continue
            
        batch_end_time = time.time()
        sample_dict.pop("evolution", None)
        sample_dict.pop("pre_evolution", None) 
        logger.info(f"{str(datetime.datetime.now())}::{Fore.MAGENTA}Ended batch of size {size} with {sample_dict} after {batch_end_time - batch_start_time} seconds {Style.RESET_ALL}")
        click.echo(f"{str(datetime.datetime.now())}::{Fore.MAGENTA}Ended batch of size {size} with {sample_dict} after {batch_end_time - batch_start_time} seconds {Style.RESET_ALL}")

    return wrapper
