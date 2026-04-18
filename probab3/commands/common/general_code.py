import os
import json
import dataclasses
import click
import zipfile
import glob
import pandas
import numpy as np
from typing import Any
from pathlib import Path
from numpyencoder import NumpyEncoder

os.makedirs('logs', exist_ok=True)
LOG_FILE_NAME = f"logs/probab3-{os.getpid()}.log"

def append_dict_to_file(file_path, output):
    out = [output]

    if os.path.exists(file_path):
        with open(file_path, "r+") as outfile:
            infile = outfile.read()
            if infile:
                infile_list = json.loads(infile)
                infile_list.append(output)
                out = infile_list
    
    with open(file_path, "w+") as outfile: 
        outfile.write(json.dumps(out, cls=EnhancedJSONEncoder))

def fast_slim_append_dict_to_file(file_path, output):
    
    if os.path.exists(file_path):
        with open(file_path, "rb+") as outfile:
            if os.stat(file_path).st_size > 0:
                outfile.seek(-1, os.SEEK_END)
                outfile.truncate()
        
        with open(file_path, "a") as outfile:
            if os.stat(file_path).st_size > 0:
                outfile.write(',')
                outfile.write(json.dumps(output, cls=SlimEnhancedJSONEncoder))
                outfile.write(']')

    else:
        out = [output]
        with open(file_path, "w+") as outfile: 
            outfile.write(json.dumps(out, cls=SlimEnhancedJSONEncoder))

def fast_slim_append_dict_to_jsonl_file(file_path, output):
    
    if os.path.exists(file_path):
        with open(file_path, "a") as outfile:
            outfile.write(json.dumps(output, cls=SlimEnhancedJSONEncoder))
            outfile.write('\n')

    else:
        with open(file_path, "w+") as outfile: 
            outfile.write(json.dumps(output, cls=SlimEnhancedJSONEncoder))
            outfile.write('\n')


def read_dict_from_file(file_path):
    is_jsonl = False
    try:
        if zipfile.is_zipfile(file_path):
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                with zip_ref.open(zip_ref.namelist()[0]) as outfile:
                    is_jsonl = outfile.name.endswith(".jsonl")
                    infile_list = pandas.read_json(outfile, lines=is_jsonl)

        else:
            infile_list = []
            with open(file_path, "r+") as outfile:
                is_jsonl = outfile.name.endswith(".jsonl")
                infile_list = pandas.read_json(outfile, lines=is_jsonl)

        return infile_list
    except Exception as err:
        click.echo(f"Got ERROR {err} when trying to read file {file_path}.")
        if not is_jsonl:
            click.echo(f"Trying to read file as JSON")
            return pandas.DataFrame(read_dict_from_json_file(file_path))
        else:
            raise err


def read_dict_from_json_file(file_path):
    try:
        if zipfile.is_zipfile(file_path):
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                with zip_ref.open(zip_ref.namelist()[0]) as outfile:
                    infile_dict = json.load(outfile)

        else:
            infile_dict = {}
            with open(file_path, "r+") as outfile:
                infile_dict = json.load(outfile)

        return infile_dict
    except Exception as err:
        click.echo(f"Got ERROR {err} when trying to read file {file_path}.")
        raise err

class IgnoreNoneValues(dict):
    def __setitem__(self, __key: Any, __value: Any) -> None:
        if __value is not None:
            return super().__setitem__(__key, __value)

class EnhancedJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if dataclasses.is_dataclass(o):
            return dataclasses.asdict(o, dict_factory=IgnoreNoneValues) 
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)

class SlimEnhancedJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if dataclasses.is_dataclass(o):
            # dont save the keys, just the values
            return tuple(dataclasses.asdict(o).values())
        return super().default(o)


def convert_zip_to_jsonl(input_file):

    new_output_file = str(Path(input_file).with_suffix("")) + ".jsonl"

    with zipfile.ZipFile(input_file, 'r') as zip_ref:
        with zip_ref.open(zip_ref.namelist()[0]) as input_file_inner:
            with open(new_output_file, 'w+') as f_out:
                # Parse the JSON file iteratively
                objects = json.loads(input_file_inner.read())

                for obj in objects:
                    # Write each object to the output file on a new line
                    f_out.write(json.dumps(obj))
                    f_out.write('\n')

    jsonl_zip_filename = str(Path(input_file).with_suffix("")) + "_jsonl.zip"

    with zipfile.ZipFile(jsonl_zip_filename, 'w', compression=zipfile.ZIP_DEFLATED) as zip_ref:
        zip_ref.write(new_output_file, arcname=os.path.basename(new_output_file))

    os.remove(new_output_file)

def convert_json_to_jsonl(input_file):
    
    new_output_file = str(Path(input_file).with_suffix("")) + ".jsonl"
    
    with open(input_file, 'r') as input_file_inner:
        with open(new_output_file, 'w+') as f_out:
            # Parse the JSON file iteratively
            objects = json.loads(input_file_inner.read())
    
            for obj in objects:
                # Write each object to the output file on a new line
                f_out.write(json.dumps(obj))
                f_out.write('\n')


def merge_jsonl_files(input_files_pattern, output_file):

    input_files = glob.glob(input_files_pattern)
    click.echo(f"Found {len(input_files)} files to merge. Saving to {output_file}")

    if len(input_files) == 0:
        click.echo(f"No files found with pattern {input_files_pattern}. Not saving anything.")
        return

    if not os.path.exists(output_file):
        with open(output_file, "w") as outfile: 
            pass

    for index, input_file in enumerate(input_files):
        click.echo(f"Processing file {index+1} of {len(input_files)}: {input_file}")

        with open(input_file, "r") as infile:
            infile_data = infile.read()

            with open(output_file, "a") as outfile: 
                outfile.write(infile_data)
    
    click.echo(f"Done processing {len(input_files)} files. Output file: {output_file}")


def save_numpy_dict_to_json_file(file_path, output):
    with open(file_path, "w+") as outfile: 
        json.dump(output, outfile, cls=NumpyEncoder, indent=4)

def get_file_lines_len(file_path):
    with open(file_path, "r") as file:
        return sum(1 for line in file)