# Cookbooks

Cookbooks are practical examples demonstrating how to use the `probab3` package for sampling 3-body properties in star clusters. They provide ready-to-use scripts and notebooks that showcase different sampling scenarios, such as equal mass and unequal mass triple systems in different environments. The cookbooks serve as starting templates for users to learn how to use `probab3` for their own sampling projects, and they can be run directly via the command line or imported as Python modules.

See the cookbook folder for available cookbooks and their descriptions.

## How to Implement

To create a functional cookbook that can be run via the `probab3 run-cookbook` command, you must implement the following functions in your Python script: run(), plot(), options(), and merge_files().

See for example the `equal_mass_in_sub_cluster.py` cookbook in the `equal_mass` folder for a template implementation of these functions.

### Required Functions

#### `run()`
The main sampling function. This performs the actual sampling of triple systems.
You can invoke any sampling function from the `probab3.commands` module, such as `sample_dbs`, `sample_dbs_in_cluster`, or `sample_dbs_with_tertiary`. The `option` parameter can be used to specify different sampling configurations.

#### `plot()`
Function to visualize the sampling results.

### Optional Functions

#### `options()`
Return a list of available options for `run` and `plot` functions.

#### `merge_files()`
Function to merge multiple data files (useful for parallel runs).

## How to Run Cookbooks

Use the `probab3 run-cookbook` command to execute cookbooks:

Basic Run:
```bash
probab3 run-cookbook --input-path equal_mass/equal_mass_in_sub_cluster.py --module-name equal_mass_in_sub_cluster
```

Run with Plotting:
```bash
probab3 run-cookbook --input-path unequal_mass/unequal_mass_bbh.py --module-name unequal_mass_bbh --plot
```

Run with Specific Option:
```bash
probab3 run-cookbook --input-path equal_mass/equal_mass_in_sub_cluster.py --module-name equal_mass_in_sub_cluster --option fast
```

Merge Files:
```bash
probab3 run-cookbook --input-path unequal_mass/unequal_mass_bbh.py --module-name unequal_mass_bbh --merge-files
```

Plot with Option:
```bash
probab3 run-cookbook --input-path unequal_mass/unequal_mass_bbh.py --module-name unequal_mass_bbh --plot --option masses
```
