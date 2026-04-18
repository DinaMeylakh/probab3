# Command Line Usage

You can use the probab3 command line tool to generate triple systems data and evolve them.

Make sure you have the command line tool installed. If not, follow the installation instructions in the [README](../README.md).

To view all available commands and their descriptions, type:

```bash
probab3 --help
```

This will return the following list of commands and their descriptions:

```s
Usage: probab3.cmd [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  plot          plot a dynamical binary sequence evolution from a file.
  run-cookbook  run a certain cookbook script.
  sample-dbs     Sample a dynamical binary sequence system in a star cluster.
  show-stats    show statistics of triple evolutions from a file.
  version       show version and exit.
```

You can then continue to inspect probab3 command options by typing:

```bash
probab3 <command> --help
```

Where `<command>` is one of the commands listed above. This will return a list of all the available options for that command and their descriptions.

## Sample an unequal mass black hole system in a two-component star cluster

Let's say you want to evolve a black hole binary in a star cluster, introducing equal mass teritiary all the time and save to "example.json" file.

You can do this by running the `sample-dbs` command. Let's inspect it's options with:

```bash
probab3 sample-dbs --help
```
This will return the following list of options and their descriptions:

```s
Usage: probab3.cmd sample-dbs [OPTIONS] COMMAND [ARGS]...

  Sample a dynamical binary sequence system in a star cluster. Choose a sampling
  method. View the help for a list of methods and the help of each method for
  more details. The arguments for this command should be provided before those
  of the specific method, in the form of: sample-dbs --arg_name=arg_value
  <method> --method_arg_name=method_arg_value.

Options:
  --output_path TEXT           output file path to store the sample.
  --Mc FLOAT                   The Cluster Mass in [MSun]
  --rh FLOAT                   The rh of the cluster in [parsec]
  --rc FLOAT                   The rc of the cluster in [parsec]
  --soft-energy-ratio INTEGER  The ratio between the initial binary energy and
                               the final binary energy that will determine if
                               the binary is ionized.
  --alpha FLOAT                The alpha parameter of the chaotic region.
  --N INTEGER                  The number of dynamical binary sequences to sample.
  --only-inverse-cdf           Only use inverse cdf to sample.
  --help                       Show this message and exit.

Commands:
  equal-mass-basic    Equal mass with mass provided.
  equal-mass-bbh      Evolve stars binary to equal mass BHs binary and...
  hard-soft           Evolve binary until it hardens or softens, use...
  unequal-mass-basic  Sample all masses from present day mass function.
  unequal-mass-bbh    Evolve stars binary to BHs binary and then evolve...
```

Say we want to use the `unequal-mass-bbh` method. We can inspect it further using the following command:
```bash
probab3 sample-dbs unequal-mass-bbh --help
```
This will return the following:

```s
Usage: probab3.cmd sample-dbs unequal-mass-bbh [OPTIONS]

  Evolve stars binary to BHs binary and then evolve it in the cluster. Sample
  all masses from present day mass function.

Options:
  -v, --verbosity LVL  Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  --help               Show this message and exit.
```

We can then run the command with the following options:

```bash
probab3 sample-dbs --output_path "unequal-masses-example.jsonl" --N 10 unequal-mass-bbh
```
This might take a while.
After the command finishes, you will have a file named "unequal-masses-example.json" in your current directory. 

This file contains the metadata of the evolutions and the statistics of the evolution ends inside the whole file.

You can view the metadata of this file by using the `show-stats` command:

```s
Usage: probab3.cmd show-stats [OPTIONS]

  show statistics of triple evolutions from a file.

Options:
  --input-path TEXT  input file path of the sample.
  --N INTEGER        The evolution number in the file to show stats for.
  --outcomes         Show all outcome stats.
  --verbose          Add more details to the output.
  --merger-only      Show only the merger stats.
  --help             Show this message and exit.
```

We can use the `--outcomes` option to view the outcomes of the evolutions in the file:
```bash
probab3 show-stats --input-path unequal-masses-example.jsonl --outcomes
```

This will return the following:

```s
({'ej_binary': 0,
  'ej_fs_collision': 0,
  'ej_fs_em': 0,
  'ej_fs_merger': 0.1,
  'ej_fs_tde': 0,
  'fs_collision': 0,
  'fs_em_merger': 0,
  'fs_merger': 0.1,
  'fs_tde': 0,
  'ims_collision': 0,
  'ims_em_merger': 0,
  'ims_merger': 0,
  'ims_tde': 0,
  'in_cluster_binary': 0.5,
  'ionized_binary': 0.3},
 10,
 {'Mc': 6.289517283862093e+35,
  'alpha': 2.5,
  'end': 'ionized_binary',
  'output_file_path': 'unequal-masses-example.json',
  'rc': 3086000000000000,
  'rejection': True,
  'rh': 30860000000000000,
  'size': 10})
```

Now you can view the evolution of the system by using the `plot` command:

```s
Usage: probab3.cmd plot [OPTIONS]

  plot a dynamical binary sequence evolution from a file.

Options:
  --input-path TEXT       input file path of the sample.
  --show-run-details      Print run details along with the plot.
  --N INTEGER             The run number to plot.
  --save-stats-path TEXT  Save the stats to a file.
  --from-stats            Plot from stats file instead of data file.
  --help                  Show this message and exit.
```

We can use the following command to plot the evolution of the 0th evolution in the file:

```bash
probab3 plot --input-path unequal-masses-example.jsonl --N 0
```

This will open a plot of the evolution of the system.

You can save the statistics of the evolution to a file by using the `--save-stats-path` option:

```bash
probab3 plot --input-path unequal-masses-example.jsonl --N 0 --save-stats-path "stats-of-ev-0.json"
```

This will save the statistics of the 0th evolution to a file named "stats-of-ev-0.json" in your current directory.

You can now plot the evolution by using the `--from-stats` option instead of going through the whole data file:

```bash
probab3 plot --from-stats --input-path "stats-of-ev-0.json"
```

This will open the same plot of the evolution of the system.
