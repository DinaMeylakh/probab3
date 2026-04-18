# Import Usage
Fist see if the package is installed by running the following command in your python environment:
```python
import probab3
```

If it returned no errors, you can now use the package in your python scripts.

For example, to sample 5 unequal mass binary black hole systems in a two-component star cluster, you can use the following code:
```python
from probab3.commands.sample_dbs import dbs_sampler_unequal_mass
from probab3.commands.pre_calculation.cluster_params.cluster_params_relations import select_rh_for_Mc
from probab3.commands.common.constants import MSun, parsec



sample_path = f"samples/example_run.jsonl"

dbs_sampler_unequal_mass.unequal_mass_bbh_sample(alpha=2.5, 
                                                 Mc=(10**6)*MSun, 
                                                 rh=select_rh_for_Mc((10**6))*parsec, 
                                                 rc=0.5*select_rh_for_Mc((10**6))**parsec,
                                                 size=5,
                                                 output_file_path=sample_path)

```

This will create a file named `example_run.jsonl` in the `samples` directory of the package. The file will contain the metadata of the evolutions and the statistics of the evolution ends inside the whole file.

