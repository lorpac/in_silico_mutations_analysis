# Analysis of in-silico amino-acid mutations of PSD95$^{\text{pdz3}}$

This repository contains the code associated with the publications:
- Classification of in-silico amino-acid mutations based on the scale of neighborhoods rearrangement they provoke: **ADD REF**
- GCAT network: **ADD REF**

## Prerequisites

The required packages are listed in the file [requirements.txt](https://github.com/lorpac/amino_acid_network/blob/master/requirements.txt). To install the requirements with [pip](https://pypi.org/project/pip/), type:

```
pip install -r requirements.txt
```

You also need Rodigo Gilardi's class [biographs](https://github.com/rodogi/biographs). Clone it and do not forget to add the class folder (`biograps/biographs`) to you Python path. Alternatively, if you are using a pip virtual environment `env`, you can manually add it to your packages by saving a copy in the directory `env/Lib/site-packages`.

If you are not familiar with virtual environments, [here](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/26/python-virtual-env/) you will find the detailed instructions on how to create and manage virtual environments with `pip`.

## Set-up the analysis

To set up the parameters of you analysis, you need to create a configuration file following the template of `config_template.json` and save it as `config.json`. In the configuration file, the following parameters must be given (in any order):

- `name`: a name of your analysis
- `pdb_wt`: either a `.pdb` or `.ent` file name found in `data/`, or a valid PDB identifier. Note: in the second case a internet connection will be needed to download the required pdb structure from the [RCSB database](http://www.rcsb.org/).
- `pdbs_wt_path`: name of the subfolder in `data/` containing the `.pdb` or `.ent` file `pdb_wt`. If null, the `pdbs_wt_path` is set to be queal to `pdbs`.
- `pdbs_mut_path`: name of the subfolder in `data/` containing the `.pdb` or `.ent` files of the mutants. The mutant files must be named using the following template: `mut_XiY.pdb` where `i` is the mutated position, `X` is the amino-acid type in the WT structure and `Y` is the amino-acid type in the mutated structure (one-letter code). Only single-amino acid mutants are allowed.
- `cutoff`: cutoff distance to determine the connectivity between amino acids.
- `dim`: type of links to consider. Valid options are: `all` (or an empty string, all links), `1D` (between first neighbors in the sequence), `2D` (between amino acids belonging to the same secondary structure element), `3D` (between amino acids belonging to the same chain but to the different secondary structure elements), `4D` (between amino acids belonging to different chains), `3-4D` (3D and 4D), `1-2D` (3D and 4D).
- `select positions`: boolean.
  - If `select positions` is true: `start` and `stop` are the sequence positions to start and stop the analysis

## Run the analysis

The code for running the analysis is provided in the Jupyter Notebook `in_silico_mutation_analysis.ipynb`.

The outputs are saved in a folder named `results/name/` where `name` is the job title defined in the configuration file.

## Authors

Lorenza Pacini - [lorpac](https://github.com/lorpac)

## How to cite

If you use this code please cite:

- Classification of in-silico amino-acid mutations based on the scale of neighborhoods rearrangement they provoke: **ADD REF**
- GCAT network: **ADD REF**
  
## Licence

This code is available under the [CeCILL](http://cecill.info/) licence. Please see `LICENCE.txt` for details.



