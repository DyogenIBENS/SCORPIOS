# SCORPiOs - Synteny-guided CORrection of Paralogies and Orthologies

[Zenodo code DOI badge] [License badge] [![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.4-brightgreen.svg)](https://snakemake.bitbucket.io)


 SCORPiOs is a synteny-guided gene tree correction snakemake pipeline. SCORPiOs builds optimized gene trees, consistent with a known WGD event, local synteny context, as well as gene sequence evolution.

If you use SCORPiOs, please cite:

*TODO citation*

## Table of content
  - [Installation](#installation)
    - [Installing conda](#conda)
    - [Installing SCORPiOs](#iscorpios)
  - [Usage](#usage)
    - [Running SCORPiOS on example data](#example)
      - [Example 1: Simple SCORPiOs run](#ex1)
      - [Example 2: SCORPiOs in iterative mode](#ex2)
    - [Running SCORPiOS on your data](#data)
    - [Understanding SCORPiOs outputs](#outputs)
      - [Basic](#basic)
      - [Advanced](#advanced)
  - [Authors](#authors)
  - [License](#license)

## <a name="installation"></a>Installation

### <a name="conda"></a>Installing conda

The Miniconda3 package management system manages all SCORPiOs dependencies, including python packages and other softwares.

To install Miniconda3:

- Download Miniconda3 installer for your system [here](https://docs.conda.io/en/latest/miniconda.html)

- Run the installation script: `bash Miniconda3-latest-Linux-x86_64.sh` or `bash Miniconda3-latest-MacOSX-x86_64.sh`, and accept the defaults

- Open a new terminal, run `conda update conda` and press y to confirm updates

### <a name="iscorpios"></a> Installing SCORPiOs

- Clone the repository and go to SCORPiOs root folder
```
git clone https://github.com/DyogenIBENS/SCORPiOs.git
cd SCORPiOs
```

- Create the main conda environment (this may take a few minutes)
```
conda env create -f envs/scorpios.yaml
```

## <a name="usage"></a> Usage

### <a name="examples"></a> Running SCORPiOs on example data

Before any SCORPiOs run, you should:
 - go to SCORPiOs root folder,
 - activate the conda environment with `conda activate scorpios`.

#### <a name="ex1"></a> Example 1: Simple SCORPiOs run

Inputs and parameters to execute SCORPiOs have to be specified in a configuration file.
An example configuration file is provided: [config_example.yaml](config_example.yaml).

This configuration file allows to execute SCORPiOs on toy example data located in [data/example/](data/example/), that you can use as reference for input formats. More details on input file formats are given in [config_example.yaml](config_example.yaml).

In brief, SCORPiOS input files are:
- A single file with the gene trees forest to correct
- A single file with the corresponding multiple alignments
- Gene coordinates files
- A species tree

The only required snakemake arguments to run SCORPiOs are `--configfile` and the `--use-conda` flag. Optionally, you can specify the number of threads via the `--cores` option. For more advanced snakemake usage, you can look at the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html).


To run SCORPiOS on example data:

```
snakemake --configfile config_example.yaml --use-conda --cores 4
```

#### <a name="ex2"></a> Example 2: SCORPiOs in iterative mode

SCORPiOs can run in iterative mode, meaning that SCORPiOS improves gene trees a first time, and then uses the corrected forest as input for a new correction run. Correcting gene trees improves orthologies accuracy, which in turn makes synteny conservation patterns more informative, allowing to better integrate it into the gene tree reconstruction. Usually, a small number of iterations (2-3) suffice to reach convergence.

To run SCORPiOs in iterative mode, you can execute the wrapper bash script `iterate_scorpios.sh`:

```
bash iterate_scorpios.sh --j=example --snake_args="--configfile config_example.yaml"
```

Command-line arguments:

```
Required
--j=JOBNAME, specifies scorpios jobname, should be the same as in config
--snake_args="SNAKEMAKE ARGUMENTS", should at minimum contain --configfile

Optional
--max_iter=MAXITER, maximum number of iterations to run, default=5.
--min_corr=MINCORR, minimum number of corrected sub-trees to continue to the next iteration, default=1.
--starting_iter=ITER, starting iteration, to resume a run at a given iteration, default=1.
```

### <a name="data"></a> Running SCORPiOs on your data
After correct formatting of your input data (see [config_example.yaml](config_example.yaml)), you have to create a new configuration file, using the provided example:

- Copy the example config file `cp config_example.yaml config.yaml`
- Open and edit `config.yaml`

If you want to check your configuration, you can execute a dry-run with `-n`.
```
snakemake --configfile config.yaml --use-conda -n
```

Finally, you can run SCORPiOs as described above.

```
snakemake --configfile config.yaml --use-conda
```

or

```
bash iterate_scorpios.sh --j=myjobname --snake_args="--configfile config.yaml"
```

### <a name="outputs"></a> Understanding SCORPiOs outputs

#### <a name="basic"></a> Basic

All outputs of a SCORPiOs run are stored in a folder named SCORPiOs_jobname, with the jobname specified in the config file.

The main output is the **SCORPiOs corrected gene trees forest**. The commands above generate `SCORPiOs_example/SCORPiOs_corrected_forest_0.nhx` for the simple run and, `SCORPiOs_example/SCORPiOs_corrected_forest_1.nhx` & `SCORPiOs_example/SCORPiOs_corrected_forest_2.nhx` for the iterative run. SCORPiOS adds correction tags in trees that allow to inspect corrections with external gene tree visualisation softwares.

Some intermediary outputs are also stored in different sub-folders (see below for a detailed description). In addition, SCORPiOs writes statistics on key steps of the workflow to the standard output. Thus, to separate output statistics from snakemake logs, you can run:

```
snakemake --configfile config_example.yaml --use-conda >out 2>err
```

or

```
bash iterate_scorpios.sh --j=example --snake_args="--configfile config_example.yaml" >out 2>err
```

#### <a name="advanced"></a> Advanced

To go beyond description statistics printed to the standard output, one might want to further investigate results of a SCORPiOS run, for one or several given gene families. This section quickly introduces key SCORPiOs concepts, in order to better understand intermediary outputs. For a detailed description of SCORPiOs, see [SCORPiOs publication](TODO link).

A gene family in SCORPiOs consists of an outgroup gene and all gene copies in WGD duplicated species. Gene families are first defined in an outgroup/duplicated species orthology table. For each family, SCORPiOs computes a synteny-derived orthology graph, then a constrained tree topology based on synteny, and finally, if necessary, a synteny-aware corrected tree. Through each of these steps, a gene family is identified by the outgroup gene name.

The Orthology table is stored in a single file in the `Families/` sub-folder. Raw synteny-predicted orthologies are stored in a single file in `Synteny/`. Predicted orthology groups based on community detection in synteny graphs are stored in a single file in `Graphs/`, along with a summary of the community detection step. Finally, `Corrections/` stores two files, one detailing trees vs synteny consistency and another with the list of successfully corrected trees.

Several tags such as the name of the corrected WGD, the outgroup species used and SCORPiOS iteration number are added to each output file, in order to identify it. When run in simple mode, the iteration number is set to 0.

Additional files can be saved if specified in the configuration file, see [config_example.yaml](config_example.yaml) for details.

## <a name="authors"></a> Authors
* [**Elise Parey**](parey@biologie.ens.fr)
* **Hugues Roest Crollius**
* **Camille Berthelot**

## <a name="license"></a> License

TODO
