# SCORPiOs - Synteny-guided CORrection of Paralogies and Orthologies

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.4-brightgreen.svg)](https://snakemake.bitbucket.io)


 SCORPiOs is a **synteny-guided gene tree correction pipeline** for clades that have undergone a whole-genome duplication event. SCORPiOs identifies gene trees where the whole-genome duplication is **missing** or **incorrectly placed**, based on the genomic locations of the duplicated genes across the different species. SCORPiOs then builds an **optimized gene tree** consistent with the known WGD event, the species tree, local synteny context, as well as gene sequence evolution.

 SCORPiOs is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. SCORPiOs takes as input either gene trees or multiple alignments, and outputs the corresponding optimized gene trees.

 For a complete description of SCORPiOs, see our preprint: https://www.biorxiv.org/content/10.1101/2020.01.30.926915v1.full

## Table of content
  - [Installation](#installation)
    - [Installing conda](#installing-conda)
    - [Installing SCORPiOs](#installing-scorpios)
  - [Usage](#usage)
    - [Setting up your working environment for SCORPiOs](#setting-up-your-working-environment-for-scorpios)
    - [Running SCORPiOs on example data](#running-scorpios-on-example-data)
      - [Example 1: Simple SCORPiOs run](#example-1-simple-scorpios-run)
      - [Example 2: Iterative SCORPiOs run](#example-2-iterative-scorpios-run)
    - [Running SCORPiOS on your data](#running-scorpios-on-your-data)
      - [Data preparation and formatting](#data-preparation-and-formatting)
      - [Preparing your configuration file](#preparing-your-configuration-file)
      - [Running SCORPiOs](#running-scorpios)
    - [Understanding SCORPiOs outputs](#understanding-scorpios-outputs)
      - [Basic](#basic)
      - [Advanced](#advanced)
    - [Complex configurations](#complex-configurations)
  - [Authors](#authors)
  - [License](#license)
  - [References](#references)

## Installation

### Installing conda

The Miniconda3 package management system manages all SCORPiOs dependencies, including python packages and other software.

To install Miniconda3:

- Download Miniconda3 installer for your system [here](https://docs.conda.io/en/latest/miniconda.html)

- Run the installation script: `bash Miniconda3-latest-Linux-x86_64.sh` or `bash Miniconda3-latest-MacOSX-x86_64.sh`, and accept the defaults

- Open a new terminal, run `conda update conda` and press `y` to confirm updates

### Installing SCORPiOs

- Clone the repository and go to SCORPiOs root folder
```
git clone https://github.com/DyogenIBENS/SCORPIOS.git
cd SCORPIOS
```

- Create the main conda environment (solving dependencies may take a while, ~ up to an hour)
```
conda env create -f envs/scorpios.yaml
```

## Usage

### Setting up your working environment for SCORPiOs

Before any SCORPiOs run, you should:
 - go to SCORPiOs root folder,
 - activate the conda environment with `conda activate scorpios`.

### Running SCORPiOs on example data

Before using SCORPiOs on your data, we recommend running a test with our example data to ensure that installation was successful and to get familiar with the pipeline, inputs and outputs.

#### Example 1: Simple SCORPiOs run

Inputs and parameters to execute SCORPiOs have to be specified in a YAML configuration file.
An example configuration file is provided: [config_example.yaml](config_example.yaml). This configuration file executes SCORPiOs on toy example data located in [data/example/](data/example/), that you can use as reference for input formats.

The only required snakemake arguments to run SCORPiOs are `--configfile` and the `--use-conda` flag. Optionally, you can specify the number of threads via the `--cores` option. For more advanced options, you can look at the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

To run SCORPiOs on example data:

```
snakemake --configfile config_example.yaml --use-conda --cores 4
```

The following outputs should be generated:
[insert outputs here]

We explain how to interpret outputs below [link to understanding SCORPiOs outputs].

#### Example 2: Iterative SCORPiOs run

SCORPiOs can run in iterative mode, meaning that SCORPiOs improves gene trees a first time, and then uses the corrected set of gene trees again as input for a new correction run. Correcting gene trees improves orthologies accuracy, which in turn makes synteny conservation patterns more informative, allowing to better integrate it into the gene tree reconstruction. Usually, a small number of iterations (2-3) suffice to reach convergence.

To run SCORPiOs in iterative mode on example data, execute the wrapper bash script `iterate_scorpios.sh`:

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

### Running SCORPiOs on your data

#### Data preparation and formatting
SCORPiOs is a flexible gene tree correction pipeline: it can either start from a set of precomputed, phylogeny-reconciled gene trees, or build one from a set of gene multiple aligments using [TreeBeST](link to Treebest paper).

SCORPiOs input files are:
- A single file with a set of phylogeny-reconciled gene trees in NHX format (extended Newick format, see [example](example)) **OR** a genes-to-species mapping file, if working from gene alignments (see [example](example))
- A single file with the corresponding gene multiple alignments in FASTA format (can be compressed with gzip or bzip2)
- Gene coordinates files for each species in BED format (see [example](example))
- A species tree in PhylTree format (see [example](example))

Detailed information on input files and formats can be found in [config_example.yaml](config_example.yaml).

#### Preparing your configuration file
Once your data is formatted correctly, you have to create a new configuration file for your SCORPiOs run, using the provided example:

- Copy the example config file `cp config_example.yaml config.yaml`
- Open and edit `config.yaml` to specify paths, files and parameters for your data

To check your configuration, you can execute a dry-run with `-n`.
```
snakemake --configfile config.yaml --use-conda -n
```

#### Running SCORPiOs
Finally, you can run SCORPiOs as described above:

```
snakemake --configfile config.yaml --use-conda
```

or in iterative mode, assuming the jobname is set to 'myjobname' in the new config file:

```
bash iterate_scorpios.sh --j=myjobname --snake_args="--configfile config.yaml"
```

### Understanding SCORPiOs outputs

#### Basic

All outputs from SCORPiOs are stored in a folder named SCORPiOs_jobname (jobname as specified in the configuration file).

The main output is the **SCORPiOs-optimized gene trees**. Gene trees are provided as a single file in NHX format. SCORPiOs tags corrected nodes in the gene trees to allow easy inspection using tree visualisation softwares. We recommand the [ETE Toolkit](http://etetoolkit.org/) for tree visualisation.

The commands above generate:
- `SCORPiOs_example/SCORPiOs_corrected_forest_0.nhx` for the simple run
- `SCORPiOs_example/SCORPiOs_corrected_forest_2_with_tags.nhx` for the iterative run.

Outputs are suffixed with a digit representing the iteration number. This number is set to 0 in simple mode and starts at 1 in iterative mode.

Some intermediary outputs are also stored in different sub-folders (see below for a detailed description). In addition, SCORPiOs writes statistics on key steps of the workflow to the standard output. Thus, to separate output statistics from snakemake logs, you can run:

```
snakemake --configfile config_example.yaml --use-conda >out 2>err
```

or

```
bash iterate_scorpios.sh --j=example --snake_args="--configfile config_example.yaml" >out 2>err
```

#### Advanced

Beyond description statistics printed to the standard output, you may want to investigate the detailed results of SCORPiOs for one or several given gene families. This section introduces a few key concepts of SCORPiOs, in order to better understand intermediary outputs.

A gene family in SCORPiOs consists of a non-duplicated outgroup gene and all potential orthologous gene copies in WGD-duplicated species, based on the uncorrected gene trees. For each family, SCORPiOs computes a synteny-derived orthology graph, then a constrained tree topology based on synteny, and finally, if necessary, a synteny-aware corrected tree. Through each of these steps, a gene family is identified by the outgroup gene name.

The orthology relationships between genes are stored in a single file in the `Families/` sub-folder. Raw synteny-predicted orthologies are stored in a single file in `Synteny/`. Predicted orthology groups based on community detection in synteny graphs are stored in a single file in `Graphs/`, along with a summary of the community detection step. Finally, `Corrections/` stores two files, one detailing trees vs synteny consistency and another with the list of successfully corrected trees.

Several tags such as the name of the corrected WGD, the outgroup species and SCORPiOs iteration number are added to each output file, in order to precisely identify outputs in case of complex configurations (see below).

Additional files can be saved if specified in the configuration file, see [config_example.yaml](config_example.yaml) for details.

### Complex configurations

SCORPiOs can correct gene trees that contain more than one whole-genome duplication event. In this case, each WGD is treated independently, starting from the more recent one (closer to the leaves) going up towards the more ancient one (closer to the root). If the WGDs are nested, the subtrees from the more recent events are ignored while correcting for the older WGD event(s), and reinserted after correction using their outgroup as a branching point.

Several WGDs can be specified in the configuration file. //Details details//

## Authors
* [**Elise Parey**](mailto:elise.parey@bio.ens.psl.eu)
* **Alexandra Louis**
* **Hugues Roest Crollius**
* **Camille Berthelot**

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
- [LICENSE GPLv3](LICENSE.txt)

## References

SCORPiOs uses the following tools to build and test gene trees:

- [ProfileNJ](https://github.com/maclandrol/profileNJ): Noutahi et al. (2016) Efficient Gene Tree Correction Guided by Genome Evolution. PLOS ONE, 11, e0159559.

- [PhyML](http://www.atgc-montpellier.fr/phyml/): Guindon et al. (2010) New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0. Syst Biol, 59, 307–321.

- [TreeBeST](https://github.com/Ensembl/treebest): Vilella et al. (2009) EnsemblCompara GeneTrees: Complete, duplication-aware phylogenetic trees in vertebrates. Genome Res., 19, 327–335.

- [CONSEL](https://github.com/shimo-lab/consel): Shimodaira and Hasegawa (2001) CONSEL: for assessing the confidence of phylogenetic tree selection. Bioinformatics, 17, 1246–1247.
