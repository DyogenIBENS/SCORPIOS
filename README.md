# SCORPiOs - Synteny-guided CORrection of Paralogies and Orthologies

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3727519.svg)](https://doi.org/10.5281/zenodo.3727519) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Snakemake](https://img.shields.io/badge/snakemake-≥6.6.1-brightgreen.svg)](https://snakemake.bitbucket.io) [![Documentation Status](https://readthedocs.org/projects/scorpios/badge/?version=latest)](https://scorpios.readthedocs.io/en/latest/?badge=latest)

 SCORPiOs is a **synteny-guided gene tree correction pipeline** for clades that have undergone a whole-genome duplication event. SCORPiOs identifies gene trees where the whole-genome duplication is **missing** or **incorrectly placed**, based on the genomic locations of the duplicated genes across the different species. SCORPiOs then builds an **optimized gene tree** consistent with the known WGD event, the species tree, local synteny context, as well as gene sequence evolution.

 SCORPiOs is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. SCORPiOs takes as input either gene trees or multiple alignments, and outputs the corresponding optimized gene trees.

 <br />


| :sparkles:  New in SCORPiOs 2.0.0: LORelEi (Lineage-specific Ohnolog Resolution Extension)|
|:---------------------------|
| SCORPiOs LORelEi analyzes sequence-synteny conflicts in gene trees and diagnose cases of delayed meiosis resolution following WGD. To learn how to use SCORPiOs and SCORPiOs LORelEi, take a look at [SCORPiOs documentation](https://scorpios.readthedocs.io/en/latest/)!   |

<br />

 ![SCORPiOs illustrated](https://github.com/DyogenIBENS/SCORPIOS/blob/master/doc/img/scorpios_illustrated.png)

If you use SCORPiOs, please cite:

Parey E, Louis A, Cabau C, Guiguen Y, Roest Crollius H, Berthelot C, Synteny-guided resolution of gene trees clarifies the functional impact of whole genome duplications, Molecular Biology and Evolution, msaa149, https://doi.org/10.1093/molbev/msaa149.

# Quick start

**Below is a quick start guide to using SCORPiOs, we recommend reading [SCORPiOs documentation](https://scorpios.readthedocs.io/en/latest/) for detailed instructions.**

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
      - [Preparing your configuration file](#preparing-your-configuration-file)
      - [Running SCORPiOs](#running-scorpios)
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

- Create the main conda environment.

  We recommend using [Mamba](https://github.com/mamba-org/mamba) for a faster installation:
  ```
  conda install -c conda-forge mamba
  mamba env create -f envs/scorpios.yaml
  ```

### Updating SCORPiOs conda environment

- As of SCORPiOs v2.0.0, the conda environment was updated and needs to be reinstalled for users who have a previous version:

  ```
  conda env remove --name scorpios
  mamba env create -f envs/scorpios.yaml
  ```

## Usage

### Setting up your working environment for SCORPiOs

Before any SCORPiOs run, you should:
 - go to SCORPiOs root folder,
 - activate the conda environment with `conda activate scorpios`.

### Running SCORPiOs on example data

Before using SCORPiOs on your data, we recommend running a test with our example data to ensure that installation was successful and to get familiar with the pipeline, inputs and outputs.

#### Example 1: Simple SCORPiOs run

SCORPiOs uses a YAML configuration file to specify inputs and parameters for each run.
An example configuration file is provided: [config_example.yaml](config_example.yaml). This configuration file executes SCORPiOs on toy example data located in [data/example/](data/example/), that you can use as reference for input formats.

The only required snakemake arguments to run SCORPiOs are `--configfile`, the `--use-conda` flag abd the `--scheduler=greedy` option. You also need to specify the number of threads via `--cores`. For more advanced options, you can look at the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

To run SCORPiOs on example data:

```
snakemake --configfile config_example.yaml --use-conda --cores 4 --scheduler=greedy
```

The following output should be generated:
`SCORPiOs_example/SCORPiOs_output_0.nhx`.

#### Example 2: Iterative SCORPiOs run

SCORPiOs can run in iterative mode: SCORPiOs improves the gene trees a first time, and then uses the corrected set of gene trees again as input for a new correction run, until convergence. Correcting gene trees improves orthologies accuracy, which in turn makes synteny conservation patterns more informative, improving the gene tree reconstructions after successive runs. Usually, a small number of iterations (2-3) suffice to reach convergence.

To run SCORPiOs in iterative mode on example data, execute the wrapper bash script `iterate_scorpios.sh`:

```
bash iterate_scorpios.sh --snake_args="--configfile config_example.yaml --cores 4 --scheduler=greedy"
```

Command-line arguments:

```
Required
--snake_args="SNAKEMAKE ARGUMENTS", should at minimum contain --configfile

Optional
--max_iter=MAXITER, maximum number of iterations to run, default=5.
--min_corr=MINCORR, minimum number of corrected sub-trees to continue to the next iteration, default=1.
--starting_iter=ITER, starting iteration, to resume a run at a given iteration, default=1.
```
The following output should be generated: `SCORPiOs_example/SCORPiOs_output_2_with_tags.nhx`.

### Running SCORPiOs on your data

#### Preparing your configuration file
To run SCORPiOs on your data, you have to create a new configuration file for your SCORPiOs run. You will need to format your input data adequately and write your configuration file, using the provided example [config_example.yaml](config_example.yaml) as a guide.

- Copy the example config file `cp config_example.yaml config.yaml`
- Open and edit `config.yaml` to specify paths, files and parameters for your data

To check your configuration, you can execute a dry-run with `-n`.
```
snakemake --configfile config.yaml --use-conda -n
```

#### Running SCORPiOs
Finally, you can run SCORPiOs as described above:

```
snakemake --configfile config.yaml --use-conda --cores 4 --scheduler=greedy
```

or in iterative mode:

```
bash iterate_scorpios.sh --snake_args="--configfile config.yaml --cores 4 --scheduler=greedy"
```

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

- [RAxML](https://github.com/stamatak/standard-RAxML): Stamatakis (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30, 1312–1313.

- [PhyML](http://www.atgc-montpellier.fr/phyml/): Guindon et al. (2010) New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0. Syst Biol, 59, 307–321.

- [TreeBeST](https://github.com/Ensembl/treebest): Vilella et al. (2009) EnsemblCompara GeneTrees: Complete, duplication-aware phylogenetic trees in vertebrates. Genome Res., 19, 327–335.

- [CONSEL](https://github.com/shimo-lab/consel): Shimodaira and Hasegawa (2001) CONSEL: for assessing the confidence of phylogenetic tree selection. Bioinformatics, 17, 1246–1247.
