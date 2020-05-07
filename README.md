# SCORPiOs - Synteny-guided CORrection of Paralogies and Orthologies

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3727519.svg)](https://doi.org/10.5281/zenodo.3727519) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.4-brightgreen.svg)](https://snakemake.bitbucket.io)


 SCORPiOs is a **synteny-guided gene tree correction pipeline** for clades that have undergone a whole-genome duplication event. SCORPiOs identifies gene trees where the whole-genome duplication is **missing** or **incorrectly placed**, based on the genomic locations of the duplicated genes across the different species. SCORPiOs then builds an **optimized gene tree** consistent with the known WGD event, the species tree, local synteny context, as well as gene sequence evolution.

 SCORPiOs is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. SCORPiOs takes as input either gene trees or multiple alignments, and outputs the corresponding optimized gene trees.

 For a complete description of SCORPiOs, see our preprint: https://www.biorxiv.org/content/10.1101/2020.01.30.926915v1.full

 ![SCORPiOs illustrated](https://github.com/DyogenIBENS/SCORPIOS/blob/master/doc/scorpios_illustrated.png)

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
      - [Tree visualization](#tree-visualization)
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

The following output should be generated:
`SCORPiOs_example/SCORPiOs_output_0.nhx`.

We explain how to interpret outputs below ([Understanding SCORPiOs outputs](#understanding-scorpios-outputs)).

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
The following output should be generated: `SCORPiOs_example/SCORPiOs_output_2_with_tags.nhx`.

### Running SCORPiOs on your data

#### Data preparation and formatting
SCORPiOs is a flexible gene tree correction pipeline: it can either start from a set of precomputed, phylogeny-reconciled gene trees, or build one from a set of gene multiple aligments using [TreeBeST](https://github.com/Ensembl/treebest). Because SCORPiOs leverages local synteny similarity, i.e evolution of neighboring genes, it requires genome-wide data.

##### Input files
SCORPiOs input files are:
- A single file with a set of phylogeny-reconciled gene trees in NHX format (extended Newick format, see [example](data/example/forest.nhx)) **OR** a genes-to-species mapping file, if working from gene alignments (see [example](data/example/genes_sp_mapping.txt))
- A single file with the corresponding gene multiple alignments in FASTA format (can be compressed with gzip) (see [example](data/example/ali.fa.gz))
- Gene coordinates files for each species in BED format (see [example](data/example/genes/genes.Danio.rerio.bed))
- A species tree in Newick format, with names of ancestral species indicated at internal nodes (see [example](data/example/species_tree.nwk)).

If starting from gene trees, SCORPiOS uses the NHX 'S' (species name) tag to build the gene-species mapping. Otherwise, it uses the gene-to-species mapping file.

More details can be found in [config_example.yaml](config_example.yaml).

##### Parameters
All parameters for a SCORPiOs run have to be indicated in a configuration file, as shown in [config_example.yaml](config_example.yaml).

One critical parameter is the positions of WGD(s) in the species tree and the species to use as outgroup. They both have to be specified together using the `WGDs` keyword. The WGD localization has to be indicated with the name of the last common ancestor of all duplicated species.

For instance, consider the simple species tree below:

*(spotted_gar, (zebrafish, (medaka, (tetraodon, fugu)Tetraodontidae)Euteleosteomorpha)Clupeocephala)Neopterygii;*

![basic tree](https://github.com/DyogenIBENS/SCORPIOS/blob/master/doc/basic_sptree.png)

As the fish WGD is located at the ancestor "Clupeocephala" and we wish to use the spotted_gar as outgroup, the following line shoud be in the configuration file:

```
WGDs:
  Clupeocephala: spotted_gar
```

See the [Complex configurations section](#complex-configurations) for a more sophisticated example.

##### Building an input dataset
If you do not have gene alignments available for your study species, we recommend [this paper](https://academic.oup.com/gigascience/article/7/3/giy005/4841850). The authors explain how Ensembl Compara groups genes in families and subsequently build multiple alignments (and gene trees). In addition, they developped a Galaxy workflow, [GeneSeqToFamily](https://github.com/TGAC/earlham-galaxytools/tree/master/workflows/GeneSeqToFamily), to interactively run each step of the pipeline.

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

The main output is the **SCORPiOs-optimized gene trees**. Gene trees are provided as a single file in NHX format. See the [next section](#tree-visualization) for explanations on how to visualize individual SCORPiOs-corrected trees.

The commands above generate:
- `SCORPiOs_example/SCORPiOs_output_0.nhx` for the simple run
- `SCORPiOs_example/SCORPiOs_output_2_with_tags.nhx` for the iterative run.

Outputs are suffixed with a digit representing the iteration number. This number is set to 0 in simple mode and starts at 1 in iterative mode.

Some intermediary outputs are also stored in different sub-folders (see the [advanced section](#advanced) for a detailed description).

In addition, SCORPiOs writes statistics on key steps of the workflow to the standard output. Thus, to separate output statistics from snakemake logs, you can run:

```
snakemake --configfile config_example.yaml --use-conda >out 2>err
```

or

```
bash iterate_scorpios.sh --j=example --snake_args="--configfile config_example.yaml" >out 2>err
```

#### Tree visualization

SCORPiOs tags corrected nodes in the gene trees to allow easy inspection using tree visualisation softwares (i.e [ETE Toolkit](http://etetoolkit.org/) or [ggtree](https://guangchuangyu.github.io/software/ggtree/)). To facilitate correction inspection, we provide a custom script that generates images for corrected trees.

With the default configuration, SCORPiOs saves individual corrected trees in `SCORPiOs_jobname/Corrections/tmp_whole_trees_0` (or `SCORPiOs_jobname/Corrections/tmp_whole_trees_i` for each iteration i in iterative mode).
Our script `scripts/trees/make_tree_images.py` generates images based on the trees saved in this folder.

For instance, after a simple SCORPiOs run on example data, the following command creates images allowing to view all corrections for the salmonids WGD :

```
python scripts/trees/make_tree_images.py -i SCORPiOs_example/Corrections/tmp_whole_trees_0 -wgd Salmonidae -outgr 'Esox.lucius,Gasterosteus.aculeatus,Oryzias.latipes' -o SCORPiOs_example/Corrections/trees_img
```

Here are example figures for a corrected tree (`SCORPiOs_example/Corrections/trees_img/img_cor_27.png`, right) and its before-correction counterpart (`SCORPiOs_example/Corrections/trees_img/img_ori_27.png`, left):

<img src="https://github.com/DyogenIBENS/SCORPIOS/blob/master/doc/example_ori_27.png" alt="tree1" width="420"/>  <img src="https://github.com/DyogenIBENS/SCORPIOS/blob/master/doc/example_cor_27.png" alt="tree2" width="420"/>

Internal nodes are colored according to convention: duplications in red, dubious duplications in cyan and speciation in blue.  Leaves of the SCORPiOs-corrected subtree, including the wgd subtree and the gene used as outgroup, are shown in the same color in the corrected and uncorrected versions. The corrected WGD node is highlighted with a bigger circle and a grey background.

For a more exhaustive description of the visualization tool please run:
```
python scripts/trees/make_tree_images.py --help
```

Finally, users can inspect original and corrected trees using the [phylo.io](https://phylo.io/) web interface. Through the compare function, original and corrected trees can be inspected side by side, with all differences highlighted.

#### Advanced

Beyond description statistics printed to the standard output, you may want to investigate the detailed results of SCORPiOs for one or several given gene families. This section introduces a few key concepts of SCORPiOs, in order to better understand intermediary outputs.

A gene family in SCORPiOs consists of a non-duplicated outgroup gene and all potential orthologous gene copies in WGD-duplicated species, based on the uncorrected gene trees. For each family, SCORPiOs computes a synteny-derived orthology graph, then a constrained tree topology based on synteny, and finally, if necessary, a synteny-aware corrected tree. Through each of these steps, a gene family is identified by the outgroup gene name.

The orthology relationships between genes are stored in a single file in the `Families/` sub-folder. Raw synteny-predicted orthologies are stored in a single file in `Synteny/`. Predicted orthology groups based on community detection in synteny graphs are stored in a single file in `Graphs/`, along with a summary of the community detection step. Finally, `Corrections/` stores two files, one detailing trees vs synteny consistency and another with the list of successfully corrected trees.

Several tags such as the name of the corrected WGD, the outgroup species and SCORPiOs iteration number are added to each output file, in order to precisely identify outputs in case of complex configurations.

Additional files can be saved if specified in the configuration file, see [config_example.yaml](config_example.yaml) for details.

### Complex configurations

SCORPiOs can correct gene trees that contain more than one whole-genome duplication event. In this case, each WGD is treated independently, starting from the more recent one (closer to the leaves) going up towards the more ancient one (closer to the root). If the WGDs are nested, the subtrees from the more recent events are ignored while correcting for the older WGD event(s), and reinserted after correction using their outgroup as a branching point.

SCORPiOs can also use more than one reference outgroup to correct gene trees. Outgroup(s), separated by commas if more than one, are to be indicated for each WGDs.

For instance, in the example [config_example.yaml](config_example.yaml), WGDs to correct are specified by:

```
WGDs:
  Clupeocephala: 'Lepisosteus.oculatus,Amia.calva'
  Salmonidae: 'Esox.lucius,Gasterosteus.aculeatus,Oryzias.latipes'
```

This specifies that gene trees have to be corrected for the teleost WGD (species below the Clupeocephala ancestor in the species tree) and for the salmonids WGD (species below the Salmonidae ancestor in the species tree). Lepisosteus.oculatus and Amia.calva should be used as outgroups to the teleost WGD and Esox.lucius, Gasterosteus.aculeatus and Oryzias.latipes as outgroups to the salmonids WGD.

Again, we refer to [config_example.yaml](config_example.yaml) for details.

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

For a tutorial on how to assemble an input dataset for SCORPiOs:

- [GeneSeqToFamily](https://github.com/TGAC/earlham-galaxytools/tree/master/workflows/GeneSeqToFamily): Thanki et al. (2018) GeneSeqToFamily: a Galaxy workflow to find gene families based on the Ensembl Compara GeneTrees pipeline. GigaScience 7.