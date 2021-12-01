.. SCORPiOs documentation master file, created by
   sphinx-quickstart on Mon May 11 13:13:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SCORPiOs documentation!
==================================


.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3727519.svg
    :target: https://doi.org/10.5281/zenodo.3727519

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0

.. image:: https://img.shields.io/badge/snakemake-≥6.6.1-brightgreen.svg
    :target: https://snakemake.bitbucket.io

.. image:: https://readthedocs.org/projects/scorpios/badge/?version=latest
   :target: https://scorpios.readthedocs.io/en/latest/?badge=latest

SCORPiOs is a **synteny-guided gene tree correction pipeline** for clades that have undergone a whole-genome duplication event. SCORPiOs identifies gene trees where the whole-genome duplication is **missing** or **incorrectly placed**, based on the genomic locations of the duplicated genes across the different species. SCORPiOs then builds an **optimized gene tree** consistent with the known WGD event, the species tree, local synteny context, as well as gene sequence evolution.


SCORPiOs is implemented as a `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ pipeline. SCORPiOs takes as input either gene trees or multiple alignments, and outputs the corresponding optimized gene trees.

For more information, you can take a look at `SCORPiOs publication <https://doi.org/10.1093/molbev/msaa149>`_.

.. image:: https://raw.githubusercontent.com/DyogenIBENS/SCORPIOS/master/doc/img/scorpios_illustrated.png


References
----------

SCORPiOs uses the following tools to build and test gene trees:

   * `ProfileNJ <https://github.com/maclandrol/profileNJ>`_: Noutahi et al. (2016) Efficient Gene Tree Correction Guided by Genome Evolution. PLOS ONE, 11, e0159559.


   * `RAxML <https://github.com/stamatak/standard-RAxML>`_: Stamatakis (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30, 1312–1313.


   * `PhyML <http://www.atgc-montpellier.fr/phyml/>`_: Guindon et al. (2010) New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0. Syst Biol, 59, 307–321.


   * `TreeBeST <https://github.com/Ensembl/treebest>`_: Vilella et al. (2009) EnsemblCompara GeneTrees: Complete, duplication-aware phylogenetic trees in vertebrates. Genome Res., 19, 327–335.


   * `CONSEL <https://github.com/shimo-lab/consel>`_: Shimodaira and Hasegawa (2001) CONSEL: for assessing the confidence of phylogenetic tree selection. Bioinformatics, 17, 1246–1247.


.. toctree::
   :caption: Quick start
   :name: quick_start
   :hidden:
   :maxdepth: 1

   getting_started_installation.rst
   getting_started_usage.rst


.. toctree::
   :caption: Input Data Preparation
   :name: input_data_preparation
   :hidden:
   :maxdepth: 1

   input_description.rst
   input_formatting.rst
   input_your_configuration_file.rst
   input_building_a_dataset.rst


.. toctree::
   :caption: Outputs Description
   :name: outputs_description
   :hidden:
   :maxdepth: 1

   output_genetrees.rst
   output_treeviz.rst
   output_advanced.rst


.. toctree::
   :caption: LORelEi (LORe Extension)
   :name: lorelei
   :hidden:
   :maxdepth: 1

   lorelei_introduction.rst
   lorelei_usage.rst
   lorelei_configuration_file.rst

.. toctree::
   :caption: Project Information
   :name: project_information
   :hidden:
   :maxdepth: 1

   scripts.rst
   project_changelog.rst
   project_info.rst
