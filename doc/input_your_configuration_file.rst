Configuration file
==================

As alluded to in the previous sections, the configuration for a SCORPiOs run has to be set in a **YAML configuration file**. YAML (standing for YAML Ain't Markup Language) is a data-serialization language with a simple syntax, commonly used for configuration files. 

Using a configuration file is the preferred way to run a Snakemake workflow such as SCORPiOs, as it allows to precisely store input paths and parameters in a static file, thus ensuring **reproducibility**.

As an example, we provide `config_example.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example.yaml>`_, which allows to run SCORPiOs on toy example `data <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/>`_. You will only need to slightly modify this example file in order to run SCORPiOs on your data.

.. warning::
	In this toy dataset, species names are appended to the gene names, for clarity only. This is **not required** and **not** used by SCORPiOs to infer the genes-species mapping.

.. code:: yaml

	#=========================== SCORPiOs EXAMPLE CONFIGURATION FILE ============================#
	##For a full description of supported settings, please take a look at SCORPiOs documentation.#


	#----------------------------------------- INPUTS -------------------------------------------#

	# INPUT1 - The gene trees to correct, as a single file in New Hampshire Extended (.nhx) format.
	trees: data/example/forest.nhx

	# INPUT2 - The multiple sequence alignments used to build the trees, as a single fasta file.
	alis: data/example/ali.fa.gz

	# INPUT3 - The genes coordinates for all duplicated species and outgroup(s).
	#one file per species, in BED (.bed) format.
	genes: data/example/genes/genes.%s.bed

	# INPUT4 - The species tree in Newick format with labelled internal nodes (ancestor names).
	species_tree: data/example/species_tree.nwk

	# Gene-to-species mapping file : a single text file with two columns: gene_name species_name.
	#genes_sp_mapping: data/example/genes_sp_mapping.txt


	#----------------------------------------- OUTPUTS ------------------------------------------#

	# Provide a job name, which will be appended to the output folder name.
	jobname: 'example'


	#---------------------------------------- PARAMETERS ----------------------------------------#

	# Whole-genome duplication(s) and outgroup(s)
	WGDs:
  		Clupeocephala: 'Lepisosteus.oculatus,Amia.calva'
  		Salmonidae: 'Esox.lucius,Gasterosteus.aculeatus,Oryzias.latipes'

	# Whether the average number of syntenic orthologs to include genes in the Orthology Table
	# should be optimized: yes ('y') or no ('n'). Default threshold value if not optimized: 2.0.
	optimize_synteny_support_threshold: 'y'

	# Whether orthologs without synteny support should be discarded from the synteny analysis:
	# yes ('y') or no ('n')
	filter_otable_nosynteny: 'n'

	# Size of the sliding window used in the pairwise synteny analysis:
	windowSize: 15

	# Cut-off on pairwise synteny similarity scores (between 0 and 1):
	cutoff: 0

	# Whether branch-lengths should be recomputed after corrections: yes ('y') or no ('n').
	brlength: 'y'

	# Any species with a poorer assembly quality that should be discarded for synteny analysis.
	## Comment out if you want to use all species in the synteny analysis.##
	lowcov_sp: 'data/example/lowcov'

	# Whether to ignore tree-synteny inconsistencies when an orthogroup graph community contains
	# only a single gene: yes ('y') or no ('n'). These are poorly-supported WGD duplication nodes.
	ignoreSingleGeneCom: 'y'

	# Whether each individual corrected tree and its non-corrected counterpart should be saved,
	# each in a .nhx file: yes ('y') or no ('n'). It facilitates direct inspection of corrections.
	save_tmp_trees: 'y'

	# Whether more detailed intermediary outputs should be saved: yes ('y') or no ('n').
	# Setting `save_subtrees_lktest` to 'y' saves, in addition to default outputs:
	# - constrained tree topologies
	# - profileNJ and TreeBeST synteny-aware trees
	# - AU-tests outputs
	save_subtrees_lktest: 'n'

	# Optionally, use spectral clustering instead of Girvan-Newman for graph community detection.
	spectral: 'n'

	# Iterative-correction related option, automatically updated by the wrapper iterate_scorpios.sh
	## DO NOT MODIFY MANUALLY even if using iterative mode.##
	current_iter: 0


	#---------------------------------------- RESSOURCES ----------------------------------------#

	# Maximum number of threads (will never use more than this number).
	# It will be restricted to the number specified via --cores (1 if --cores is not invoked).
	ncores: 14

	# Memory (--buffer_size) parameter for a bash sort. If decreased, more /tmp space will be used.
	buffer_size: 10G

	# Use a parallelization scheme specific to large jobs: yes ('y') or no ('n').
	parallel_scheme_large_job: 'n'

	# Limit number of cores for the branch length computation (after all corrections).
	## Uncomment to reduce RAM usage.##
	#limit_threads_for_branch_lengths: 37

We detail each of the settings in the next section.

Supported settings
------------------

Input data
^^^^^^^^^^

The gene trees to correct (INPUT1)
"""""""""""""""""""""""""""""""""""
**Optional, can be replaced by alternative_INPUT1.** Trees should be provided as a single file in New Hampshire Extended (.nhx) format. Please refer to the :ref:`Data file formats` section for file format details.

Example:

.. code:: yaml

	trees: data/example/forest.nhx

.. important::
	If you want to build the trees from gene sequence alignments using `TreeBeST <https://github.com/Ensembl/treebest>`_, you should remove or comment out the :code:`tree` entry.


The multiple sequence alignments (INPUT2)
""""""""""""""""""""""""""""""""""""""""""
**Required.** Multiple sequence alignments used to build the trees, as a single file in fasta (.fa) format. The file can be gzipped (.gz) or not. Please refer to the :ref:`Data file formats` section for file format details.
Example:

.. code:: yaml

	alis: data/example/ali.fa.gz

The genes coordinates (INPUT3)
"""""""""""""""""""""""""""""""
**Required.** The genes coordinates for all duplicated species and outgroup(s), one file per species, in BED (.bed) format. Files can be bzipped2 (.bz2). Please refer to the :ref:`Data file formats` section for file format details.

Example:

.. code:: yaml

	genes: data/example/genes/genes.%s.bed

The species tree (INPUT4)
""""""""""""""""""""""""""
**Required.** The species tree in Newick format (.nwk) with labelled internal nodes (ancestor names). Please refer to the :ref:`Data file formats` section for file format details.

Example:

.. code:: yaml

	species_tree: data/example/species_tree.nwk

The gene-to-species mapping (alternative_INPUT1)
"""""""""""""""""""""""""""""""""""""""""""""""""
**Optional, can be replaced by INPUT1.** Gene-to-species mapping file : a single text file with two columns: gene_name; species_name. Please refer to the :ref:`Data file formats` section for file format details.

Example:

.. code:: yaml

	genes_sp_mapping: data/example/genes_sp_mapping.txt

..  important::

	You should use the :code:`genes_sp_mapping` entry **only** if you wish to build starting trees from gene sequence alignments with `TreeBeST <https://github.com/Ensembl/treebest>`_.


Outputs
^^^^^^^
Unique jobname
""""""""""""""
**Required.** A (descriptive) job name, which will be appended to the output folder name. All results will be stored in the output folder :code:`SCORPiOs_jobname/`. This allows to invoke different SCORPiOs runs (e.g with different input data or parameters).

Example:

.. code:: yaml

	jobname: 'example'

..  tip::
	Using the example, the corrected gene trese file will be: :code:`SCORPiOs_example/SCORPiOs_output_0.nhx`.


Parameters
^^^^^^^^^^

Whole-genome duplication(s) and outgroup(s)
"""""""""""""""""""""""""""""""""""""""""""

**Required.** Each WGD event in the species tree should be indicated via the name of the ancestor of all duplicated species. Then, for each WGD, provide one or several outgroup species to use as reference in the synteny analysis. Any non-duplicated species can be used as outgroup, but phylogenetically close outgroup should be preferred as synteny with duplicated species will be more conserved. Multiple reference outgroups can be provided as a comma-separated list. For an illustrated explanation on how to specify the duplicated ancestor, please see the "Data preparation and formatting" section.

Example:

.. code:: yaml

	WGDs:
  		Clupeocephala: 'Lepisosteus.oculatus,Amia.calva'
  		Salmonidae: 'Esox.lucius,Gasterosteus.aculeatus,Oryzias.latipes'


Synteny threshold optimization
""""""""""""""""""""""""""""""
**Optional (default='n').** Whether the minimum required number of syntenic orthologs to include genes as potential orthologs should be optimized: yes ('y') or no ('n'). Default value if the threshold is not optimized is 2.0.


Example:

.. code:: yaml

	optimize_synteny_support_threshold: 'y'

Filter orthologs based on synteny
"""""""""""""""""""""""""""""""""
**Optional (default='n').** Whether phylogenetic orthologs without synteny support should be discarded from the synteny analysis: yes ('y') or no ('n').

Example:

.. code:: yaml

	filter_otable_nosynteny: 'n'


Sliding window size
"""""""""""""""""""
**Optional (default=15).** Size of the sliding window used in the pairwise synteny analysis.

Example:

.. code:: yaml

	windowSize: 15


Cut-off on :math:`{\Delta}S` score
""""""""""""""""""""""""""""""""""
**Optional (default=0).** Cut-off on pairwise synteny similarity scores (float between 0 and 1).

Example:

.. code:: yaml

	cutoff: 0

Branch-lengths computation after correction
""""""""""""""""""""""""""""""""""""""""""""
**Optional (default='y').** Whether branch-lengths should be recomputed after subtree corrections: yes ('y') or no ('n').

Example:

.. code:: yaml

	brlength: 'y'

Lower-quality genome assemblies
"""""""""""""""""""""""""""""""
**Optional.** A file listing species with a poorer assembly quality that should be discarded for synteny analysis. You should still provide their genes coordinate files. 

Example:

.. code:: yaml

	lowcov_sp: 'data/example/lowcov'

..  note::

	Comment out or remove the :code:`lowcov_sp` entry if you want to use all species in the synteny analysis.

Poorly-supported WGD duplication nodes
""""""""""""""""""""""""""""""""""""""
**Optional (default='y').** Whether to ignore tree-synteny inconsistencies when an orthology graph community contains only a single gene: yes ('y') or no ('n'). These are poorly-supported WGD duplication nodes.

Example:

.. code:: yaml

	ignoreSingleGeneCom: 'y'

Save individual correction tree files
"""""""""""""""""""""""""""""""""""""
**Optional (default='n').** Whether each individual corrected tree and its non-corrected counterpart should be saved, each in a .nhx files: yes ('y') or no ('n').

..  tip::

	This facilitates direct inspection of corrections.

Example:

.. code:: yaml

	save_tmp_trees: 'y'

Save additional intermediary outputs
"""""""""""""""""""""""""""""""""""""
**Optional (default='n').** Whether more detailed intermediary outputs should be saved: yes ('y') or no ('n'). Setting :code:`save_subtrees_lktest` to 'y' saves, in addition to default outputs:

	- constrained tree topologies

	- profileNJ and TreeBeST synteny-aware trees

	- AU-tests outputs

.. note::

	A description of all intermediary outputs can be found in the "Outputs description" chapter.

Example:

.. code:: yaml

	save_subtrees_lktest: 'n'

Spectral clustering
"""""""""""""""""""""
**Optional (default='n').** Use spectral clustering instead of Girvan-Newman for graph community detection. On large graphs, spectral clustering is computationally more efficient. Consider using it if your dataset contains many duplicated species.

Example:

.. code:: yaml

	sapectral: 'y'

Iterative-correction related option
"""""""""""""""""""""""""""""""""""
Iterative-correction related option, (automatically updated by the wrapper :code:`iterate_scorpios.sh`). 0 if SCORPiOs is run in simple mode, current iteration otherwise.

..  warning::

	Do not modify manually, even if using iterative mode.

Example:

.. code:: yaml

	current_iter: 0


Computational ressources
^^^^^^^^^^^^^^^^^^^^^^^^

Threads
"""""""
**Required.** Maximum number of threads. SCORPiOs will never, in any case, use more than this number, nor more than the number of threads specified via :code:`--cores` (1 if :code:`--cores` is not invoked). In other words, the number of threads will always be min(:code:`ncores`, :code:`--cores`).

Example:

.. code:: yaml

	ncores: 14

Memory for bash sort
""""""""""""""""""""
**Required.** Memory (:code:`--buffer_size`) parameter for a bash sort. If decreased, more :code:`/tmp` space will be used.

Example:

.. code:: yaml

	buffer_size: 10G

Parallelization scheme
"""""""""""""""""""""""
**Optional (default='n').** Use a parallelization scheme specific to large jobs: yes ('y') or no ('n'). If the number of duplicated species is large (~ >25), the default parallelization scheme is slow in snakemake. Setting :code:`parallel_scheme_large_job` to 'y' will greatly reduce computation time.

Example:

.. code:: yaml

	parallel_scheme_large_job: 'n'

Parallel jobs for branch-length computation (soon deprecated)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
**Optional**. Limit the number of cores for the branch length computation (after all corrections). Recomputing branch lengths can be RAM intensive for large trees (SCORPiOs uses TreeBeST PhyML here). To use less RAM, you may want to reduce the number of parallel jobs.

Example:

.. code:: yaml

	limit_threads_for_branch_lengths: 10
