Data preparation
================


SCORPiOs is a flexible gene tree correction pipeline: it can either start from a set of precomputed, phylogeny-reconciled gene trees, or build one from a set of gene multiple aligments using `TreeBeST <https://github.com/Ensembl/treebest>`_. 

If you do not have gene trees or gene alignments readily available for your study species, please refer to the :ref:`Building a dataset` section.

.. warning::
	Because SCORPiOs leverages local synteny similarity, i.e evolution of neighboring genes, it requires **genome-wide data**.

Input files
-----------

SCORPiOs requires four input files, which are:

1. A set of phylogeny-reconciled gene trees as a single file in NHX format (extended Newick format, see `example <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/forest.nhx>`_ and :ref:`description<Gene tree file>`). 

**OR**

1. (bis) A genes-to-species mapping file, if starting the process from gene alignments only (see `example <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/genes_sp_mapping.txt>`_ and :ref:`description<Genes-to-species mapping file>`).

.. |br| raw:: html

    <br />

2. The gene multiple alignments corresponding to the gene trees, as a single file in FASTA format (can be compressed with gzip, see `example <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/ali.fa.gz>`_ and :ref:`description<Gene multiple alignment file>`).

.. |br| raw:: html

    <br />

3. Gene coordinates files for each species in BED format (see `example <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/genes/genes.Danio.rerio.bed>`_ and :ref:`description<Gene coordinates files>`).

.. |br| raw:: html

    <br />

4. A species tree in NEWICK format, with names of ancestral species indicated at internal nodes (see `example <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/species_tree.nwk>`_ and :ref:`description<Species tree file>`).


For a detailed description of expected formats please refer to the Data formatting section.

.. note::
	If starting from gene trees, SCORPiOS uses the NHX :code:`S` (species name) tag to build the gene-species mapping. Otherwise, it uses the gene-to-species mapping file.


Parameters
----------

All parameters for a SCORPiOs run have to be indicated in the YAML configuration file, as shown in `config_example.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example.yaml>`_.

A critical parameter is the position(s) of WGD(s) in the species tree and the species to use as outgroup(s). They both have to be specified together using the :code:`WGDs` keyword. The WGD position should be indicated with the name of the last common ancestor of all duplicated species.

For instance, consider the species tree below:

::

	(spotted_gar, (zebrafish, (medaka, (tetraodon, fugu)Tetraodontidae)Euteleosteomorpha)Clupeocephala)Neopterygii;

.. image:: https://raw.githubusercontent.com/DyogenIBENS/SCORPIOS/master/doc/img/basic_sptree.png

The fish WGD occurred in the "Clupeocephala" ancestor, and we wish to use the spotted_gar as outgroup. This should be specified in the configuration file as:

.. code-block:: yaml

	WGDs:
  	  Clupeocephala: spotted_gar

For a detailed description of all parameters available in SCORPiOs please refer to the :ref:`Configuration file` section.


Complex configurations
----------------------

SCORPiOs can correct gene trees that contain **more than one whole-genome duplication event**.

In this case, each WGD is treated independently, starting from the more recent one (closer to the leaves) going up towards the more ancient one (closer to the root). If the WGDs are nested, the subtrees from the more recent events are ignored while correcting for the older WGD event(s), and reinserted after correction using their outgroup as a branching point.

SCORPiOs can also use more than one reference outgroup to correct gene trees. Outgroup(s), separated by commas if more than one, are to be indicated for each WGDs.

For instance, in the example `config_example.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example.yaml>`_, WGDs to correct are specified by:

.. code-block:: yaml

	WGDs:
  	  Clupeocephala: Lepisosteus.oculatus,Amia.calva
  	  Salmonidae: Esox.lucius,Gasterosteus.aculeatus,Oryzias.latipes

This specifies that gene trees have to be corrected for the teleost WGD (species below the Clupeocephala ancestor in the `species tree <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/species_tree.nwk>`_) and for the salmonids WGD (species below the Salmonidae ancestor in the `species tree <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/species_tree.nwk>`_). Lepisosteus.oculatus and Amia.calva should be used as outgroups to the teleost WGD and Esox.lucius, Gasterosteus.aculeatus and Oryzias.latipes as outgroups to the salmonid WGD.
