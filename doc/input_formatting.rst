Data file formats
=================

Input files should typically follow the general conventions for their file format. Here is a more detailed overview of what SCORPiOs expects to find within each input file.


Gene tree file
--------------

All gene trees should be listed into a single NHX (New Hampshire Extended) file, separated by '//'. 

The gene trees should be phylogeny-reconciled: internal nodes should be tagged with the :code:`D` attribute to specify if they correspond to a speciation or a duplication (e.g. :code:`D=N` or :code:`D=Y`). Leaves should be tagged using the :code:`S` attribute to indicate the species (e.g. :code:`S=Danio.rerio`). Optionally, the :code:`DD` and :code:`DCS` attributes can help flag dubious duplications at internal nodes (:code:`DD=Y` or :code:`DCS=0` for dubious duplications).

See an example gene tree file `here <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/forest.nhx>`_.

Obtaining phylogeny-reconciled trees
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Phylogenetic reconciliation compares gene trees to the species history, to annotate whether internal nodes in the gene tree correspond to speciation or duplication events.
For many taxa, precomputed phylogeny-reconciled gene trees can be downloaded from public comparative genomics databases, such as Ensembl.

If phylogeny-reconciled gene trees are not available for your study species, we suggest running SCORPiOs from gene alignments. SCORPiOs will then build an initial set of gene trees from the alignments using TreeBeST, and will optimize them using synteny information.

Alternatively, you can use TreeBeST to reconcile NEWICK gene trees computed with another program to the species phylogeny.

.. tip ::

	TreeBeSt is one of SCORPiOs dependencies, therefore you can invoke it without having to re-install it. Simply make sure you have activated scorpios conda environnment before running the command below (activation command: :code:`conda activate scorpios`).

Given a species tree in the newick format and a gene tree in newick format that you want to convert to (.nhx) phylogeny-reconciled format. You will then need to:

* append "_"+speciesname to gene names in your gene tree

* run :code:`treebest sdi` to reconcile the gene tree:

.. prompt:: bash

	treebest sdi my_gene_tree.nwk -s my_species_tree.nwk > my_gene_tree.nhx

..  note::

	Unreconciled tree formats (such as Newick) will raise an error when running SCORPiOs. If you do not have reconciled trees available, we recommend using gene alignments as your primary SCORPiOs input.
	
Converting to NHX format
^^^^^^^^^^^^^^^^^^^^^^^^

If you already have phylogeny-reconciled gene trees, but they are not in NHX format, you will need to convert them.
We recommend the `Phylo module <https://biopython.org/wiki/Phylo>`_ in Biopython, which will handle a variety of tree formats.


Gene multiple alignment file
-----------------------------

The multiple sequence alignments used to build the trees should be provided as a single file in fasta (.fa) format. These should be nucleotide sequence alignments (CDS alignments, back-translated from protein alignments). The file can be gzipped (.gz) or not. 
Alignments should be separated by '//' and should appear in the same order as their corresponding trees. Gene names should match those used in the trees.

This file can also be used as the primary input for SCORPiOs if phylogeny-reconciled trees are not available. SCORPiOs will then use these alignments to build the initial set of trees.

See an example gene multiple alignment file `here <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/ali.fa.gz>`_.

Building a gene sequence multiple alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Precomputed gene sequence alignments can be downloaded from a variety of public sources (databases, publications) for many species sets. However, if you cannot find precomputed gene alignments that suit your analysis, you will need to start the work from scratch and build your own gene alignments before you can use SCORPiOs. We provide a dedicated outline here: :ref:`Building a dataset`.

.. warning::
	Building gene multiple alignments is non-trivial and we recommend you use precomputed alignments if possible, if you are not familiar with this task.

Gene coordinates files
----------------------

All genes and their coordinates should be provided as a separate file per species in BED format. Files can be bzipped2 (.bz2). Files must be a minimal BED format with 4 columns: chromosome; start; end; gene_name. 
All file names should follow similar conventions, to be retrieved with a regular expression. Gene names should match those used in the trees and alignments.

See an example gene coordinates file `here <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/genes/genes.Danio.rerio.bed>`_.

Alternatively, genes coordinate file can be provided in 'dyogen' format: a tab-separated file with 5 columns (chromosome; start; end; strand; gene_name).

See an example gene coordinates file in 'dyogen' format `here <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example2/genes/genes.Danio.rerio.list>`_.

Genes-to-species mapping file
-----------------------------

This is an alternative input to the gene tree file, when phylogeny-reconciled gene trees are unavailable and gene alignments are used as the primary input for SCORPiOs.

Correspondances between gene IDs and species names should be provided as a single text file with two columns: gene_name; species_name. Genes from the same family should appear consecutively in the file. Genes families should be separated by '//'. Families should appear in the same order as their corresponding alignment in the alignments file . Gene names and species names should be the same as in the alignment and species tree, respectively.

See an example genes-to-species mapping file `here <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/genes_sp_mapping.txt>`_.

Species tree file
-----------------

The species tree in NEWICK format, with names of ancestral species indicated at internal nodes. The species tree should contain all species included in the gene trees. Species names should not contain underscores '_'. For optimal tree reconstruction with `TreeBeST <https://github.com/Ensembl/treebest>`_, the tree should not contain polytomies.

See an example species tree file `here <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/species_tree.nwk>`_.