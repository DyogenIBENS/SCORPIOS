Building a dataset
===================

If you do not have gene alignments available for your study species, you will need to start by building an input dataset for SCORPiOs. This implies that you will need to:

1. List the duplicated study species that you want to include in your dataset, while taking care to also select adequate non-duplicated outgroups. You can then generate the corresponding species phylogeny.

.. |br| raw:: html

    <br />

2. Extract all genes along with their CDS nucleotide sequences, for all your study species.

.. |br| raw:: html

    <br />

3. Group genes of all species in families of homologous genes, i.e genes that descend from a common ancestral gene copy.

.. |br| raw:: html

    <br />

4. Build multiple sequence alignments for each gene family.

.. |br| raw:: html

    <br />

5. Optionally, build a gene tree for each family. Alternatively, SCORPiOs can build starting trees for you. In the latter case, you would need to write a file giving the gene to species correspondence.


We've listed above the main steps required to build an input dataset for SCORPiOs. However, assembling genome-wide phylogenetic datasets is a complex process, frequently refined by all leading comparative genomic databases. As a state-of-the-art reference, we recommend the `GeneSeqToFamily paper <https://academic.oup.com/gigascience/article/7/3/giy005/4841850>`_, which describes how Ensembl Compara builds gene families, alignment and gene trees for a given set of species. In addition, the authors developped a Galaxy workflow `GeneSeqToFamily <https://github.com/TGAC/earlham-galaxytools/tree/master/workflows/GeneSeqToFamily>`_, to interactively run each step of the pipeline.

Reference
---------

For a tutorial on how to assemble an input dataset for SCORPiOs:

- `GeneSeqToFamily <https://github.com/TGAC/earlham-galaxytools/tree/master/workflows/GeneSeqToFamily>`_: Thanki et al. (2018) GeneSeqToFamily: a Galaxy workflow to find gene families based on the Ensembl Compara GeneTrees pipeline. GigaScience 7.