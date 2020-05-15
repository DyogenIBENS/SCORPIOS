Intermediary outputs
====================

Beyond description statistics printed to the standard output and final corrected trees, you may want to investigate step-by-step results of SCORPiOs for one or several specific gene families.

.. important::
	A gene family in SCORPiOs consists of a non-duplicated outgroup gene and all potential orthologous gene copies in WGD-duplicated species, based on the uncorrected gene trees. For each family, SCORPiOs computes a synteny-derived orthology graph, then a constrained tree topology based on synteny, and finally, if necessary, a synteny-aware corrected tree. Through each of these steps, a gene family is identified by the outgroup gene name.

.. tip::
	Several suffixes such as the name of the corrected WGD, the outgroup species and SCORPiOs iteration number are added to each output file, in order to precisely identify outputs, even in case of complex configurations.

Comprehensive list of orthologs
--------------------------------
The **orthology relationships** between genes of **duplicated species and outgroup** are stored in a single file, whose name starts with :code:`OrthoTable_`, and located inside the :code:`Families/` sub-folder. This table retains all gene copies since the ingroup/outgroup speciation node, as well as any other homologs with a loosely similar syntenic context.

In the example run, one such file is:

 :code:`SCORPiOs_example/Families/OrthoTable_Salmonidae_Esox.lucius_0`.

The three first columns of the file describe outgroup genes (chromosome, index of the gene on the chromosome, gene name). Other columns gives the predicted orthologous genes in duplicated species, with the same information.

Pairwise synteny orthology predictions
--------------------------------------
Raw **synteny-predicted orthologies** amongst duplicated species are stored in a single file, whose name starts with :code:`Sorted_SyntenyOrthoPred_`, and located in :code:`Synteny/`.

In the example run, one such file is:

 :code:`SCORPiOs_example/Synteny/Sorted_SyntenyOrthoPred_Salmonidae_Esox.lucius_0.gz`.

It is a 4-columns gunzipped (.gz) file, giving orthologous genes predicted between duplicated species, after the pairwise synteny analysis. The first columns shows a gene in a duplicated species 1, the second gives its predicted ortholog in duplicated species 2, the third gives the associated :math:`{\Delta}S` synteny score and the fourth the outgroup gene name. 

.. note::
	In this file, species names are appended to the gene names.

Orthogroups in synteny graphs
------------------------------
Predicted **orthogroups** based on community detection in **synteny graphs** are stored in a single file (:code:`GraphsOrthogroups_`) in :code:`Graphs/`, along with a summary of the community detection step (:code:`Summary_`). 

In the example run, the following file gives predicted orthogroups:

 :code:`SCORPiOs_example/Graphs/GraphsOrthogroups_Clupeocephala_Lepisosteus.oculatus_0`.

The first column gives the name of the outgroup gene with an appended "a" or "b" letter to uniquely identify each the two post-WGD orthogroups. Other columns gives the duplicated species gene members.

In addition, :code:`SCORPiOs_example/Graphs/Summary_Clupeocephala_Lepisosteus.oculatus_0` is a simple 3-columns table describing the community detection step. The outgroup gene is indicated in the first column, followed by the algorithm used for community detection and the number of graph edges removed in the second and third columns, respectively.


Subtree corrections
-------------------

Correction summary
^^^^^^^^^^^^^^^^^^^

The :code:`Corrections/` folder stores two files, one detailing **trees vs synteny consistency** and another with the list of **successfully corrected subtrees**.

In the example run, the following file gives an inconsistency summary (with respect to the Clupeocephala WGD):

 :code:`SCORPiOs_example/Corrections/Trees_summary_Clupeocephala_0`.

In addition, the following file lists corresponding accepted corrections:

 :code:`SCORPiOs_example/Corrections/Accepted_Trees_Clupeocephala_0`.



Subtree corrections (additional)
--------------------------------
Additional files can be saved if specified in the configuration file (see the configuration keyword :code:`save_subtrees_lktest`).


Constrained tree topologies
^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Constrained tree topologies** are stored in the :code:`Trees/ctrees_0/` folder (:code:`Trees/ctrees_i/` for each iteration i in iterative mode).

In the example, one constrained topology file is:

 :code:`SCORPiOs_example/Trees/ctrees_0/Clupeocephala/C_102697250_Lepisosteus.oculatus.nh`

This file gives the constrained tree topology for the gene family identified by the outgroup gene :code:`102697250_Lepisosteus.oculatus`, in the newick format.

profileNJ and TreeBeST solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Synteny-aware trees built with `**ProfileNJ** <https://github.com/maclandrol/profileNJ>`_ and `**TreeBeST phyml** <https://github.com/Ensembl/treebest>`_, using the **constrained tree topology**, are stored in the :code:`Corrections/PolyS_0/` and :code:`Corrections/TreeB_0/` folders, respectively.


In the example, one profileNJ tree file is:

 :code:`SCORPiOs_example/Corrections/PolyS_0/Clupeocephala/102697250_Lepisosteus.oculatus.nh`.

Trees are in the newick format.

..  note::

	SCORPiOs does not build a TreeBeST tree if the profileNJ solution is accepted. In this case, TreeBeST tree files will be empty.

Likelihood au-tests
^^^^^^^^^^^^^^^^^^^^
Output of the likelihood au-tests are stored in the :code:`Corrections/Res_polylk_0/` and :code:`Corrections/Res_treeBlk_0/` folders. These are direct outputs from the `CONSEL <https://github.com/shimo-lab/consel>`_ software.

In the example, the following file gives **au-test likelihood tests** results for the **original subtree** vs the corresponding synteny-aware tree resolved with **profileNJ**:

 :code:`SCORPiOs_example/Corrections/Res_polylk_0/Clupeocephala/Res_102697250_Lepisosteus.oculatus.txt`

Similarly, files in the :code:`SCORPiOs_example/Corrections/Res_polylk_0/Clupeocephala/` stores comparisons of **original subtree vs TreeBeST phyml** solution.

..  note::

	Au-test result files for TreeBeST solutions will be empty if the profileNJ solution was accepted.




