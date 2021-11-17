LORelEi configuration file
===========================

As an example, we provide `config_example2_lorelei.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example2_lorelei.yaml>`_, which allows to run SCORPiOs LORelEi on toy example `data <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example2/>`_. The input data are exactly the same input data as required for SCORPiOs, with the same formatting requirements. The LORelEi configuration file only specify additional parameters specific to the LORe analyses.

You will only need to slightly modify the example file in order to run SCORPiOs LORelEi on your data. Note that you still need to write a configuration file for the main SCORPiOs workflow, as detailed in SCORPiOs :ref:`Configuration file` section.

.. code:: yaml

    #======================== SCORPiOs LORelEi EXAMPLE CONFIGURATION FILE =========================#
    ## For a full description of supported settings, please take a look at SCORPiOs documentation.##

    #------------------------------------ REQUIRED PARAMETERS -------------------------------------#
    #                                          (all modes)                                         #

    # Configuration file for SCORPiOs gene tree correction run.
    scorpios_config: "config_example2.yaml"

    # SCORPiOs LORelEi mode, can be 'diagnostic', 'likelihood_tests'.
    mode: "diagnostic"

    # Duplicated genome to plot conflicted gene families or LORe/AORe families on.
    # The corresponding gene coordinates file will be found using the path given in scorpios config.
    dup_genome: "Oryzias.latipes"


    #------------------------------------ OPTIONAL PARAMETERS -------------------------------------#
    #                                          (all modes)                                         #

    # Optional, SCORPiOs iteration to use (0 for a simple run, 1 recommended in iterative SCORPiOs).
    # Default is 0, so you can omit it in normal mode.
    #iter: 0

    # Optional, WGD to consider, you can omit it if you ran SCORPiOs on this WGD only.
    #lore_wgd: "Osteoglossocephalai"

    # Optional, to append a lorelei jobname to the scorpios jobname,
    # (useful to run different lorelei jobs on the same SCORPiOs main job).
    #jname: "myrun1"


    #----------------------------------- Diagnostic PARAMETERS ------------------------------------#

    # Optional, outgroup OR ancestral reconstruction to use as a proxy to the pre-WGD karyotype.
    # You can omit it if you used only 1 outgroup for SCORPiOs correction and want to use it here.
    pre_dup_proxy:
      # use_outgr: "Lepisosteus.oculatus"
      use_anc:  "data/example2/preduplication_ancgenes.tsv"

    #------------------------------------ LK-TESTS PARAMETERS -------------------------------------#

    #Optional, set raxml random seed (for reproducibility)
    #raxml_seed: 1234

    # Optional, to prune gene trees and retain only a subset of the duplicated species
    # sp_to_keep: "Arapaima.gigas Danio.rerio Scleropages.formosus Oryzias.latipes"


.. tip:: Since SCORPiOs LORelei requires few arguments, you may want to omit the configuration file and provide arguments in the command line directly (see the :ref:`Example 2 <Example 2: SCORPiOs LORelEi likelihood-tests mode>` in LORelEi usage instructions). The benefit of using a configuration file is that it keeps a trace of the configuration for reproducibility purposes.

We detail all available settings below.


Supported settings
-------------------

Required parameters (all modes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SCORPiOs configuration file
""""""""""""""""""""""""""""
**Required.** The configuration file for the main scorpios job, see the example `config_example2.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example2.yaml>`_. Please refer to the :ref:`Configuration file` section for a complete description.

Example:

.. code:: yaml

    scorpios_config: "config_example2.yaml"

.. important::
    In order for LORelEi to run, save_subtrees_lktest should always be set to 'y' in the main SCORPiOs configuration file. In addition, when using multiple outgroup genomes to correct gene trees for the WGD of interest, these outgroup species must be monophyletic for SCORPiOs LORelEi to run correctly.

LORelEi mode
"""""""""""""
**Required.** Should be either 'diagnostic' or 'likelihood_tests'.

Example:

.. code:: yaml

    mode: "likelihood_tests"

Duplicated genome for plot
"""""""""""""""""""""""""""
**Required.** Name of the modern genome to plot the localization of conflicted gene families or LORe/AORe families on.

It should be the same name as in the corresponding genes coordinates file. For instance, here, giving Oryzias.latipes will tell LORelEi to fetch :code:`data/example2/genes/genes.Oryzias.latipes.list` (the path to genes coordinates files is extracted from SCORPiOs configuration file, here `config_example2.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example2.yaml>`_).

Chromosomes names should be distinguishable from scaffolds in the BED file. LOReLei identifies genomic regions as chromosomes if their name is in one of the following formats: 'chr1', '1', 'LG1', 'HiC_scaffold_1' or 'group1', numbers can also be roman numerals (chrI, chrII or groupI, groupII etc).

Example:

.. code:: yaml

    dup_genome: "Oryzias.latipes"

Optional parameters (all modes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SCORPiOs iteration
"""""""""""""""""""
**Optional (default=0).** SCORPiOs iteration to use for the LORe analyses (0 for a simple SCORPiOs run, 1 recommended in iterative SCORPiOs).

Example:

.. code:: yaml

    iter: 1

WGD to run LORelEi on
""""""""""""""""""""""
**Optional.** By default, if you ran SCORPiOs to correct gene trees for a single WGD event, LORelEi will run on this same WGD. LORelEi will throw an error if you run SCORPiOs to correct gene trees for several WGDs and do not provide this parameter in the LORelEi config.

Example:

.. code:: yaml

    lore_wgd: "Osteoglossocephalai"

LORelEi jobname
""""""""""""""""
**Optional.** If provided, will append the jobname suffix to the LORelEi output folder. For instance, in the example config `config_example2_lorelei.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example2_lorelei.yaml>`_, adding :code:`jname: "myrun1"` will change the output folder from :code:`SCORPiOs_LORelEi_example2/` to :code:`SCORPiOs_LORelEi_example2_myrun1/`. This is useful if you want to re-run a LORelEi job in the same mode and on the same SCORPiOs data, but with different LORelEi parameters.

Example:

.. code:: yaml

    jname: "myrun1"

Parameters for the diagnostic mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Proxy for the ancestral pre-WGD karyotype
""""""""""""""""""""""""""""""""""""""""""
**Optional.** By default, LORelEi will attempt to use the outgroup genome used for SCORPiOs correction as a proxy for the ancestral karyotype. However, LORelEi will throw an error if you used multiple outgroups for SCORPiOs correction: in this case you need to specify which outgroup to use.

If you are using an outgroup genome, chromosomes names should be distinguishable from scaffolds in the BED file. LOReLei identifies genomic regions as chromosomes if their name is in one of the following formats: 'chr1', '1', 'LG1', 'HiC_scaffold_1' or 'group1', numbers can also be roman numerals (chrI, chrII or groupI, groupII etc).

Example (using an outgroup):

.. code:: yaml

    pre_dup_proxy:
        use_outgr: "Lepisosteus.oculatus"

Alternatively, you can provide an inferred ancestral karyotype, as a 3-columns tab-delimited file: ancestral gene unique identifier (in any format), descending modern genes (space-separated), inferred ancestral post-duplication chromosome (see the `example <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example2/preduplication_ancgenes.tsv>`_).

Example (using ancestral reconstruction):

.. code:: yaml

    pre_dup_proxy:
        use_anc:  "data/example2/preduplication_ancgenes.tsv"

Parameters for the likelihood-test mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RAxML random number seed
""""""""""""""""""""""""""
**Optional (default=1234).** Random seed to initiate RAxML tree search. Different re-runs with the same random seed and gene families alignments will always generate the same resulting trees. We recommend to perform re-runs with different random seed to assess the stability of the results.

Example:

.. code:: yaml

    raxml_seed: 41

Restricted list of duplicated species
""""""""""""""""""""""""""""""""""""""
**Optional.** Space-separated list of duplicated genomes to retain in the gene trees for the likelihood-tests. This can be useful if you have many genomes and want to make a first try with a smaller set. You should keep at least one species from each of the two main clades after the first speciation event.

Example:

.. code:: yaml

    sp_to_keep: "Arapaima.gigas Danio.rerio Scleropages.formosus Oryzias.latipes"