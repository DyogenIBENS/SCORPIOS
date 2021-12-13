LORelEi usage instructions
==========================

SCORPiOs LORelEi is implemented as an extension to the main SCORPiOs pipeline, from which it depends as a snakemake subworkflow. LORelEi can be run after a previously completed SCORPiOs job or be directly invoked, in which case the main SCORPiOs job will be automatically run first (except for more advanced usage, see the :ref:`Example 3 <Example 3: SCORPiOs iterative and LORelEi>` below).

In all cases, you will need to prepare two configuration files: one for the main SCORPiOs job and a second for LORelEi. Note, however, that the configuration file for LORelEi is a lot simpler than the one for SCORPiOs, and that it requires no additional input data.

.. warning::
	LORelEi was introduced in SCORPiOs v2.0.0. The implementation of subworkflows required a more recent version of Snakemake, which was updated from version 5.5.4 to 6.6.1 in `SCORPiOs conda environment <https://github.com/DyogenIBENS/SCORPIOS/blob/master/envs/scorpios.yaml>`_. Thus, users having a previous version of the environment need to update it (as explained :ref:`here <Updating SCORPiOs conda environment>`). The update of Snakemake also implies that additional command-line arguments are required to run SCORPiOs, see the :ref:`updated usage instructions<Usage instructions>`.

Running SCORPiOs LORelEi on example data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An example configuration file for SCORPiOs LORelEi is provided: `config_example2_lorelei.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example2_lorelei.yaml>`_. This configuration file executes SCORPiOs LORelEi on toy example `data <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example2/>`_. We explain how to format your own configuration file in the next chapter (see :ref:`LORelEi configuration file`).


As introduced in the previous chapter, two modes are available to run SCORPiOs LORelEi:
:code:`diagnostic` (:ref:`Example 1 <Example 1: SCORPiOs LORelEi diagnostic mode>`) and :code:`likelihood_tests` (:ref:`Example 2 <Example 2: SCORPiOs LORelEi likelihood_tests mode>`). We also show how to use SCORPiOs iterative correction in conjonction with LORelEi (:ref:`Example 3 <Example 3: SCORPiOs iterative and LORelEi>`).


.. important::
    Remember that you should go to the SCORPiOs root folder and activate the conda environment with the command :code:`conda activate scorpios` before running SCORPiOs and/or LORelEi.


Example 1: SCORPiOs LORelEi diagnostic mode
--------------------------------------------

To run SCORPiOs LORelEi on example data in :code:`diagnostic` mode:

.. prompt:: bash

	 snakemake -s scorpios_lorelei.smk --configfile config_example2_lorelei.yaml --use-conda --cores 4 --scheduler=greedy

Outputs for the main SCORPiOs job should be generated in :code:`SCORPiOs_example2/`, while LORelEi results are stored in :code:`SCORPiOs_LORelEi_example2/`.

The following LORelEi output figures should be generated: :code:`SCORPiOs-LORelEi_example2/diagnostic/seq_synteny_conflicts_by_homeologs.svg` and :code:`SCORPiOs-LORelEi_example2/diagnostic/seq_synteny_conflicts_on_genome.svg`.


Example 2: SCORPiOs LORelEi likelihood_tests mode
--------------------------------------------------

If you ran the example in :code:`diagnostic` mode previously, the SCORPiOs main job will not be re-run and LORelEi will re-use pre-computed outputs from the :code:`SCORPiOs_example2/` folder.

The example configuration file contains all necessary arguments to also run LORelEi in :code:`likelihood_tests` mode, you only need to update the :code:`mode` parameter, either in the configuration file or directly in the command-line as follows:

.. prompt:: bash

	 snakemake -s scorpios_lorelei.smk --configfile config_example2_lorelei.yaml --config mode=likelihood_tests --use-conda --cores 4 --scheduler=greedy

Since the configuration is very simple here, you can even omit the configuration file and only provide the three only required arguments in the command-line:

.. prompt:: bash

     snakemake -s scorpios_lorelei.smk --config scorpios_config=config_example2.yaml mode=likelihood_tests dup_genome=Oryzias.latipes --use-conda --cores 4 --scheduler=greedy

The following LORelEi outputs should be generated: :code:`SCORPiOs-LORelEi_example2/lktests/lore_aore_on_genome.svg` (figure) and :code:`SCORPiOs-LORelEi_example2/lktests/lore_aore_summary.tsv` (summary of LORe and AORe gene families).

Example 3: SCORPiOs iterative and LORelEi
------------------------------------------

To run LORelEi in conjonction with SCORPiOs iterative gene tree correction, you will need to run SCORPiOs iterative correction first and then LORelei, specifying the iteration you want to analyze sequence-synteny conflicts on. We recommend using iteration 1 (or 2) of an iterative run for LORelEi, since the number of gene trees considered for correction by SCORPiOs - and thus by LORelEi afterwards - typically decreases a lot in later iterations.

.. prompt:: bash
     
     bash iterate_scorpios.sh --snake_args="--configfile config_example2.yaml --cores 4 --scheduler=greedy"
     snakemake -s scorpios_lorelei.smk --configfile config_example2_lorelei.yaml --config iter=1 --use-conda --cores 4 --scheduler=greedy

The following LORelEi outputs should be generated: :code:`SCORPiOs-LORelEi_example2/diagnostic/seq_synteny_conflicts_by_homeologs.svg` and :code:`SCORPiOs-LORelEi_example2/diagnostic/seq_synteny_conflicts_on_genome.svg`. You can change the :code:`jname` parameter to not overwrite previous results (see :ref:`LORelEi configuration file`).

Running SCORPiOs LORelEi on your data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Like for SCORPiOs, you have to create a new configuration file to run LORelEi on your own data. You can use the example configuration file as a guide to write your own (see :ref:`LORelEi configuration file`) and then run:

.. prompt:: bash

	 snakemake -s scorpios_lorelei.smk --configfile config_lorelei.yaml --use-conda --cores 4 --scheduler=greedy