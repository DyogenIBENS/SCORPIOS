Corrected gene trees
====================

All outputs from SCORPiOs are stored in a folder named :code:`SCORPiOs_jobname/` (with the jobname as specified in the configuration file).

The main output is the **SCORPiOs-optimized gene trees**. Gene trees are provided as a single file in NHX format. We explain in the next section how to visualize SCORPiOs corrections.

The commands above generate:

    - :code:`SCORPiOs_example/SCORPiOs_output_0.nhx` for the simple run,
    - :code:`SCORPiOs_example/SCORPiOs_output_2_with_tags.nhx` for the iterative run.

Outputs are suffixed with a digit representing the iteration number. This number is set to 0 in simple mode and starts at 1 in iterative mode.

..  note::
	NHX tags are added to the corrected gene trees:

	- leaves of corrected subtrees are tagged with :code:`CORR_ID_WGD`,

 	- leaves rearranged to reinsert subtrees are tagged with :code:`MOVED_ID_WGD`;

 	where WGD will be replaced by the name of the corrected WGD and members of the same corrected subtree will be given the same ID.

 	For instance, two leaves, with the tag :code:`CORR_ID_Clupeocephala=1` for one and :code:`CORR_ID_Clupeocephala=2` for the other, belong to two different subtrees corrected for the Clupeocephala WGD.

 	In iterative mode, correction tags are additionaly suffixed with the iteration number of the first correction (i.e :code:`CORR_ID_Clupeocephala_1=1`).

..  warning::

	Please note that :code:`MOVED_ID_` is currently **not** supported in the final output in iterative mode. However, they **are** stored in individual, by iteration, corrected trees (see the entry :code:`save_tmp_trees` of the configuration file), in the exact same format as in the non-iterative run. Corrections at each iteration **can** be visualized with SCORPiOs custom tree visualization script, using all available drawing options (see :ref:`Tree visualization`).


You will also see that some intermediary outputs are stored in different sub-folders of :code:`SCORPiOs_jobname/`. Please see the section on intermediary outputs for a detailed description.

In addition, SCORPiOs writes statistics on key steps of the workflow to the standard output. Thus, to separate output statistics from snakemake logs, you can run:

.. prompt:: bash

	snakemake --configfile config_example.yaml --use-conda >out 2>err

or

.. prompt:: bash

	bash iterate_scorpios.sh --j=example --snake_args="--configfile config_example.yaml" >out 2>err




