Usage instructions
==================

.. important::
	Before using SCORPiOs, you should go to the SCORPiOs root folder and activate the conda environment with the command :code:`conda activate scorpios`.

Running SCORPiOs on example data
--------------------------------

We recommend running a test with our example data to ensure that installation was successful and to get familiar with the pipeline, inputs and outputs.

SCORPiOs uses a YAML configuration file to specify inputs and parameters for each run. An example configuration file is provided: `config_example.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example.yaml>`_. This configuration file executes SCORPiOs on toy example `data <https://github.com/DyogenIBENS/SCORPIOS/blob/master/data/example/>`_, that you can use as reference for input formats. We explain how to format your own configuration file and input files in more details in the next chapter (see :ref:`Data file formats` and :ref:`Configuration file`).

Here, we present the main commands to run SCORPiOs.

Example 1: Simple SCORPiOs run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
The only required Snakemake arguments to run SCORPiOs are :code:`--configfile` and the :code:`--use-conda` flag. Optionally, you can specify the number of threads via the :code:`--cores` option. For more advanced options, you can look at the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_.

To run SCORPiOs on example data, go to the SCORPiOs root folder and run:

.. prompt:: bash

	 snakemake --configfile config_example.yaml --use-conda --cores 4

The following output should be generated: :code:`SCORPiOs_example/SCORPiOs_output_0.nhx`.

Example 2: Iterative SCORPiOs run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SCORPiOs can run in iterative mode: SCORPiOs improves the gene trees a first time, and then uses the corrected set of gene trees again as input for a new correction run, until convergence. Correcting gene trees improves orthologies accuracy, which in turn makes synteny conservation patterns more informative, improving the gene tree reconstructions after successive runs. Usually, a small number of iterations (2-3) suffice to reach convergence.

To run SCORPiOs in iterative mode on example data, execute the wrapper bash script :code:`iterate_scorpios.sh` as follows:

.. prompt:: bash

	 bash iterate_scorpios.sh --snake_args="--configfile config_example.yaml"


The following output should be generated: :code:`SCORPiOs_example/SCORPiOs_output_2_with_tags.nhx`.

Command-line arguments for :code:`iterate_scorpios.sh`
""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Required:**

--snake_args=snakemake_arguments  Snakemake arguments, should at minimum contain :code:`--configfile`.

**Optional:**

--max_iter=maxiter  Maximum number of iterations to run (default=5).

--min_corr=mincorr  Minimum number of corrected subtrees to continue to the next iteration (default=1).

--starting_iter=iter  Starting iteration, to resume a run at a given iteration (default=1).


Running SCORPiOs on your data
-----------------------------

To run SCORPiOs on your data, you have to create a new configuration file for your SCORPiOs run. You will need to format your input data adequately and write your configuration file, using the provided example `config_example.yaml <https://github.com/DyogenIBENS/SCORPIOS/blob/master/config_example.yaml>`_ as a guide.

* Copy the example config file :code:`cp config_example.yaml config.yaml`
* Open and edit :code:`config.yaml` to specify paths, files and parameters for your data

To check your configuration, you can execute a dry-run with :code:`-n`.

.. prompt:: bash

	 snakemake --configfile config.yaml --use-conda -n

Finally, you can run SCORPiOs as described above:

.. prompt:: bash

	 snakemake --configfile config.yaml --use-conda

or in iterative mode:

.. prompt:: bash

	 bash iterate_scorpios.sh --snake_args="--configfile config.yaml"
