Installation
============

SCORPiOs is implemented as a `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ pipeline. Snakemake is a python-based language to build scalable and reproducible workflows. We take advantage of Snakemake's integration with the package manager `Conda <https://docs.conda.io/en/latest/>`_ to ship all SCORPiOs dependencies. The following instructions will help you get a running copy of the pipeline and set up your environnement.

Installing conda
----------------

The Conda package management system manages all SCORPiOs dependencies, including python packages and other software.

To install Conda:

* Download Miniconda3 installer for your system `here <https://docs.conda.io/en/latest/miniconda.html>`_

.. |br| raw:: html

    <br />

* Run the installation script: :code:`bash Miniconda3-latest-Linux-x86_64.sh` or :code:`bash Miniconda3-latest-MacOSX-x86_64.sh`, and accept the defaults

.. |br| raw:: html

    <br />

* Open a new terminal, run :code:`conda update conda` and press :code:`y` to confirm updates


Installing SCORPiOs
-------------------

* Clone the repository:

.. prompt:: bash

	git clone https://github.com/DyogenIBENS/SCORPIOS.git

* Go to SCORPiOs root folder:

.. prompt:: bash

	cd SCORPIOS

* Create the main conda environment. We recommend using `Mamba <https://quantstack.net/mamba.html>`_ for a faster installation:

.. prompt:: bash

	conda install -c conda-forge mamba
	mamba env create -f envs/scorpios.yaml

* **Alternatively,** you can use conda directly (solving dependencies may take a while, ~ up to an hour):

.. prompt:: bash

	conda env create -f envs/scorpios.yaml


.. note:: Once the conda environnment is successfully created, the installation process is complete. You can proceed to the next section and test your installation on example data. Before running SCORPiOs, remember to activate the conda environment with :code:`conda activate scorpios`.

Reference
----------

- `Snakemake: <https://snakemake.readthedocs.io/en/stable/>`_ Köster and Rahmann (2012) Snakemake - A scalable bioinformatics workflow engine. Bioinformatics, 28, 2520–2522.