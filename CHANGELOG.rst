All notable changes to SCORPiOs, added after the first released version v1.0.0, will be documented here.

UNRELEASED
-----------

If I release it, SCORPiOs will need to be run with the --scheduler=greedy snakemake args (default scheduler in 6.6.1 is unstable --> the workflow gets stuck)...

Added
^^^^^
- New option to recompute branch-lengths with RAxML after all subtrees corrections (instead of treebest phyml).
- Removed the deprecated buffer_size argument from the configuration file.

Changed
^^^^^^^
- Updated SCORPiOs conda environment to snakemake version 6.6.1


[Version 1.3.0] - 23/11/2020 - `v1.3.0 <https://github.com/DyogenIBENS/SCORPIOS/tree/v1.3.0>`_
-------------------------------------------

Further workflow updates to improve computational efficiency and simplify its usage. Note that the main conda environnment has been updated and needs to be reinstalled if you have the previous version.

Added
^^^^^
- Default output trees in iterative mode have now tags for corrected internal WGD nodes (in addition to corresponding descending leaves).

Changed
^^^^^^^
- **RAxML** replaces PhyML to compute site likelihood (CONSEL input for likelihood AU-tests).
- Simplified usage for iterative correction: the wrapper script directly parses the YAML configuration file so that providing the jobname is no longer required (see the updated :ref:`usage instructions <Example 2: Iterative SCORPiOs run>`).


Version 1.2.0 - 18/10/2020 - `v1.2.0 <https://github.com/DyogenIBENS/SCORPIOS/tree/v1.2.0>`_
-------------------------------------------
 
Minor workflow updates to improve computational efficiency and scalability to large datasets.
 
Added
^^^^^
- Option to perform community detection with spectral clustering instead of Girvan-Newman, for improved computational effciency on large datasets. (See the :ref:`documentation<Spectral clustering>` for usage instructions).

Changed
^^^^^^^
- treebest distmat now replaces fastdist to build input distance matrices for profileNJ.
- :ref:`Installation instructions <Installing SCORPiOs>` now recommend mamba for a faster dependency solving process.

Version 1.1.0 - 19/05/2020 - `v1.1.0 <https://github.com/DyogenIBENS/SCORPIOS/tree/v1.1.0>`_
-------------------------------------

Developed a new visualization tool to inspect tree corrections.

Added
^^^^^
- New tool to generate .png images for corrected trees.
- New html documentation.

Fixed
^^^^^
- Fixed high memory usage in subtrees reinsertion step.
- Minor python code speed-ups.
