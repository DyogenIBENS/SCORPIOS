All notable changes to SCORPiOs, added after the first released version v1.0.0, will be documented here.

[Version 1.3.0] - UNRELEASED - LINK
-------------------------------------------

Further workflow updates to improve computational efficiency and simplify its usage. Note that the main conda environnment has been updated and needs to be reinstalled.

Added
^^^^^
- Default output trees in iterative mode have now tags for corrected internal WGD nodes (in addition to corresponding descending leaves).

Changed
^^^^^^^
- **RAxML** replaces PhyML to compute site likelihood (CONSEL input for likelihood AU-tests). *TODO branch-length computation in a later release*
- Improved conda envs definition: new specific environments for groups of rules, simplification of the master env. *TODO remove that : I want to keep treebest in the env*
- Simplified usage for iterative correction: the wrapper script directly parses the YAML configuration file so that providing the jobname is no longer required (see the updated :ref:`usage instructions <Example 2: Iterative SCORPiOs run>`).


[Version 1.2.0] - UNRELEASED - LINK
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
