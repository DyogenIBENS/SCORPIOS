All notable changes to SCORPiOs, added after the first released version v1.0.0, will be documented here.

[Version 1.3.0] - UNRELEASED - LINK
-------------------------------------------

Replaced PhyML by RAxML for improved computational efficiency on large datasets.

Changed
^^^^^^^
- **RAxML** replaces PhyML to compute site likelihood (CONSEL input for likelihood AU-tests) *and TODO branch-length computation*
- Improved conda environnment definition: new specific conda environments for groups of rules, simplification of the master scorpios env and installation instructions now recommend mamba.
- Simplified iterative correction: *TODO* the wrapper script now parses the YAML configuration file, see new usage instructions.
- *TODO* New option to tag corrected WGD nodes instead of leaves in the final output file (or both).

[Version 1.2.0] - UNRELEASED - LINK
-------------------------------------------
 
Workflow updates to improve computational efficiency and scalability to large datasets.
 
Added
^^^^^
- Option to perform community detection with spectral clustering instead of Girvan-Newman, for improved computational effciency on large datasets. (See the :ref:`documentation<Spectral clustering>` for usage instructions.)

Changed
^^^^^^^
- treebest distmat now replaces fastdist to build input distance matrix for profileNJ

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
