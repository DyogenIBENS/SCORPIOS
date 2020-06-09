All notable changes to SCORPiOs, added after the first released version v1.0.0, will be documented here.

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