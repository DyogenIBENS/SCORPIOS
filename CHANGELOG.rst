All notable changes to SCORPiOs, added after the first released version v1.0.0, will be documented here.

[Version 1.1.1] - UNRELEASED
-----------------------------
 
Workflow updates to improve computational efficiency and scalability to large datasets.
 
Added
^^^^^
- Option to perform community detection with spectral clustering instead of Girvan-Newman, for improved computational effciency on large datasets.

    See [config_example.yaml](config_example.yaml) for usage instructions.

    The master conda environnment has been updated and now includes scikit-learn to perform spectral clustering. Users with a previous version of scorpios conda env should either install the updated env [envs/scorpios.yaml](envs/scorpios.yaml) or install scikit-learn in their previous version.
 
Changed
^^^^^^^
- treebest distmat now replaces fastdist to build input distance matrix for profileNJ

Version 1.1.0 - 19/05/2020 - `v1.1.0 <https://github.com/DyogenIBENS/SCORPIOS/tree/v1.1.0>`_
-------------------------------------

Developped a new visualization tool to inspect tree corrections.

Added
^^^^^
- New tool to generate .png images for corrected trees.
- New html documentation.

Fixed
^^^^^
- Fixed high memory usage in subtrees reinsertion step.
- Minor python code speed-ups.