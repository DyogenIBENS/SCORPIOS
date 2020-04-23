
# Change Log

All notable changes to SCORPiOs workflow, added after the first released version v1.0.0, will be documented in this file.
  
## [1.1.0] - UNRELEASED
 
Workflow updates to improve computational efficiency and scalability to large datasets.
 
### Added
- Option to perform community detection with spectral clustering instead of Girvan-Newman, for improved computational effciency on large datasets.

    See [config_example.yaml](config_example.yaml) for usage instructions.

    The master conda environnment has been updated and now includes scikit-learn to perform spectral clustering. Users with a previous version of scorpios conda env should either install the updated env [envs/scorpios.yaml](envs/scorpios.yaml) or install scikit-learn in their previous version.
 
### Changed
- treebest distmat now replaces fastdist to build input distance matrix for profileNJ
 
### Fixed
- Fixed elevated memory usage in the subtree regrafting step
- Various python code speed-ups