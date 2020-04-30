
# Change Log

All notable changes to SCORPiOs workflow, added after the first released version v1.0.0, will be documented in this file.
  
## [1.1.0] - UNRELEASED

Workflow updates to improve computational efficiency and scalability to large datasets. More efficient tools and algorithms have been selected for the most expensive tasks. Changes concerns mainly the implementation of the steps and do not significantly affect the results.
 
### Added
- Option to perform community detection with spectral clustering instead of Girvan-Newman, for improved computational effciency on large datasets.

    See [config_example.yaml](config_example.yaml) for usage instructions.
 
### Changed
- treebest distmat now replaces fastdist to build input distance matrix for profileNJ
- raxml now replaces phyml to compute site-likelihoods (CONSEL input for likelihood AU-tests)
- improved conda environnment definition: new specific conda environments for groups of rules and simplification of the master scorpios env.
 
### Fixed
- Fixed elevated memory usage in the subtree regrafting step
- Various python code speed-ups