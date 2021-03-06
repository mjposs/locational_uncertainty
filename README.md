# Locational uncertainty

This repository contains the algorithms used to solve robust optimization problems studied in the paper

*Marin Bougeret, Jérémy Omer, and Michael Poss: Optimization problems in graphs with locational uncertainty. Available on [hal](https://hal.archives-ouvertes.fr/hal-03331166)*

The following algorithms are available:
* The exact branch-and-cut algorithm  described in Section 4 and Algorithm 1 of the manuscript, see function `exact()`
* A heuristic algorithm based on deterministic distances, as described in Section 5 of the manuscript, see function `heuristic_determinisitic()`. The deterministic costs can be based on maximum distances, average distances or distances between the centers of the uncertainty sets.
* A heuristic algorithm based on the two-stage reformulation and affine decision rules approximation proposed by [ZRRH21](https://doi.org/10.1287/ijoc.2020.1025 "Robust optimization for models with uncertain second-order cone and semidefinite programming constraints.") and described in Appendix H of the manuscript see function `heuristic_adr`

## Guide

The repository contains the following Julia files:
* **algo.jl**: Contains all algorithms
* **data.jl**: Contains functions parsing data files and creating instances using generators
* **process_data.jl**: Contains functions parsing result files to create the plots used in the manuscript
* **run.jl**: Read/Create the instances and calls the algorithms
* **X_and_cost.jl**: Contains the functions creating the constraints specific to each application `build_IP_model()`, and the computation of the objective function `c()`
* **dmax_study.jl**: the functions therein were used for the search of approximation ratios, they are not relevant for the actual solution of locational uncertainty problems

The code currently contains two applications: 
* A Steiner Tree problem.
* A strategic facility.
Additional applications can be added by creating the corresponding functions.

To run the code, you should create the result folders (res/UFLP, res/Steiner/P6E/, and res/Steiner/small_i/) and execute the file `run.jl`.
