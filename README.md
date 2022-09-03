Data and scripts for geometric feature sizes
============================================

Copyright (C) 2022 [Parker
Edwards](https://sites.nd.edu/parker-edwards/)


External requirements
---------------------

1. [Julia](https://julialang.org/downloads/>). Written and tested with version 1.4.1.
2. [Jupyter notebook](https://jupyter.org/install)
3. [HomologyInferenceWithWeakFeatureSize.jl](https://github.com/P-Edwards/HomologyInferenceWithWeakFeatureSize.jl). 
4. (Optional) If you would like to re-compute bottleneck correspondence solutions, [Bertini](https://bertini.nd.edu/), parallel version recommended.


Setup
-----
1. Clone this repository
2. In a terminal, boot Julia, go to package mode with `]`. There are now two options
3. A - If you are not worried about e.g. `HomotopyContinuation.jl` in your main Julia environment, you can just run `add https://github.com/P-Edwards/HomologyInferenceWithWeakFeatureSize.jl`
3. B - If you are concerned about this. First run `activate .`. Then run `add https://github.com/P-Edwards/HomologyInferenceWithWeakFeatureSize.jl`. 
4. (Optional) If you would like to run the `butterfly_lfs_sparse` example, also `add Distances, PersistenceDiagrams, ProgressBars, JLD, CSV, Measures, Plots, DataFrames`
5. B - With option B, a `Project.toml` file was created. Copy this into any directory with a notebook you would like to run yourself.

Description and Usage
---------------------
This repository contains data and scripts accompanying the manuscript [Computing geometric feature sizes for algebraic manifolds](a/link/to/arxiv). The directory structure is sorted by example. To recompute Julia-based filtering steps with precomputed solutions, install the first three requirements and: 

1. Start jupyter notebook via terminal with `jupyter notebook`
2. Navigate in the browser tab opened to a notebook file in one of the provided directories. 

If you would like to check the algebraic computations, each directory for an example which uses Bertini includes a directory `bertini_inputs`. Each leaf folder contains a file labelled `input`, and you can run Bertini, say with 12 threads, via `mpirun -np 12 bertini input`. Follow directions in the notebooks to filter these new solutions.


Large data files
----------------
Some data files are too large to store on GitHub. They are instead available for download at [https://notredame.box.com/v/feature-size-large-data](https://notredame.box.com/v/feature-size-large-data). For each directory at that link, place the data in the "paper_data" folder for this repository's correspondingly named directory.


Computational warning
---------------------
Algebraic computations for 3- and 4-bottlenecks for the examples in R^3, i.e. the quartic surface and torus-Clebsch intersection curve, are expensive. The 4-bottlenecks in particular can take a week or more to re-compute with Bertini. 


License
-------
This repository is distributed under GPLv3. 