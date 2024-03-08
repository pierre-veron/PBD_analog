# Compare the gamma statistics from PBD trees vs BD trees
## Path
``` 
PBD_bayesian_fit_CBD/8-gamma_tree_shape
```

## Overview
The $\gamma$ statistics is a summary of the tree shape. If $\gamma > 0$, the internal nodes are closer to the tips, and if $\gamma < 0$ the internal nodes are closer to the root, compared to a pure birth tree.

> Pybus, O. G. and Harvey, P. H. (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. _Proceedings of the Royal Society of London. Series B. Biological Sciences_, **267**, 2267--2272. DOI: [10.1098/rspb.2000.1278](https://doi.org/10.1098/rspb.2000.1278)


## Files
* `0-get_gamma.R`: parse the directory where the reconstructed trees from PBD simulation are stored and calculate their $\gamma$ statistics. 
* `1-simul_BD_get_gamma.R`: simulate birth-death trees with the equivalent birth-death rates from the PBD parameters used in the simulations and calculate their $\gamma$ statistics.
* `2-plot_gamma.ipynb`: plot the distributions of the $\gamma$ statistics from PBD trees _vs_ BD trees and compare them.

## Output 
* `0-get_gamma.R` generates a table of the $\gamma$ statistics in `strees_stats.csv`
* `1-simul_BD_get_gamma.R` generates the same table but for trees simulated under a BD model in `simulated_BD_trees_stats.csv`
* `2-plot_gamma.ipynb` generates the table of statistical tests `Gamma_PBD_BD_MannWhitney_U_test.csv` and the figure with boxplots `CBD_gamma_vs_BD_gamma.pdf`. 

## History 
* `2024/03/08` run the pipeline for the output of the simulation with cluster ID `12152`, Git version `b2cb4e1`. 