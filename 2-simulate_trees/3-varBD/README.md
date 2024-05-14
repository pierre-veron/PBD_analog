# Simulate time-dependant birth-death (varBD) trees
## Path
```
2-simulate_trees/3-varBD
```

## Overview
Simulate the trees in the time-dependant BD model with rates calculated from the PBD parameters used in the pipeline `1-PBD` ; and calculates some statistics on these trees. 

## Files 
* `sim_parameters.csv`: a table of the PBD parameters used in the simulations
* `1-predict_variable_BD_rates.ipynb`: notebook to calculate the time-dependant BD rates based on the PBD parameters defined in the previous file. This script generates the output:
    * `variable_BDrates.npy` a `numpy` array containing the values of the rates for each time-step and each combination of parameters defined in `sim_parameters.csv`
* `2-Compare_PBD_varBD_trees.Rmd`: a notebook to generate the trees based on the time-dependant BD model and calculate some statistics on those trees. 

## Output 
The generated trees and statistics are stored in 
```
simulations_output/3-varBD
```