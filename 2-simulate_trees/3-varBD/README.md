# Simulate time-dependant birth-death (varBD) trees
## Path
```
2-simulate_trees/3-varBD
```

## Overview
Simulate the trees in the time-dependant BD model with rates calculated from the PBD parameters used in the pipeline `1-PBD` ; and calculates some statistics on these trees. 

## Files 
* `sim_parameters.csv`: a table of the PBD parameters used in the simulations
* `1-predict_variable_BD_rates.ipynb`: notebook to calculate the time-dependant BD rates based on the PBD parameters defined in the previous file. This script generates the output (in the directory `simulations_output/3-varBD`):
    * `variable_BDrates.npy` a `numpy` array containing the values of the rates for each time-step and each combination of parameters defined in `sim_parameters.csv`
* `2-simulate_varBD.Rmd` and `2bis-simulate_varBD.Rmd`: simulate the trees under the time-dependent BD rates with the rates calculated previously. The first script is for the cases where all the PBD rates vary independently, and the second script is for the simplified case where $\lambda_1 = \lambda_3$ and $\mu_1 = \mu_2$.

## Output 
The generated trees, data and statistics are stored in 
```
simulations_output/3-varBD
```