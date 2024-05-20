# Simulate PBD process and calculate some statistics
## Path 
```
2-simulate_trees/1-PBD
```

## Overview
This pipeline runs the PBD model for several replicates and combination of parameters. Scripts 0 and 1 are intended to be run on a cluster. 

## Files
* `0-install_packages.R`: install the required packages on the used environment
* `1-simulate_with_PBD.R`: simulate the trees and infer the BD rates. This script only runs 1 replicate for several values for the parameters $\lambda_1, \lambda_2, \lambda_3, \mu_1, \mu2$. Takes external arguments:
    * id of the replicate (change the seed and name of the stored files)
    * path to store the output files. 
    
    For instance 
```
Rscript ./1-simulate_with_PBD.R 1 /data/biodiv/pveron/PBD_analog/out
```
* `1-wrapper.sh`: a wrapper to call the script `1-simulate_with_PBD.R` in the right environment and with the right arguments
* `1-submission.sub`: submission file to run the whole process on the cluster, with several replicates, each replicates changes the id that is passed as argument to the script `1-simulate_with_PBD.R`.
* `2-tidy.py`: script to put all the replicates in a unique `.csv` file. 
* `3-get_stats.R`: calculate statistics on the trees. 

## Output 
The generated trees and statistics are stored in 
```
simulations_output/1-PBD
```

## History 
* `2024-05-19`: Run simulations on the cluster for `1-simulate_with_PBD.R` (cluster ID `12946`, Git version ` d5a5323`)