# Simulate constant birth-death (CBD) trees
## Path
```
2-simulate_trees/2-CBD
```

## Overview
Simulate CBD trees with rates calculated with the same parameters as for `1-PBD` (that is why this pipeline must be run before), _ie_:
$$ \hat\lambda = (1-\pi) \lambda_1$$
and 
$$ \hat\mu = \mu_1$$
with 
$\pi = \frac{\lambda_2 + \lambda_3 + \mu_2}{2\lambda_3} \left(  1 - \sqrt{ 1 - 4\frac{\lambda_3 \mu_2}{(\lambda_2 + \lambda_3 + \mu_2)^2}}  \right)$. 

## File 
* `1-simul_BD_get_gamma.R`: simulate the BD process and calculate the gamma statistics.

## Output 
The generated trees and statistics are stored in 
```
simulations_output/2-CBD
```