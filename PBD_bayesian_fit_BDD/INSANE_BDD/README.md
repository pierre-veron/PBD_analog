# scripts to fit a birth-death-diffusion (BDD) model to trees simulated from the PBD model
## Path
``` 
PBD_bayesian_fit_BDD/INSANE_BDD
```

## Overview

## Files 
* `0-create_list_trees.py`: generates the list of trees names by parsing the directory where the reconstructed trees are saved.
* `1-ConstantPBDtrees_BDD.jl`: fit a BDD model to a tree passed in argument 
* `1-ConstantPBDtrees_BDD.submit`: submission file to run the inference on each tree in the file `list_trees.txt`
* `Source_INSANE.jl`: build a custom version of Insane to cope with exponential priors

## Run 
Clone the version of `Tapestree` on the cluster. 
```
git clone -b insane https://github.com/Jeremy-Andreoletti/Tapestree.jl.git
```
In `Source_INSANE.jl`, adapt the variable called `Tapestree_path` to make it point to the directory where the repository `Tapestree.jl`was cloned. 

To run the inference on a specific tree, simply run 
```
julia 1-ConstantPBDtrees_BDD.jl <tree_name>
```

To run the inference on multiple trees, run the script to parse the folder where the reconstructed trees are stored. 
```
python 0-create_list_trees.py
# will create or overwrite the file list_trees.txt
```
Then run on the cluster 
```
condor_submit 1-ConstantPBDtrees_BDD.submit
```
This will submit 1 job per line in `list_trees.txt`.

## Outputs 
* `outputs/<X>-BDD_ConstantPBDtrees_<Y>iter_seed<Z>.log`: MCMC samples of the parameters of the BDD model, with `X`=tree name, `Y`= number of iterations, `Z`= random seed used to run the MCMC. 
* `outputs/<X>-BDD_ConstantPBDtrees_<Y>iter_seed<Z>.txt`: MCMC samples of the augmented trees, with the corresponding rates. 

No image or animation is generated since the corresponding lines in `1-ConstantPBDtrees_BDD.jl` is commented. 