# ------------------------------------------------------------------------------
# simulate_with_PBD.R
# Simulates protracted birth death (PBD) trees with varying parameters and runs 
# inference of simple birth death (BD) models on the reconstructed species trees.
#
# Author: Pierre Veron, pierre.veron.2017[at]polytechnique.org 
# Created February 2024
# ------------------------------------------------------------------------------

# First argument : index of the simulation (integer)
# Second argument : path to the output directory
args_input = commandArgs(trailingOnly=TRUE)

# Initialization
outdir <- args_input[2]
library(PBD)
library(ape)


# Setting the default parameters 
param_default <- c(0.5, 1.0, 0.4, 0.2, 0.2) # PBD parameters 
param_n <- 5 # Number of values tested for each parameter
param_ranges <- list( # All values tested for each parameter
  c(10^-1, 10^-0.5, 1.0, 10^0.2, 10^0.4),
  10^(seq(from = -1, to = 1, length.out = param_n)),
  c(10^-1, 10^-0.5, 1.0, 10^0.2, 10^0.4),
  10^(seq(from = -1, to = 1, length.out = param_n)),
  10^(seq(from = -1, to = 1, length.out = param_n))
)

age <- 15 # Duration simulations

i_tree <- as.numeric(args_input[1])
print(i_tree)
set.seed(91 + 6 * i_tree)


# Simulation and inference 
simul_infer <- lapply(1:length(param_default), function(i_param) {
  print(paste0("Starting parameter ", i_param))
  param <- param_default
  
  as.data.frame(t(sapply(1:param_n, function(i_param_var) {
    
    param[i_param] <- param_ranges[[i_param]][i_param_var] * param_default[i_param]
    
    
    fname <- paste0("-par", i_param, "-var", i_param_var, "-rep", i_tree)
    
    # Simulate and save lineage tree
    phyPBD <- list()
    i_trials <- 0
    while (length(phyPBD) == 0 & i_trials < 100) {
      tryCatch({
        phyPBD <- PBD::pbd_sim(param, age = age, soc = 1)
      }, error= function(e) {})
      i_trials <- i_trials + 1 
    }
    if (length(phyPBD) > 0) {
      saveRDS(phyPBD$igtree.extinct, paste0(outdir, "/igtree.extinct", fname, ".rds"))
      
      # Get species tree 
      phySPE <- phyPBD$stree_random
      ape::write.tree(phySPE, paste0(outdir, "/stree_random", fname, ".nwk"))
      
      c("param_vary" = i_param, 
        "i_param_var" = i_param_var,
        "replicate" = i_tree, 
        "PBD.l1" = param[1], 
        "PBD.l2" = param[2],
        "PBD.l3" = param[3], 
        "PBD.mu1" = param[4],
        "PBD.mu2" = param[5],
        "SR" = length(phySPE$tip.label), 
        "PBD.LR.all" = length(phyPBD$igtree.extinct$tip.label),
        "PBD.LR.extant" = length(phyPBD$igtree.extant$tip.label))
    } else {
      c("param_vary" = i_param, 
        "i_param_var" = i_param_var,
        "replicate" = i_tree, 
        "PBD.l1" = param[1], 
        "PBD.l2" = param[2],
        "PBD.l3" = param[3], 
        "PBD.mu1" = param[4],
        "PBD.mu2" = param[5],
        "SR" = NA, 
        "PBD.LR.all" = NA,
        "PBD.LR.extant" = NA)
    }
    
    
  })))
})

# concatenate all 
df_simul_infer <- Reduce(function(x, y) merge(x, y, all=TRUE), simul_infer)  

# save results
write.csv(df_simul_infer, file = paste0(outdir, "/all_simulations_inference-rep-", i_tree, ".csv"))
