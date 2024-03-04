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
library(coda)
library(RPANDA)
library(diversitree)


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
    phyPBD <- PBD::pbd_sim(param, age = age)
    saveRDS(phyPBD$igtree.extinct, paste0(outdir, "/igtree.extinct", fname, ".rds"))
    
    # Get species tree 
    phySPE <- phyPBD$stree_random
    ape::write.tree(phySPE, paste0(outdir, "/stree_random", fname, ".nwk"))
    
    if (length(phySPE$tip.label) < 3) {
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
        "PBD.LR.extant" = length(phyPBD$igtree.extant$tip.label), 
        "ML.l" = NA,
        "ML.mu" = NA,
        "MCMC.l.mean" = NA,
        "MCMC.l.median" = NA,
        "MCMC.l.sd" = NA,
        "MCMC.l.q25" = NA,
        "MCMC.l.q75" = NA,
        "MCMC.l.ess" = NA,
        "MCMC.mu.mean" = NA,
        "MCMC.mu.median" = NA,
        "MCMC.mu.sd" = NA,
        "MCMC.mu.q25" = NA,
        "MCMC.mu.q75" = NA,
        "MCMC.mu.ess" = NA,
        "MCMC.div.mean" = NA,
        "MCMC.div.median" = NA,
        "MCMC.div.sd" = NA,
        "MCMC.div.q25" = NA,
        "MCMC.div.q75" = NA,
        "MCMC.turnov.mean" = NA,
        "MCMC.turnov.median" = NA,
        "MCMC.turnov.sd" = NA,
        "MCMC.turnov.q25" =NA,
        "MCMC.turnov.q75" = NA
      )
    } else {
      # Infer BD rates
      lik_bd <- diversitree::make.bd(phySPE)
      fit_bd <- diversitree::find.mle(lik_bd, x.init = c(param[1], param[4]), 
                                      method="nlminb", lower=0, condition.surv=TRUE)
      prior_bd <- diversitree::make.prior.exponential(1)
      samples <- diversitree::mcmc(lik_bd, fit_bd$par, nsteps=5000, prior=prior_bd,
                                   lower=c(0, 0), upper=c(Inf, Inf), w=c(.1, .1),
                                   fail.value=-Inf, print.every=0)
      
      # remove burnin, extract rates
      lambda <- samples$lambda
      mu <- samples$mu  
      lambda_cv <- lambda[(0.1 * length(lambda)):length(lambda)]
      mu_cv <- mu[(0.1 * length(mu)):length(mu)]
      div_cv <- lambda_cv - mu_cv
      turnov_cv <- mu_cv / lambda_cv
      
      # save mcmc samples
      write.csv(samples, paste0(outdir, "/mcmc", fname, ".csv"))
      
      # return as one line the simulation
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
        "PBD.LR.extant" = length(phyPBD$igtree.extant$tip.label),
        "ML.l" = fit_bd$par['lambda'][[1]],
        "ML.mu" = fit_bd$par['mu'][[1]],
        "MCMC.l.mean" = mean(lambda_cv), 
        "MCMC.l.median" = median(lambda_cv),
        "MCMC.l.sd" = sd(lambda_cv), 
        "MCMC.l.q25" = quantile(lambda_cv,0.25)[[1]],
        "MCMC.l.q75" = quantile(lambda_cv, 0.75)[[1]], 
        "MCMC.l.ess" = coda::effectiveSize(lambda_cv)[[1]],
        "MCMC.mu.mean" = mean(mu_cv), 
        "MCMC.mu.median" = median(mu_cv),
        "MCMC.mu.sd" = sd(mu_cv), 
        "MCMC.mu.q25" = quantile(mu_cv,0.25)[[1]],
        "MCMC.mu.q75" = quantile(mu_cv, 0.75)[[1]], 
        "MCMC.mu.ess" = coda::effectiveSize(mu_cv)[[1]],
        "MCMC.div.mean" = mean(div_cv), 
        "MCMC.div.median" = median(div_cv),
        "MCMC.div.sd" = sd(div_cv), 
        "MCMC.div.q25" = quantile(div_cv,0.25)[[1]],
        "MCMC.div.q75" = quantile(div_cv, 0.75)[[1]], 
        "MCMC.turnov.mean" = mean(turnov_cv), 
        "MCMC.turnov.median" = median(turnov_cv),
        "MCMC.turnov.sd" = sd(turnov_cv), 
        "MCMC.turnov.q25" = quantile(turnov_cv,0.25)[[1]],
        "MCMC.turnov.q75" = quantile(turnov_cv, 0.75)[[1]])
    }
    
    
  })))
})

# concatenate all 
df_simul_infer <- Reduce(function(x, y) merge(x, y, all=TRUE), simul_infer)  

# save results
write.csv(df_simul_infer, file = paste0(outdir, "/all_simulations_inference-rep-", i_tree, ".csv"))
