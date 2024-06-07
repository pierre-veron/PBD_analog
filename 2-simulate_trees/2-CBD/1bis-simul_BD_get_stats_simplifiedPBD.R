library(ape)
library(treestats)
library(TreeSim)

#setwd("C:\Users\pveron\Documents\GitHub\PBD_analog")
#setwd("Users/Jeremy/Nextcloud/Recherche/1_Methods/PBD_analog")

set.seed(394)

age <- 15
n_trees  <- 500
n_trials <- 100
n_par <- 3
n_val <- 5

equivalent_bd_rates <- function(param) {
  l1 <- param[1]
  l2 <- param[2]
  l3 <- param[3]
  m1 <- param[4]
  m2 <- param[5]
  p <- 0.5*(l2+l3+m2)/l3 * (1-sqrt(1-4*l3*m2/((l2+l3+m2)^2)))
  l <- (1-p)*l1
  m <- m1
  rates <- c(l, m)
  names(rates) <- NULL
  rates
}

simul_infer <- read.csv("simulations_output/1-PBD/all_simulations_inference_simplified_PBD.csv")
param_PBD_names <-  paste0("PBD.", c("l1", "l2", "l3", "mu1", "mu2"))
sim_pars <- simul_infer[simul_infer$replicate==0, c("param_vary", "i_param_var", param_PBD_names)]
write.csv(sim_pars, "2-simulate_trees/3-varBD/sim_parameters_simplified.csv")

# Case of simulations where only one parameter varies
#if ("param_vary" %in% colnames(simul_infer)) {
#  col_to_keep <- c("param_vary", "i_param_var", "replicate", param_PBD_names)
#} else {
#  col_to_keep <- c("i_param_var", "replicate", param_PBD_names)
#}

sim_tree_BD <- function(lambda, mu){
  tree <- sim.bdsky.stt(n = 0, lambdasky = lambda, deathsky = mu, timesky = 0.0, sampprobsky = 0, rho = 1, timestop = age)[[1]]
  
  # Single tip
  if (is.numeric(tree) && tree == 1){
    #branch <- list(edge =  matrix(c(2, 1), ncol = 2), edge.length = age, tip.label = "t1", Nnode = 1, root.edge = 0)
    #class(branch) <- "phylo"
    #return(branch)
    return(-1)
  }
  
  return(tree)
}

trees <- apply(sim_pars[param_PBD_names], 1, function(param_PBD){
  eq_bd <- equivalent_bd_rates(param_PBD)
  lapply(1:n_trees, function(n){
      i <- 0
      while (i<n_trials){
        if (eq_bd[1] - eq_bd[2] > log(10000)/age) return (NULL) ## Remove trees whose expected size is too large
        tree <- sim_tree_BD(lambda = eq_bd[1], mu = eq_bd[2])
        if (is.phylo(tree)) return (tree)
        i <- i+1
      }
  })
})

Ntips <- sapply(trees, sapply, function(result)if(is.null(result)){NA}else{length(result$tip.label)})

trees_gammma <- sapply(trees, sapply, function(result){
  if (is.null(result)){
    return(NaN)
  }else{
    return(ape::gammaStat(result))
  }
})

trees_stairs2 <- sapply(trees, sapply, function(result){
  if (is.null(result)){
    return(NaN)
  }else{
    return(treestats::stairs2(result))
  }
})

trees_beta <- sapply(trees, sapply, function(result){
  if (is.null(result)){
    return(NaN)
  }else{
    return(treestats::beta_statistic(result))
  }
})

trees_stats <- data.frame(param_vary = rep(rep(1:n_par, each=n_val), n_trees), 
                          i_param_var = rep(1:n_val, n_trees), 
                          replicate = rep(0:(n_trees-1), each=nrow(sim_pars)),
                          gamma = as.numeric(t(trees_gammma)),
                          stairs2 = as.numeric(t(trees_stairs2)),
                          beta = as.numeric(t(trees_beta)),
                          SR = as.integer(t(Ntips)))

write.csv(trees_stats, "simulations_output/2-CBD/simulated_BD_trees_stats_simplified_PBD.csv")

for (i in 1:length(trees)){
  trees_nonNull <- Filter(Negate(is.null), trees[[i]])
  if (length(trees_nonNull) > 0)  write.tree(trees_nonNull, paste0("simulations_output/2-CBD/trees_simplified/trees_row_", i))
}
