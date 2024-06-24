#!/usr/bin/env Rscript

# Load required libraries
library(RPANDA)
library(ape)
library(phytools)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide par_row as input argument.")
}
par_row <- as.integer(args[1])

# Function to slice the tree
slice_tree <- function(tree, time_stop, tot_time){
  sliced_tree <- tree
  sliced_sub_trees <- treeSlice(sliced_tree, slice = tot_time - time_stop, trivial = TRUE)
  
  for (i in 1:length(sliced_sub_trees)){
    if (Ntip(sliced_sub_trees[[i]]) > 1){
      sliced_tree <- drop.tip(sliced_tree, tip = sliced_sub_trees[[i]]$tip.label[2:Ntip(sliced_sub_trees[[i]])])
    }
  }
  
  for (i in which(node.depth.edgelength(sliced_tree) > (tot_time - time_stop))){
    sliced_tree$edge.length[which(sliced_tree$edge[,2] == i)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2] == i)] - time_stop
  }
  return(sliced_tree)
}

# Function to compute equivalent BD rates
equivalent_bd_rates <- function(param) {
  l1 <- param[1]
  l2 <- param[2]
  l3 <- param[3]
  m1 <- param[4]
  m2 <- param[5]
  p <- 0.5 * (l2 + l3 + m2) / l3 * (1 - sqrt(1 - 4 * l3 * m2 / ((l2 + l3 + m2)^2)))
  l <- (1 - p) * l1
  m <- m1
  rates <- c(l, m)
  names(rates) <- NULL
  return(rates)
}

# Read parameters
pars_df <- read.csv("simulations_output/1-PBD/all_simulations_inference.csv")
pars <- pars_df[par_row,]

# Load the tree
tree_file <- sprintf("simulations_output/1-PBD/trees/stree_random-par%s-var%s-rep%s.nwk", 
                     pars$param_vary, pars$i_param_var, pars$replicate)
if (!file.exists(tree_file)) {
  stop(paste("Tree file not found:", tree_file))
}

PBD_tree <- read.tree(tree_file)
tot_time <- max(node.age(PBD_tree)$ages)

# Compute and print equivalent BD rates
if (nrow(pars) > 0) {
  eq_BD <- equivalent_bd_rates(pars[c("PBD.l1", "PBD.l2", "PBD.l3", "PBD.mu1", "PBD.mu2")])
  print(paste("Equivalent BD rates:", paste(eq_BD, collapse = ", ")))
} else {
  print("No matching parameter row found")
}

# Fit the constant-rate birth-death model
f.lamb <- function(t, y) { y[1] }
f.mu <- function(t, y) { y[1] }
lamb_par <- c(0.5)
mu_par <- c(0.2)

fit_full <- fit_bd(PBD_tree, tot_time, f.lamb, f.mu, lamb_par, mu_par)

split_times <- c(2**((1:6)/7)-1, 2**((0:10)/5))
#split_times <- 1
fits_truncated <- sapply(split_times, function(time_stop) {
  sliced_tree <- slice_tree(PBD_tree, time_stop, tot_time)
  result_cst <- fit_bd_in_past(sliced_tree, tot_time, time_stop, f.lamb, f.mu, 
                               desc = Ntip(PBD_tree), tot_desc = Ntip(PBD_tree), lamb_par, mu_par, dt = 1e-3)
  return(result_cst)
})

result <- rbind(split_times=c(NA,0,split_times), 
                cbind(c("equiv. BD", unlist(eq_BD)), 
                      fit_full[c("model", "lamb_par", "mu_par")], 
                      fits_truncated[c("model", "lamb_par", "mu_par"),]))

# Save the result to a file
write.csv(result, sprintf("simulations_output/4-Truncated_PBD/Fits/FitBD_results_par%s_var%s_rep%s.csv",
                          pars$param_vary, pars$i_param_var, pars$replicate))

