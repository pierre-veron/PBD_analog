# Get the gamma statistics from a directory of trees 
library(ape)

setwd("C:/Users/pveron/Output_clusters/PBD_analog/12152")

n_par <- 5
n_val <- 5
n_rep <- 200

list_out <- lapply(1:n_par, function(i_par) {
  list_out <- lapply(1:n_val, function(i_val) {
    out <- t(sapply(0:(n_rep-1), function(i_rep) {
      fname <- paste0("stree_random-par", i_par, "-var", i_val, "-rep", i_rep, ".nwk")
      
      tree <- ape::read.tree(fname)
      c("param_vary" = i_par,
        "i_param_var" = i_val, 
        "replicate" = i_rep, 
        "gamma" = ape::gammaStat(tree))
    }))
    out
  })
  
  
  Reduce(function(x, y) merge(x, y, all=TRUE), list_out)
})


df <- Reduce(function(x, y) merge(x, y, all=TRUE), list_out)


write.csv(df, "strees_stats.csv")