setwd("C:/Users/pveron/Documents/Github/PBD_analog")


## Packages
library(ape)
library(treesliceR)
library(epm)


dir <- "simulations_output/3-varBD/trees/"

n_config <- 25
pb <- txtProgressBar(max = n_config, style = 3)

sapply(1:n_config, function(i) {
  trees <- ape::read.tree(paste0(dir, "sim_varBD_pars", i, ".trees"))
  if(class(trees) == "phylo") { # only one tree in the list
    print(i)
    trees <- list(trees)
  }
  lapply(1:length(trees), function(i_tree) {
    phy <- trees[[i_tree]]
    dr <- epm::DRstat(phy)
    write.csv(dr, file = paste0(dir, "/DR/config-", i, "-tree-", i_tree, ".csv"))
  })
  setTxtProgressBar(pb, i)
})
close(pb)
  