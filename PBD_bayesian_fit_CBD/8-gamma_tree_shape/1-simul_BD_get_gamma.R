library(ape)

setwd("C:/Users/pveron/Output_clusters/PBD_analog/12152")
age <- 15

equivalent_bd_rates <- function(param) {
  l1 <- param[1]
  l2 <- param[2]
  l3 <- param[3]
  m1 <- param[4]
  m2 <- param[5]
  D <- sqrt((l2+l3)^2 + 2*(l2-l3)*m2+m2^2)
  phi <- l2-l3+m2
  p <- 0.5*(l2+l3+m2)/l3 * (1-sqrt(1-4*l3*m2/((l2+l3+m2)^2)))
  den <- 1 + (m1 + (1-p)*l1) * 2 * log(2*D/(D+phi)) / (D-phi) 
  l <- (1-p)*l1 / den
  m <- m1 / den
  rates <- c(l, m)
  names(rates) <- NULL
  rates
}

simul_infer <- read.csv("all_simulations_inference.csv")

param_PBD_names <-  paste0("PBD.", c("l1", "l2", "l3", "mu1", "mu2"))
col_to_keep <- c("param_vary", "i_param_var", "replicate", param_PBD_names)


tree_stats_df <- as.data.frame(t(sapply(rownames(simul_infer), function(rw) {
  param_PBD <- unlist(simul_infer[rw, param_PBD_names])
  eq_bd <- equivalent_bd_rates(param_PBD)
  tree <- ape::rbdtree(eq_bd[1], eq_bd[2], age)
  c(unlist(simul_infer[rw, col_to_keep]),
    "gamma" = ape::gammaStat(tree),
    "SR" = length(tree$tip.label),
    "equiv_birth" = eq_bd[1],
    "equiv_death" = eq_bd[2])
})))

write.csv(tree_stats_df, "simulated_BD_trees_stats.csv")