library(ape)

setwd("C:\Users\pveron\Documents\GitHub\PBD_analog")

set.seed(394)

age <- 15

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

simul_infer <- read.csv("simulations_output/PBD/all_simulations_inference.csv")

param_PBD_names <-  paste0("PBD.", c("l1", "l2", "l3", "mu1", "mu2"))


# Case of simulations where only one parameter varies
if ("param_vary" %in% colnames(simul_infer)) {
  col_to_keep <- c("param_vary", "i_param_var", "replicate", param_PBD_names)
} else {
  col_to_keep <- c("i_param_var", "replicate", param_PBD_names)
}




tree_stats_df <- as.data.frame(t(sapply(rownames(simul_infer), function(rw) {
  param_PBD <- unlist(simul_infer[rw, param_PBD_names])
  eq_bd <- equivalent_bd_rates(param_PBD)
  tryCatch(
    {
    tree <- ape::rbdtree(eq_bd[1], eq_bd[2], age, eps = 1e-04)
    fname <- paste0("CBD_tree_sim_", rw, "_b_", eq_bd[1], "_d_", eq_bd[2],  ".nwk")
    ape::write.tree(tree, fname)
    out <- c(unlist(simul_infer[rw, col_to_keep]),
      "gamma" = ape::gammaStat(tree),
      "SR" = length(tree$tip.label),
      "equiv_birth" = eq_bd[1],
      "equiv_death" = eq_bd[2])
    out
    },
    error = function(e) {
      out <- c(unlist(simul_infer[rw, col_to_keep]),
               "gamma" = NA,
               "SR" = NA,
               "equiv_birth" = eq_bd[1],
               "equiv_death" = eq_bd[2])
      out
    }
  )
})))

write.csv(tree_stats_df, "simulations_output/stats.csv")