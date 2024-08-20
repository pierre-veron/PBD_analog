setwd("C:/Users/Utilisateur/Github/PBD_analog")


## Packages
library(ape)
library(treesliceR)
library(epm)

npar <- 5
nvar <- 5
nrep <- 500
df <- data.frame(par = rep(1:npar, each = nvar * nrep),
                 var = rep(1:nvar, each = nrep, times = nvar),
                 rep = rep((1:nrep)-1, times = npar*nvar))
df$fname <- paste0("stree_random-par", df$par, "-var", df$var, "-rep", df$rep, ".nwk")

dir <- "simulations_output/1-PBD/trees/"

pb <- txtProgressBar(max = nrow(df), style = 3)

drstats <- t(sapply(1:nrow(df), function(i) {
  filename <- paste0(dir, df[i, 'fname'])
  phy <- ape::read.tree(filename)

  #dr <- treesliceR::DR(phy)
  dr <- epm::DRstat(phy)
  # Save dr 
  write.csv(dr, file = paste0(dir, "DR/", df[i, 'fname'], ".csv"))
  setTxtProgressBar(pb, i)
  c(DR.mean =  mean(dr), DR.sd = sd(dr))
}))


close(pb)

df$DR.mean <- drstats[,'DR.mean']
df$DR.sd <- drstats[,'DR.sd']