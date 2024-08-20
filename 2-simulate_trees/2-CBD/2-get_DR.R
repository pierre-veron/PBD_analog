setwd("C:/Users/pveron/Documents/Github/PBD_analog")


## Packages
library(ape)
library(treesliceR)
library(epm)


dir <- "simulations_output/2-CBD/trees/"
lf <- list.files(dir, pattern = "*.nwk")
pb <- txtProgressBar(max = length(lf), style = 3)

drstats <- t(sapply(1:length(lf), function(i) {
  filename <- lf[i]
  filepath <- paste0(dir, filename)
  phy <- ape::read.tree(filepath)

  #dr <- treesliceR::DR(phy)
  dr <- epm::DRstat(phy)
  # Save dr 
  write.csv(dr, file = paste0(dir, "DR/", filename, ".csv"))
  setTxtProgressBar(pb, i)
  c(DR.mean =  mean(dr), DR.sd = sd(dr))
}))


close(pb)

df$DR.mean <- drstats[,'DR.mean']
df$DR.sd <- drstats[,'DR.sd']