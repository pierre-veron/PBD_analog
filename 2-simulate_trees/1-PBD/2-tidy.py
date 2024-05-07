import pandas as pd
import numpy as np
from matplotlib.cbook import boxplot_stats
import json

# Adapt this section
outdir = "C:/Users/pveron/Output_clusters/PBD_analog/12152"
n_replicates = 200
mcmc_size = 5000

list_df = []

np.random.seed(266)

for i in range(n_replicates):
    list_df.append(pd.read_csv(outdir + "/all_simulations_inference-rep-{}.csv".format(i)))

df = pd.concat(list_df, ignore_index=True)

df.to_csv(outdir + "/all_simulations_inference.csv")

print("Concatenated {} replicates.".format(n_replicates))
