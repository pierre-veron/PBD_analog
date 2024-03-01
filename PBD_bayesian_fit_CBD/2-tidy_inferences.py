import pandas as pd
import numpy as np
from matplotlib.cbook import boxplot_stats
import json

# Adapt this section
outdir = "C:/Users/pveron/Output_clusters/PBD_analog/12152"
two_tips_posterior = "INSANE_CBD_2tips_trees/outputs/CBD_2tips_50000000iter_rootAge15_seed0.log"
n_replicates = 200
mcmc_size = 5000

list_df = []

np.random.seed(266)

for i in range(n_replicates):
    list_df.append(pd.read_csv(outdir + "/all_simulations_inference-rep-{}.csv".format(i)))

df = pd.concat(list_df, ignore_index=True)

df.to_csv(outdir + "/all_simulations_inference.csv")

print("Concatenated {} replicates.".format(n_replicates))


# Load the posteriors for the 2-tips trees
two_tips = pd.read_table(two_tips_posterior, sep = "\t")

# Summarize all MCMCs
burnin = 500
n_param_var = int(max(df["i_param_var"]))
par_names_PBD = ["PBD." + s for s in ['l1','l2','l3', 'mu1', 'mu2']]
keys_keep = ['param_vary', 'i_param_var'] + par_names_PBD

list_df = []

for i in range(len(par_names_PBD)):
    for i_param_var in range(1, 1+n_param_var):

        par_simul = (df.loc[(df.param_vary == i + 1) & (df.i_param_var == i_param_var) & (df.replicate == 1)]).to_dict("records")[0]
        df_simul = df.loc[(df.param_vary == i + 1) & (df.i_param_var == i_param_var)]
        par_simul = {k:par_simul[k] for k in keys_keep}

        birth_all, death_all = [],[]

        # Load MCMC samples
        n_trees_non_trivial = 0

        for i_tree in range(n_replicates):
            sr = list(df_simul.loc[df_simul.replicate == i_tree, "SR"])[0]
            if sr > 2:
                fname = "mcmc-par{}-var{}-rep{}.csv".format(i+1, i_param_var, i_tree)
                samples = pd.read_csv(outdir + "/" + fname)
                samples = samples.loc[burnin:]
                n_trees_non_trivial += 1
            else:
                subsampling = np.random.choice(two_tips.index, size = mcmc_size - burnin, replace = False)
                samples = two_tips.loc[subsampling]
                
            birth_all += list(samples["lambda"])
            death_all += list(samples["mu"])

        birth_all, death_all = np.array(birth_all), np.array(death_all)
        div_all = birth_all - death_all
        turnov_all = death_all / birth_all

        # save statistics for all type of rate 
        rate_names = ["l", "mu", "div", "turnov"]
        rate_data = [birth_all, death_all, div_all, turnov_all]
        par_simul["n_trees_non_trivial"] = n_trees_non_trivial

        for name, data in zip(rate_names, rate_data):
            par_simul["allMCMC."+name+".mean"] = np.mean(data)
            par_simul["allMCMC."+name+".median"] = np.median(data)
            par_simul["allMCMC."+name+".sd"] = np.std(data)
            par_simul["allMCMC."+name+".q25"] = np.quantile(data, 0.25)
            par_simul["allMCMC."+name+".q75"] = np.quantile(data, 0.75)
            par_simul["allMCMC."+name+".q025"] = np.quantile(data, 0.025)
            par_simul["allMCMC."+name+".q975"] = np.quantile(data, 0.975)
            par_simul["allMCMC."+name+".min"] = np.min(data)
            par_simul["allMCMC."+name+".max"] = np.max(data)

            # boxplot data 
            box = boxplot_stats(data)[0]
            box['fliers'] = list(box['fliers'])
            with open(outdir + "/boxplot-"+name+"-par{}-var{}.json".format(i+1, i_param_var), 'w') as f: 
                json.dump(box, f)
        
        list_df.append(par_simul)

df = pd.DataFrame(list_df)

df.to_csv(outdir + "/summary_all_MCMC.csv")

print("Summarized all replicates of MCMC samples.")