import pandas as pd
import numpy as np


outdir = "C:/Users/pveron/Output_clusters/PBD_analog/12038"
n_replicates = 100

list_df = []

for i in range(n_replicates):
    list_df.append(pd.read_csv(outdir + "/all_simulations_inference-rep-{}.csv".format(i)))

df = pd.concat(list_df, ignore_index=True)

df.to_csv(outdir + "/all_simulations_inference.csv")

print("Concatenated {} replicates.".format(n_replicates))


# Summarize all MCMCs
burnin = 500
n_param_var = int(max(df["i_param_var"]))
par_names_PBD = ["PBD." + s for s in ['l1','l2','l3', 'mu1', 'mu2']]
keys_keep = ['param_vary', 'i_param_var'] + par_names_PBD

list_df = []

for i in range(len(par_names_PBD)):
    for i_param_var in range(1, 1+n_param_var):

        par_simul = (df.loc[(df.param_vary == i + 1) & (df.i_param_var == i_param_var) & (df.replicate == 1)]).to_dict("records")[0]
        par_simul = {k:par_simul[k] for k in keys_keep}

        birth_all, death_all = [],[]

        # Load MCMC samples
        n_trees_mcmc = 0
        for i_tree in range(n_replicates):
            fname = "mcmc-par{}-var{}-rep{}.csv".format(i+1, i_param_var, i_tree+1)
            try:
                samples = pd.read_csv(outdir + "/" + fname)
                samples = samples.loc[burnin:]
            except:
                pass
            else:
                n_trees_mcmc += 1
                birth_all += list(samples["lambda"])
                death_all += list(samples["mu"])

        birth_all, death_all = np.array(birth_all), np.array(death_all)

        par_simul["n_trees_MCMC"] = n_trees_mcmc

        par_simul["allMCMC.l.mean"] = np.mean(birth_all)
        par_simul["allMCMC.l.median"] = np.median(birth_all)
        par_simul["allMCMC.l.sd"] = np.std(birth_all)
        par_simul["allMCMC.l.q25"] = np.quantile(birth_all, 0.25)
        par_simul["allMCMC.l.q75"] = np.quantile(birth_all, 0.75)
        par_simul["allMCMC.l.min"] = np.min(birth_all)
        par_simul["allMCMC.l.max"] = np.max(birth_all)

        par_simul["allMCMC.mu.mean"] = np.mean(death_all)
        par_simul["allMCMC.mu.median"] = np.median(death_all)
        par_simul["allMCMC.mu.sd"] = np.std(death_all)
        par_simul["allMCMC.mu.q25"] = np.quantile(death_all, 0.25)
        par_simul["allMCMC.mu.q75"] = np.quantile(death_all, 0.75)
        par_simul["allMCMC.mu.min"] = np.min(death_all)
        par_simul["allMCMC.mu.max"] = np.max(death_all)

        div_all = birth_all - death_all
        turnov_all = death_all / birth_all

        par_simul["allMCMC.div.mean"] = np.mean(div_all)
        par_simul["allMCMC.div.median"] = np.median(div_all)
        par_simul["allMCMC.div.sd"] = np.std(div_all)
        par_simul["allMCMC.div.q25"] = np.quantile(div_all, 0.25)
        par_simul["allMCMC.div.q75"] = np.quantile(div_all, 0.75)
        par_simul["allMCMC.div.min"] = np.min(div_all)
        par_simul["allMCMC.div.max"] = np.max(div_all)

        par_simul["allMCMC.turnov.mean"] = np.mean(turnov_all)
        par_simul["allMCMC.turnov.median"] = np.median(turnov_all)
        par_simul["allMCMC.turnov.sd"] = np.std(turnov_all)
        par_simul["allMCMC.turnov.q25"] = np.quantile(turnov_all, 0.25)
        par_simul["allMCMC.turnov.q75"] = np.quantile(turnov_all, 0.75)
        par_simul["allMCMC.turnov.min"] = np.min(turnov_all)
        par_simul["allMCMC.turnov.max"] = np.max(turnov_all)

        list_df.append(par_simul)

df = pd.DataFrame(list_df)

df.to_csv(outdir + "/summary_all_MCMC.csv")

print("Summarized all replicates of MCMC samples.")