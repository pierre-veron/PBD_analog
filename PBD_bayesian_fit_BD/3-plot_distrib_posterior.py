# Plot distribution of posterior samples from MCMC for each tree

# Import 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
plt.style.use("customstyle")

import sys
sys.path.append("../modules")

import PBD_analog

colors = ["#" + x for x in ["000000","cf5c36","bcd696","985f99","9684a1"]]


# Load the results of the inference
outdir = "C:/Users/pveron/Output_clusters/PBD_analog/12038"

simul_infer = pd.read_csv(outdir + "/all_simulations_inference.csv")
simul_infer["combinaison"] = 10*simul_infer.param_vary + simul_infer.i_param_var

simul_infer_g = simul_infer.groupby("combinaison")
mean_simul = simul_infer_g.mean()
std_simul = simul_infer_g.std()


# Set different bins for birth/div rate and death/turnover rate
bins1 = np.linspace(0, 2, 40) # birth rates
bins2 = np.linspace(0, 6, 40) # death rates
bins3 = np.linspace(-1, 1, 40) # diversification rates
bins4 = np.linspace(0, 5, 40) # turnover

bins1_c = (bins1[1:] + bins1[:-1]) / 2 
bins2_c = (bins2[1:] + bins2[:-1]) / 2
bins3_c = (bins3[1:] + bins3[:-1]) / 2 
bins4_c = (bins4[1:] + bins4[:-1]) / 2

burnin = 500

bins = [bins1, bins2, bins3, bins4]

# Make the plots 
for param_vary in range(1, 6):
    for i_param_var in range(1, 6):

        n_rep = np.sum((simul_infer.param_vary == param_vary) & (simul_infer.i_param_var == i_param_var))

        par_simul = (simul_infer.loc[(simul_infer.param_vary == param_vary) & (simul_infer.i_param_var == i_param_var) & (simul_infer.replicate == 1)]).to_dict("records")[0]

        all_infer = [[],[],[],[]]

        # Load MCMC samples
        mcmc = []
        no_mcmc = []
        for i_tree in range(n_rep):
            fname = "mcmc-par{}-var{}-rep{}.csv".format(param_vary, i_param_var, i_tree+1)
            try:
                samples = pd.read_csv(outdir + "/" + fname)
                samples = samples.loc[burnin:]
            except:
                no_mcmc.append(i_tree + 1)
            else:
                mcmc.append(samples)



        fig, axes = plt.subplots(2, 2, figsize = (10, 10))
        axes = axes.flatten()
        (ax1, ax2, ax3, ax4) = axes
        # calculate estimates
        est_l, est_mu = PBD_analog.analog_BD_rates(l1 = par_simul["PBD.l1"], 
                                        l2 = par_simul["PBD.l2"],
                                        l3 = par_simul["PBD.l3"], 
                                        m1 = par_simul["PBD.mu1"],
                                        m2 = par_simul["PBD.mu2"])

        # Prior
        prior = scipy.stats.expon(scale = 1) # attention scale = 1/lambda in scipy!
        ax1.fill_between(bins1, -9, -9 + 5*prior.pdf(bins1), color = colors[2], label = "Prior")
        ax2.fill_between(bins2, -9, -9 + 5*prior.pdf(bins2), color = colors[2])


        for i in range(len(mcmc)):
            l = mcmc[i]["lambda"]
            hst = np.histogram(l, bins = bins1, density= True)
            shift = 3*i
            ax1.fill_between(bins1_c, shift , shift + hst[0], alpha = 0.7, color = colors[3], ec = None)
            all_infer[0] += list(l)

            mu = mcmc[i]["mu"]
            hst = np.histogram(mu, bins = bins2, density= True)
            shift = 7*i
            ax2.fill_between(bins2_c, shift , shift + hst[0], alpha = 0.7, color = colors[3], ec = None)
            all_infer[1] += list(mu)

            div = mcmc[i]["lambda"] - mcmc[i]["mu"]
            hst = np.histogram(div, bins = bins3, density= True)
            shift = 3*i
            ax3.fill_between(bins3_c, shift , shift + hst[0], alpha = 0.7, color = colors[3], ec = None)
            all_infer[2] += list(div)

            turnov = mcmc[i]["mu"] / mcmc[i]["lambda"]
            hst = np.histogram(turnov, bins = bins4, density= True)
            shift = 7*i
            ax4.fill_between(bins4_c, shift , shift + hst[0], alpha = 0.7, color = colors[3], ec = None)
            all_infer[3] += list(turnov)
            
        ax1.fill_between([],[],[], alpha = 0.7, color = colors[3], ec = None, label = "Posterior")

        ax1.axvline(est_l, ls = "--", color = colors[1], lw = 2, label = "Estimate")
        ax2.axvline(est_mu, ls = "--", color = colors[1], lw = 2)
        ax3.axvline(est_l - est_mu, ls = "--", color = colors[1], lw = 2)
        ax4.axvline(est_mu/est_l, ls = "--", color = colors[1], lw = 2)

        ax1.set_xlabel("Speciation rate")
        ax2.set_xlabel("Extinction rate")
        ax3.set_xlabel("Diversification rate")
        ax4.set_xlabel("Turnover rate")

        for i in range(4):
            axes[i].set_yticks([])
            axes[i].set_xlim(bins[i][0], bins[i][-1])

            ylim_lock = axes[i].get_ylim()
            axes[i].errorbar([np.mean(all_infer[0])], [ylim_lock[0]],  
                             xerr = [np.std(all_infer[i])], marker = "o", 
                             color = colors[3], capsize = 3, label = "Mean, std of fit", 
                             ls = "")
        
        ax1.legend(loc = "upper right")

        fig.suptitle("$\\lambda_1 = {:.3f}, \\lambda_2 = {}, \\lambda_3 = {}, \\mu_1 = {}, \\mu_2 = {}$".format(par_simul["PBD.l1"], 
                                                                                                        par_simul["PBD.l2"], 
                                                                                                        par_simul["PBD.l3"], 
                                                                                                        par_simul["PBD.mu1"],
                                                                                                        par_simul["PBD.mu2"]));
        plt.savefig("../fig/PBD_bayesian_fit_BD/distrib_posterior/par{}-var{}.png".format(param_vary, i_param_var))
        plt.close()