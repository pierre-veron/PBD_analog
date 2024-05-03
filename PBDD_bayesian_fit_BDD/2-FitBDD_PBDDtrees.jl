# using Pkg; ]add Tapestree#insane
#using Tapestree

# Needs to use the custom INSANE build to set gamma priors
#include("/data/biodiv/pveron/PBD_analog/PBD_bayesian_fit_BDD/INSANE_BDD/Source_INSANE.jl");
include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Random: seed!
using Plots
using DelimitedFiles

## Read in data
i = ARGS[1]
seed_nb = ARGS[2]
path = "/Volumes/data/PBD_analog/PBDD_bayesian_fit_BDD/"
# tree = read_newick("/data/biodiv/pveron/PBD_analog/trees12152/$(tree_name)")
# tree = read_newick(path*"datasets_PBDD/dataset$(i)/tree_reconstructed.tre")
# # tree = read_newick("/Volumes/data/PBD_analog/PBDD_bayesian_fit_BDD/datasets_PBDD/dataset$(i)/tree_reconstructed.tre")

# ## Sample a tree in the MCMC trace
# # seed_nb = rand(1:10000)
# seed_nb = parse(Float64, i)
# seed!(seed_nb)

# ### BDD inference ###

veryShortMCMC = false
shortMCMC = false

# λa_prior = (1.5, 1.0)
# μa_prior = (1.5, 1.0)
# α_prior  = (0.0, 0.5)
# σλ_prior = (5.0, 0.5)
# σμ_prior = (5.0, 0.5)
niter    = veryShortMCMC ? 10 : (shortMCMC ? 50_000 : 2_000_000)
# nthin    = veryShortMCMC ? 1 : (shortMCMC ? niter+1 : niter+1)
# nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 100_000)
# nflushθ  = Int64(ceil(niter/20_000))
# nflushΞ  = Int64(ceil(niter/100))
ofile    = "dataset$(i)-FitBDD_PBDDtrees_$(niter)iter_seed$(seed_nb)"
# isdir("outputs/") || mkdir("outputs/")
# ϵi       = 0.2
# λi       = NaN
# μi       = NaN
# αi       = 0.0
# σλi      = 0.1
# σμi      = 0.1
# pupdp    = (0.02, 0.1, 0.01, 0.1, 0.2)
# # pupdp    = (0.0, 0.1, 0.01, 0.1, 0.2)
# δt       = 1e-3
# survival = true
# mxthf    = Inf
# prints   = 5
# stnλ     = 0.5
# stnμ     = 0.5
# tρ       = Dict("" => 1.0)


# seed!(seed_nb); insane_gbmbd(tree::sT_label,
#                              λa_prior = λa_prior,
#                              μa_prior = μa_prior,
#                              α_prior  = α_prior,
#                              σλ_prior = σλ_prior,
#                              σμ_prior = σμ_prior,
#                              niter    = niter,
#                              nthin    = nthin,
#                              nburn    = nburn,
#                              nflushθ  = nflushθ,
#                              nflushΞ  = nflushΞ,
#                              ofile    = "outputs/"*ofile,
#                              ϵi       = ϵi,
#                              λi       = λi,
#                              μi       = μi,
#                              αi       = αi,
#                              σλi      = σλi,
#                              σμi      = σμi,
#                              pupdp    = pupdp,
#                              δt       = δt,
#                              survival = survival,
#                              mxthf    = mxthf,
#                              prints   = prints,
#                              stnλ     = stnλ,
#                              stnμ     = stnμ,
#                              tρ       = tρ)

out_trees = iread(path*"outputs/$(ofile).txt")
# tree_PBDD = iread(path*"datasets_PBDD/dataset$(i)/tree_full.txt")[1]
rtree_PBDD = iread(path*"datasets_PBDD/dataset$(i)/tree_reconstructed.txt")[1]

# ### Plot trees and rates ###

# ENV["GKSwstype"] = "nul"
# # isdir(path*"Animations/") || mkdir(path*"Animations/")
# isdir(path*"Images/") || mkdir(path*"Images/")
# gr(dpi=400, size=(500,300))
# # anim_tree = @animate for tree_i in out_trees
# #  plot(tree_i, shownodes=(true, true, true), showda=true)
# # end
# # mp4(anim_tree, path*"Animations/$(ofile)_anim_tree.mp4", fps=5)

# # gr(dpi=400, size=(500,300))
# # anim_tree = @animate for tree_i in out_trees
# #  plot(tree_i, b, shownodes=(true, true, true))
# # end
# # mp4(anim_tree, path*"Animations/$(ofile)_anim_tree_λ.mp4", fps=5)

# plot(ltt.(out_trees), 0.05) ; plot!(ltt(tree), scale=:identity)
# png(path*"Images/$(ofile)_LTT.png")

# plot(tree_PBDD)
# png(path*"Images/$(ofile)_PBDDtree.png")
# plot(tree_PBDD, b, ylab="Initiation rate", left_margin=7Plots.mm)
# png(path*"Images/$(ofile)_PBDDtree_b.png")
# plot(tree_PBDD, c, ylab="Completion rate", left_margin=7Plots.mm)
# png(path*"Images/$(ofile)_PBDDtree_c.png")
# plot(tree_PBDD, d, ylab="Death rate", left_margin=7Plots.mm)
# png(path*"Images/$(ofile)_PBDDtree_d.png")

# plot(tree_PBDD, b, δt, q0=Float64[], ylab="b(t)", linecolor="darkblue", fillcolor="skyblue", label=["50% IQR" "Mean b(t)"])
# plot!(tree_PBDD, c, δt, q0=Float64[], ylab="λ(t)", linecolor="green", fillcolor="lightgreen", label=["50% IQR" "Mean c(t)"])
# plot!(out_trees, b, δt, ylab="BDD λ(t)", label=["" "" ""])
# plot!(legend=:outerright, legendcolumns=1, legendtitle="PBDD rates", legendtitlefontsize=8)
# png(path*"Images/$(ofile)_BDDλTT_PBDDbTTλTT.png")

# plot(tree_PBDD, d, δt, q0=Float64[], linecolor="darkred", fillcolor="pink", label=["50% IQR" "Mean μ(t)"])
# plot!(out_trees, d, δt, ylab="BDD μ(t)", label=["" "" ""])
# plot!(legend=:outerright, legendcolumns=1, legendtitle="PBDD rate", legendtitlefontsize=8)
# png(path*"Images/$(ofile)_BDDμTT_PBDDμTT.png")

### Save tip rates ###

isdir(path*"Rates") || mkdir(path*"Rates")
writedlm(path*"Rates/dataset$(i)_tipb.csv", tip_rates(rtree_PBDD, b))
writedlm(path*"Rates/dataset$(i)_tipc.csv", tip_rates(rtree_PBDD, c))
writedlm(path*"Rates/dataset$(i)_tipd.csv", tip_rates(rtree_PBDD, d))
writedlm(path*"Rates/$(ofile)_tipλ.csv", tip_rates.(remove_extinct.(out_trees), b))
writedlm(path*"Rates/$(ofile)_tipμ.csv", tip_rates.(remove_extinct.(out_trees), d))




