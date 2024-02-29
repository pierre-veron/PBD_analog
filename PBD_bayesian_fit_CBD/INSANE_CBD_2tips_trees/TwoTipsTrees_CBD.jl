# using Pkg; ]add Tapestree#insane
using Tapestree
using Random: seed!
using Plots

root = 8
tree = sT_label(sT_label(Float64(root), "t1"), sT_label(Float64(root), "t2"), 0.0, "")

## Sample a tree in the MCMC trace
seed_nb = 0
seed!(seed_nb)

### CBD inference ###
veryShortMCMC = false
shortMCMC = false

λ_prior  = (1.0, 1.0)
μ_prior  = (1.0, 1.0)
niter    = veryShortMCMC ? 10 : (shortMCMC ? 100_000 : 50_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? 100 : 5_000)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 10_000 : 50_000)
nflush   = nthin
ofile    = "CBD_2tips_$(niter)iter_rootAge$(root)_seed$(seed_nb)" ; isdir("outputs/") || mkdir("outputs/")
ϵi       = 0.4
λi       = NaN
μi       = NaN
pupdp    = (0.2,0.2,0.2)
survival = true
mxthf    = Inf
tρ       = Dict("" => 1.0)

seed!(seed_nb); out = insane_cbd(tree::sT_label, 
		                         λ_prior   = λ_prior,
		                         μ_prior   = μ_prior,
		                         niter     = niter,
		                         nthin     = nthin,
		                         nburn     = nburn,
		                         nflush    = nflush,
		                         ofile     = "outputs/"*ofile,
		                         ϵi        = ϵi,
		                         λi        = λi,
		                         μi        = μi,
		                         pupdp     = pupdp,
		                         survival  = survival,
		                         mxthf     = mxthf,
		                         tρ        = tρ)

out_trees = out[2]
#out_trees = iread("outputs/$(ofile).txt")

gr(dpi=400, size=(500,300))
anim_tree = @animate for tree_i in out_trees[ Integer.(round.(collect(range(1,length(out_trees), 100))))]
  plot(tree_i, shownodes=(true, true, true), showda=true)
end
isdir("Animations/") || mkdir("Animations/") ; mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

plot(last(ltt.(out_trees), 500), 0.05) ; plot!(ltt(tree), scale=:identity)
isdir("Images/") || mkdir("Images/") ; png("Images/$(ofile)_LTT.png")
