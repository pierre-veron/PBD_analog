# using Pkg; ]add Tapestree#insane
using Tapestree
using Random: seed!
using Plots

## Read in data
tree_nb = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])
tree = read_newick("", true)

## Sample a tree in the MCMC trace
seed_nb = 0
seed!(seed_nb)

### BDD inference ###

veryShortMCMC = false
shortMCMC = false

λa_prior = (1.0, 1.0)
μa_prior = (1.0, 1.0)
α_prior  = (0.0, 1.0)
σλ_prior = (3.0, 0.5)
σμ_prior = (3.0, 0.5)
niter    = veryShortMCMC ? 10 : (shortMCMC ? 50_000 : 20_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? 10 : 10_000)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 1_000_000)
nflush   = nthin
ofile    = "BDD_Cetaceans_$(niter)iter_seed$(seed_nb)"
ϵi       = 0.2
λi       = NaN
μi       = NaN
αi       = 0.0
σλi      = 0.1
σμi      = 0.1
pupdp    = (0.01, 0.01, 0.1, 0.2)
δt       = 1e-3
survival = true
mxthf    = Inf
prints   = 5
tρ       = Dict("" => 1.0)


seed!(seed_nb); insane_gbmbd(tree::sT_label,
                             λa_prior = λa_prior,
                             μa_prior = μa_prior,
                             α_prior  = α_prior,
                             σλ_prior = σλ_prior,
                             σμ_prior = σμ_prior,
                             niter    = niter,
                             nthin    = nthin,
                             nburn    = nburn,
                             nflush   = nflush,
                             ofile    = "outputs/"*ofile,
                             ϵi       = ϵi,
                             λi       = λi,
                             μi       = μi,
                             αi       = αi,
                             σλi      = σλi,
                             σμi      = σμi,
                             pupdp    = pupdp,
                             δt       = δt,
                             survival = survival,
                             mxthf    = mxthf,
                             prints   = prints,
                             tρ       = tρ)

out_trees = iread("outputs/$(ofile).txt")[1:Int64(niter/nthin)]
tree_nb = findall(x->x==treelength(remove_extinct(out_trees[1])), treelength.(sT_label.(remove_fossils.(trees_rec))))[1]
tree = sT_label(remove_fossils(trees_rec[tree_nb]))

ENV["GKSwstype"] = "nul"
gr(dpi=400, size=(500,300))
anim_tree = @animate for tree_i in out_trees[ Integer.(round.(collect(range(1,length(out_trees), 100))))]
  plot(tree_i, shownodes=true, tip=true, speciation=true)
end
mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

function explλ(x) return(exp.(lλ(x))) end
function explμ(x) return(exp.(lμ(x))) end
gr(dpi=400, size=(500,300))
anim_tree = @animate for tree_i in out_trees[ Integer.(round.(collect(range(1,length(out_trees), 100))))]
  plot(tree_i, explλ, clim=(0,1.0), shownodes=true, tip=true, speciation=true)
end
mp4(anim_tree, "Animations/$(ofile)_anim_tree_λ.mp4", fps=5)

plot(rand(out_trees, 500), explλ, δt, ylab="λ(t)")
png("Images/$(ofile)_λTT.png")

plot(rand(out_trees, 500), explμ, δt, ylab="μ(t)")
png("Images/$(ofile)_μTT.png")

plot(ltt.(rand(out_trees, 1000)), 0.05) ; plot!(ltt(tree), scale=:identity)
png("Images/$(ofile)_LTT.png")



