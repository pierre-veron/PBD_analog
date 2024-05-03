# using Pkg; ]add Tapestree#insane
#using Tapestree

# Needs to use the custom INSANE build to set gamma priors
#include("/data/biodiv/pveron/PBD_analog/PBD_bayesian_fit_BDD/INSANE_BDD/Source_INSANE.jl");
include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Random: seed!
using Plots
using Distributions
using DataFrames

ntrees = 100  # number of trees

## Simulate datasets to assess model calibration
# Parameters
tor = 10.0      # time of origin

# Prior distributions
Gamma15_1  = Gamma(1.5,1)
# Gamma12_02  = Gamma(1.2,0.2)
# Norm01_005 = Normal(-0.1, 0.05)
InvG5_05  = InverseGamma(5.0, 0.05)

# Simulate trees
trees  = iTpbd[]
rtrees = iTpbd[]
pars   = DataFrame(b=Float64[], λ=Float64[], μ=Float64[], αb=Float64[],
                   αλ=Float64[], σb=Float64[], σλ=Float64[], σμ=Float64[])
i = 1
s = 0
count_limit = 0
while i <= ntrees
    global i
    global s
    global count_limit
    # Draw parameters from prior distributions
    @show s
    seed!(s) ; s+=1
    # b0 = rand(Gamma12_05)
    b0 = 1.0
    λ0 = rand(Gamma15_1)
    # μ0 = rand(Gamma12_02)
    μ0 = 0.7
    # αb = rand(Norm01_005)
    # αλ = rand(Norm01_005)
    αb = αλ = -0.05
    # α  = 0.0
    σb = sqrt(rand(InvG5_05))
    σλ = sqrt(rand(InvG5_05))
    σμ = sqrt(rand(InvG5_05))
    
    parms = round.((b0,λ0,μ0,αb,αλ,σb,σλ,σμ); sigdigits=3)
    @show parms

    nlim = 20_000
    tree = sim_gbmpbd(tor, b0=b0, λ0=λ0, μ0=μ0, αb=αb, αλ=αλ, σb=σb, σλ=σλ, σμ=σμ, nlim=nlim)
    if nnodes(tree) >= nlim
        @warn "Too many nodes ($(nnodes(tree)))"
        count_limit += 1
        continue
    end
    tree_rmi = remove_incipient(iTpbd(tree))

    # Condition on the survival of both crown lineages
    if ntipsalive(tree_rmi)>1
        if ntipsalive(tree_rmi.d1)>0 && ntipsalive(tree_rmi.d2)>0
            if ntips(tree_rmi)>1000
                @warn "Too many tips ($(ntips(tree_rmi)))"
                count_limit += 1
            else 
                @show tree_rmi
                push!(trees, tree)
                recTree = remove_extinct(tree_rmi)
                fixtree!(recTree)
                push!(rtrees, recTree)
                push!(pars, (b0, λ0, μ0, αb,αλ, σb, σλ, σμ))
                i += 1
            end
        end
    end
    count_limit<100 || break
end
count_limit
trees
rtrees

# Save trees and true parameters
for i in Base.OneTo(ntrees)
    @show i
    path = "/Volumes/data/PBD_analog/PBDD_bayesian_fit_BDD/datasets_PBDD/dataset$i/"
    isdir(path) || mkpath(path)

    iwrite(trees[i], path*"tree_full")
    iwrite(rtrees[i], path*"tree_reconstructed")
    write_newick(trees[i], path*"tree_full")
    write_newick(rtrees[i], path*"tree_reconstructed")

    open(path*"params.csv", "w") do io
        writedlm(io, [names(pars)])
        writedlm(io, [pars[i,:]])
    end
end



