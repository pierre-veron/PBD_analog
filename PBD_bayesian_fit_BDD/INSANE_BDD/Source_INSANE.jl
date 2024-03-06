#=

source file for testing INSANE locally

Jérémy Andréoletti

V(°-°V)

Created 28 09 2021
=#

using Random: randexp, randn!, shuffle!
using SpecialFunctions: loggamma
using SpecialFunctions: erf
using DelimitedFiles: writedlm
using ProgressMeter: Progress, next!
using Statistics: quantile, mean, median
using LoopVectorization: @turbo
using PlotUtils: cgrad, palette
using RecipesBase
using Parsers: parse as Pparse
using Distributions: Poisson, Uniform

const accerr = √eps()
const epochs = Float64[720.0, 635.0, 538.8, 521.0, 509.0, 497.0, 485.4, 470.0, 458.4, 443.8, 443.8, 419.2, 393.3, 382.7, 358.9, 323.2, 298.9, 273.01, 259.51, 251.902, 247.2, 237, 201.3, 174.1, 163.5, 145.0, 100.5, 66.0, 56.0, 33.9, 23.03, 5.333, 2.58]

Tapestree_path = "/users/biodiv/pveron/PBD_analog/PBD_bayesian_fit_BDD/Tapestree.jl/"

include(Tapestree_path*"src/Utils.jl")
include(Tapestree_path*"src/utils/density_functions.jl")

files = readdir(Tapestree_path*"src/insane/"; join=true);
files = [Tapestree_path*m.match for m in match.(r"src/insane/.*jl", files) if m != nothing]
append!(files, readdir(Tapestree_path*"src/utils/"; join=true))
include(Tapestree_path*"src/insane/iTree.jl")
include(Tapestree_path*"src/insane/iTreeX.jl")
include(Tapestree_path*"src/insane/iB.jl")
include.(files)

const iTd = Dict{String, DataType}("sTpb"  => sTpb,
                                   "sTbd"  => sTbd,
                                   "sTfbd" => sTfbd,
                                   "iTpb"  => iTpb,
                                   "iTce"  => iTce,
                                   "iTct"  => iTct,
                                   "iTbd"  => iTbd,
                                   "iTfbd" => iTfbd)