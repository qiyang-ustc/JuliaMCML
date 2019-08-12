using Random
using LinearAlgebra
using Statistics
using DelimitedFiles

#----Global Constant Parameters Settings
include("global.jl")
#-------------------
#dispatched project-depended files:
#include("utils.jl")
include("MCStatistics.jl")
include("MCModel.jl")
include("method.jl")

# print("target seed=",random_seed,'\n')
parallel_tempering_with_reweighting([0.5 0.75 1.0 1.5 2.0 3.5 5 10 15 20 30 40 50 60 80 100][:])

# for beta in [0.5 0.75 1.0 1.5 2.0 3.5 5 10 15 20 30 40 50 60 80 100]
#     Random.seed!(random_seed)
#     print("beta=",beta,)
#     print(" seed=",random_seed,'\n')
#     mcml(beta)
#     print("\n\n")
# end

# wanglandau(0.0001)
