using Random
using LinearAlgebra
using Statistics

#----Global Constant Parameters Settings
include("global.jl")
#-------------------
#dispatched project-depended files:
#include("utils.jl")
include("MCStatistics.jl")
include("MCModel.jl")
include("method.jl")

for beta in [0.5 2 5 10 15 20 30 40 50 60 80 100]
    print('\nbeta=',beta,)
    print('seed=',random_seed)
    mcml(beta)
    print('\n')
end
# wanglandau(0.0001)

