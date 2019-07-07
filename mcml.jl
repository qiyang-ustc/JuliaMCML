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

# mcml(30.0)
wanglandau(0.0001)

