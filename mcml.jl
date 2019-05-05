using Random
using LinearAlgebra
abspath("D:\\IOPHY\\jmc\\")
#----Global Constant Parameters Settings
Random.seed!(12345)
const SpinValueType = Int8
const layers = 3
const dim = 3
const Nsamp = 1000
const Nblck = 1024
function create_bit_verctor(dim::Int)
    temp_array = Array{SpinValueType}(UndefInitializer(),dim,2^dim)
    for i in 1:1:2^dim
        temp = i
        for j in 1:1:dim     
            temp_array[j,i] = (temp % 2)
            temp >>= 1
        end
    end
    return 2*temp_array.-1
end
vectors = create_bit_verctor(dim)

 #----Global Variables

 #number of A,B,Vectors, and normailzation coefficient
const NA = layers*dim^2 
const NB = layers*dim
const NV = 2^dim
const NC = 1/(NV*dim)
const Jcp = collect(Float64,0:1:1)
#-------------------
#dispatched project-depended files:
include("utils.jl")
include("MCStatistics.jl")
include("MCModel.jl")
#-------------------
target = Target(dim,layers)
model = Model(dim,layers)
#---Variables
#---Statistics
quantity = zeros(Float64,Nsamp)
observables = zeros(Float64,Nblck)
energy_spectrum = Spectrum(2^dim*dim)
#---Monte Carlo Simulation----
for iblck in 1:1:Nblck
    for isamp in 1:1:Nsamp
                

    end
end

