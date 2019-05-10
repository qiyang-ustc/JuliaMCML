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
const CUT_ENERGY = 100  # the value in WangLandau for cut off low degeneracy state
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
function mcml(factor::Float64)
#-------------------
target = Target(dim,layers)
model = Model(dim,layers)
vectors = create_bit_verctor(dim)
#---Variables
#---Statistics
quantity = zeros(Float64,Nsamp)
observables = zeros(Float64,Nblck)
energy_spectrum = Spectrum(2^dim*dim+1)
set_ground!(energy_spectrum,0)
#---Monte Carlo Simulation----
energy = Int(cal_energy(model,target,vectors))
for iblck in 1:1:Nblck
    for isamp in 1:1:Nsamp
        index = random_index(model)
        picker!(model,index)
        new_energy =  Int(cal_energy(model,target,vectors))
        if rand()<exp(get_factor(energy_spectrum,energy)-get_factor(energy_spectrum,new_energy)) #accept the move 
            add_factor!(energy_spectrum,new_energy,factor)
            energy = new_energy
        else 
            picker!(model,index)
            add_factor!(energy_spectrum,energy,factor)
        end
    end
end
process_spectrum(energy_spectrum,-100.0)
normalize_factor(energy_spectrum)
@show exp.(energy_spectrum.factor)

end

mcml(0.01)