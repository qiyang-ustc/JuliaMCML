using Random
using LinearAlgebra

#----Global Constant Parameters Settings
include("global.jl")
#-------------------
#dispatched project-depended files:
#include("utils.jl")
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
set_target!(target)
#TODO: = []

#---Monte Carlo Simulation----
energy = Int(cal_energy(model,target,vectors))
for iblck in 1:1:Nblck
    for isamp in 1:1:Nsamp
        index = random_index(model)
        picker!(model,index)
        new_energy =  Int(cal_energy(model,target,vectors))
        if rand()<exp(get_factor(energy_spectrum,energy)-get_factor(energy_spectrum,new_energy)) #accept the move 
            energy = new_energy
        else 
            picker!(model,index)
        end
        add_factor!(energy_spectrum,energy,factor)
    end
end
process_spectrum(energy_spectrum,-100.0)
normalize_factor(energy_spectrum)

energy = collect(0.0:1.0:2.0^dim*dim)

for i in 1:1:length(energy_spectrum.factor)
   print(energy[i],"    ",exp.(energy_spectrum.factor[i]),"\n") 
end

function check_converge(model::Model,target::Target,energy_spectrum::Spectrum,vectors)
    temp_Spectrum = deepcopy(energy_spectrum)
    for i in 1:1:temp_Spectrum.length
        temp_Spectrum.factor[i] = 0
    end
    energy = Int(cal_energy(model,target,vectors))
    for i in 1:1:LENGTH_OF_CONVERGENCE_CHECK
        index = random_index(model)
        picker!(model,index)
        new_energy =  Int(cal_energy(model,target,vectors))
        if rand()<exp(get_factor(energy_spectrum,energy)-get_factor(energy_spectrum,new_energy)) #accept the move 
            energy = new_energy
            add_factor!(temp_Spectrum,energy,1.0)
        else 
            add_factor!(temp_Spectrum,energy,1.0)
            picker!(model,index)
        end
    end
    return temp_Spectrum.factor/maximum(temp_Spectrum.factor)
end
@show check_converge(model,target,energy_spectrum,vectors)
end
mcml(0.0000001)