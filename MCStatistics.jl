mutable struct Spectrum
    length::Int
    count ::Array{Int128,1}
    factor ::Array{Float64,1}
    level ::Array{Int,1}
    ground ::Int
    function Spectrum(length::Int)
        count = zeros(Int128,length)
        factor = zeros(Float64,length)
        level = zeros(Float64,length)
        ground = 0
        new(length,count,factor,level,ground)
    end
end

function set_ground!(spectrum::Spectrum,g::Int)
    spectrum.ground = g-1
    return 0
end
function add_count!(spectrum::Spectrum,energy::Int)
    spectrum.count[-spectrum.ground+energy]+=1
    return 0
end

function get_count(spectrum::Spectrum,energy::Int)
    return spectrum.count[-spectrum.ground+energy]
end
function add_factor!(spectrum::Spectrum,energy::Int,temp::Float64)
    return spectrum.factor[-spectrum.ground+energy]+=temp
end 

function get_factor(spectrum::Spectrum,energy::Int)
    return spectrum.factor[-spectrum.ground+energy]
end 

function normalize_factor(spectrum::Spectrum)
    max_factor = maximum(spectrum.factor)
    spectrum.factor = spectrum.factor .- max_factor
    spectrum.factor = exp.(spectrum.factor)
    spectrum.factor = log.(spectrum.factor/sum(spectrum.factor)) #normalization
    beta = 1.0

    return 1
end

function process_spectrum(spectrum::Spectrum,ground_value::Real)
    max_factor = maximum(spectrum.factor)
    for i in 1:spectrum.length
        if spectrum.factor[i]<(ground_value+max_factor)
            spectrum.factor[i]-=Inf
        end
    end
    return 1
end

function statistics(observables::Array{Float64,2})
#observables = zeros(Float64,Nblck,Nobs)
    for iobs in 1:1:Nobs
        ave = mean(observables[:,iobs])
        sva = sqrt(var(observables[:,iobs])/(Nblck-1))
        col = cov(observables[1:Nblck-1,iobs],observables[2:Nblck,iobs])/var(observables[:,iobs])
        print(iobs,'\t',ave,'\t',sva,'\t',col,'\n')
    end
end