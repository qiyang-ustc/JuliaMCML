function wanglandau(factor::Float64)
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
# set_target!(target)                 # IMPORTANT LINE!!!! If you want to use random target,comment this line.

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

function mcml(Jcp::Float64)
    #-------------------
    target = Target(dim,layers)
    model = Model(dim,layers)
    vectors = create_bit_verctor(dim)
    #---Variables
    #---Statistics
    quantity = zeros(Float64,Nsamp,Nmea)
    observables = zeros(Float64,Nblck,Nobs)
    #include("p_method.jl")
    
    # set_target!(target) # NEVER MISS THIS LINE!!!!!!!!!!
    #p_method
    # This jl file include an easy - optional way to measure p
    p = zeros(Int128,layers+1,layers+1,NV,NV)
    weight = [2^(i-1) for i = 1:1:dim]


function vector2number(v::Array{SpinValueType,1},weight::Array{Int,1})
    return sum(div.((v.+1),2).*weight)+1
end

function number2vector!(v::Array{SpinValueType,1},n::Int)
    temp = n-1
    for j in 1:1:dim
        v[j] = (temp % 2)
        temp >>= 1
    end
    v .*= 2
    v .-= 1
end

function measure_p!(model::Model,target::Target,p::Array{Int128,4},weight::Array{Int,1})
    #TODO: This function has not been tested
    v = zeros(SpinValueType,dim)
    temp = zeros(Int,1+layers)
    for i in 1:1:NV
        v = number2vector!(v,i)
        temp[1] = i
        for j in 1:1:layers
            v = model.A[:,:,j]*v .+ model.B[:,j] .+ target.C[:,j]
            v = SpinValueType.(sign.(v))
            temp[1+j] = vector2number(v,weight)
        end 
        for j in 1:1:layers+1
            for k in 1:1:layers+1
                p[j,k,temp[j],temp[k]]+=1
            end
        end
    end
end

function print_p(temp::Array{Int128,4})
    #TODO: This function is only avaliable for n = 3 test!
        #TODO:This function need to be test
    #---- in order to save memory
    p = zeros(Float64,layers+1,layers+1,NV,NV)
    for i in 1:1:layers+1
        for j in 1:1:layers+1
            # print(temp[i,j,:,:])
            # print(sum(temp[i,j,:,:]))
            p[i,j,:,:] = temp[i,j,:,:]./sum(temp[i,j,:,:])
        end
    end
    s01 = -sum(log.(p[1,2,:,:].+1E-100).*p[1,2,:,:])
    s02 = -sum(log.(p[1,3,:,:].+1E-100).*p[1,3,:,:])
    s03 = -sum(log.(p[1,4,:,:].+1E-100).*p[1,4,:,:])
    s12 = -sum(log.(p[2,3,:,:].+1E-100).*p[2,3,:,:])
    s13 = -sum(log.(p[2,4,:,:].+1E-100).*p[2,4,:,:])
    s23 = -sum(log.(p[3,4,:,:].+1E-100).*p[3,4,:,:])
    s0 = -sum(log.(sum(p[1,2,:,:],dims=2).+1E-100).*sum(p[1,2,:,:],dims=2))
    s1 = -sum(log.(sum(p[1,2,:,:],dims=1).+1E-100).*sum(p[1,2,:,:],dims=1))
    s2 = -sum(log.(sum(p[3,4,:,:],dims=2).+1E-100).*sum(p[3,4,:,:],dims=2))
    s3 = -sum(log.(sum(p[3,4,:,:],dims=1).+1E-100).*sum(p[3,4,:,:],dims=1))

    I01 = s0 + s1 - s01
    I02 = s0 + s2 - s02
    I03 = s0 + s3 - s03
    I12 = s1 + s2 - s12
    I13 = s1 + s3 - s13
    I23 = s2 + s3 - s23

    print(I01)
    print('\t')
    print(I02)
    print('\t')
    print(I03)
    print('\t')
    print(I12)
    print('\t')
    print(I13)
    print('\t')
    print(I23)
    print('\t')
end
#end p_method

    #---Monte Carlo Simulation----
    energy = cal_energy(model,target,vectors)*NC #use NC to normalize energy to [0,1]
    for itoss in 1:1:Ntoss
        for isamp in 1:1:Nsamp
            index = random_index(model)
            picker!(model,index)
            new_energy = cal_energy(model,target,vectors)*NC
            if rand()<exp(Jcp*(energy-new_energy))
                energy = new_energy
            else 
                picker!(model,index)
            end
        end
    end
    
    io = open("temp.dat","w")
    for iblck in 1:1:Nblck
        for isamp in 1:1:Nsamp
            index = random_index(model)
            picker!(model,index)
            new_energy = cal_energy(model,target,vectors)*NC
            if rand()<exp(Jcp*(energy - new_energy)) #accept the move 
                energy = new_energy
            else 
                picker!(model,index)
            end
            quantity[isamp,1] = energy      
            quantity[isamp,2] = energy^2
            quantity[isamp,3] = sum(model.A[:,:,1])
            quantity[isamp,4] = sum(model.A[:,:,2])
            quantity[isamp,5] = sum(model.A[:,:,3])
            quantity[isamp,6] = sum(model.A[:,:,1])^2
            quantity[isamp,7] = sum(model.A[:,:,2])^2
            quantity[isamp,8] = sum(model.A[:,:,3])^2
            measure_p!(model,target,p,weight)
        end
        function normalize(quantity::Array{Float64,2},observables::Array{Float64,2},iblck::Int64)
            for imea in 1:Nmea
                observables[iblck,imea] = mean(quantity[:,imea])
            end
            # observables[iblck,9] = observables[iblck,2] - observables[iblck,1]^2
        end
        # print(mean(quantity[:,1]),'\n')
        normalize(quantity, observables, iblck)
    end
    close(io)
    statistics(observables)
    print_p(p)
end