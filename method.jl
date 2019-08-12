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
    
    set_target!(target) # NEVER MISS THIS LINE!!!!!!!!!! 
    #set_target(target): this line will call function set_target!(target::Target) in MCModel.jl
    #So that the target will be set by configuration in target.jl

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

function print_p_for_tempering(temp::Array{Int128,4})
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

function parallel_tempering(Jcp::Array{Float64,1})
    #-------------------
    n_system = size(Jcp)[1] 
    #We need a common target
    target = Target(dim,layers)
    #We need n_system different system/model
    model = [Model(dim,layers) for i in 1:1:n_system]
    vectors = create_bit_verctor(dim)
    #---Statistics
    sample_weight = zeros(Int,Nsamp) #This quantity is used for reweight measuremt
    quantity = zeros(Float64,Nsamp,Nmea,n_system)
    observables = zeros(Float64,Nblck,Nobs,n_system)
    #include("p_method.jl")
    #set_target!(target) # NEVER MISS THIS LINE!!!!!!!!!! 
    #set_target(target): this line will call function set_target!(target::Target) in MCModel.jl
    #So that the target will be set by configuration in target.jl

    #p_method
    # This jl file include an easy - optional way to measure p
    p = zeros(Int128,layers+1,layers+1,NV,NV,n_system)
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
    function measure_p!(model::Model,target::Target,p::Array{Int128,5},imodel::Int,weight::Array{Int,1})
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
                    p[j,k,temp[j],temp[k],imodel]+=1
                end
            end
        end
    end

    
    #end settings
    #---Monte Carlo Simulation----
    #---Initialize energy and thermalization----
    energy = [cal_energy(model[imodel],target,vectors)*NC for imodel in 1:1:n_system]#use NC to normalize energy to [0,1]
    for imodel in 1:1:n_system
        # normal updates
        for itoss in 1:1:Ntoss
            for isamp in 1:1:Nsamp
                index = random_index(model[imodel])
                picker!(model[imodel],index)
                new_energy = cal_energy(model[imodel],target,vectors)*NC
                if rand()<exp(Jcp[imodel]*(energy[imodel]-new_energy))
                    energy[imodel] = new_energy
                else 
                    picker!(model[imodel],index)
                end
            end
        end
        #parallel tempering
        #we do not need parallel tempering in thermalization
    end
    # io = open("temp.dat","w")
    temp_sum = 0
    for iblck in 1:1:Nblck
        for isamp in 1:1:Nsamp
            for imodel in 1:1:n_system
                index = random_index(model[imodel])
                picker!(model[imodel],index)
                new_energy = cal_energy(model[imodel],target,vectors)*NC
                if rand()<exp(Jcp[imodel]*(energy[imodel] - new_energy)) #accept the move 
                    energy[imodel] = new_energy
                    if imodel == n_system
                        temp_sum = temp_sum + 1
                    end
                else 
                    picker!(model[imodel],index)
                end
                quantity[isamp,1,imodel] = energy[imodel]      
                quantity[isamp,2,imodel] = energy[imodel]^2
                quantity[isamp,3,imodel] = sum(model[imodel].A[:,:,1])
                quantity[isamp,4,imodel] = sum(model[imodel].A[:,:,2])
                quantity[isamp,5,imodel] = sum(model[imodel].A[:,:,3])
                quantity[isamp,6,imodel] = sum(model[imodel].A[:,:,1])^2
                quantity[isamp,7,imodel] = sum(model[imodel].A[:,:,2])^2
                quantity[isamp,8,imodel] = sum(model[imodel].A[:,:,3])^2
                measure_p!(model[imodel],target,p,imodel,weight) #imodel modified
            end
            #tempering
            for imodel in 1:1:(n_system-1)
                if rand()<exp((Jcp[imodel]- Jcp[imodel+1])*(energy[imodel]-energy[imodel+1]))
                    temp_model = model[imodel]
                    model[imodel] = model[imodel+1]
                    model[imodel+1] = temp_model
                    #----
                    temp_energy = energy[imodel] 
                    energy[imodel] = energy[imodel+1]
                    energy[imodel+1] = temp_energy
                end
            end
        end
        function normalize!(quantity::Array{Float64,3},observables::Array{Float64,3},iblck::Int64) 
            for imodel in n_system
                for imea in 1:Nmea
                    # print(observables[iblck,imea,imodel],'\n')
                    observables[iblck,imea,imodel] = mean(quantity[:,imea,imodel])
                end
            end
                # observables[iblck,9] = observables[iblck,2] - observables[iblck,1]^2
        end
        # normalize!(quantity, observables, iblck)
        for imodel in 1:n_system
            for imea in 1:Nmea
                observables[iblck,imea,imodel] = mean(quantity[:,imea,imodel])
            end
        end
    end
    # close(io)
    print(temp_sum/Nblck/Nsamp,'\n')
    for i in 1:1:n_system
        print("beta=",Jcp[i],' ',"seed=",random_seed,'\n')
        statistics(observables[:,:,i])
        print_p_for_tempering(p[:,:,:,:,i])
        print("\n\n")
    end
end

function parallel_tempering_with_reweighting(Jcp::Array{Float64,1})
    #-------------------
    n_system = size(Jcp)[1] 
    #We need a common target
    target = Target(dim,layers)
    #We need n_system different system/model
    model = [Model(dim,layers) for i in 1:1:n_system]
    vectors = create_bit_verctor(dim)
    #---Statistics
    sample_weight = zeros(Int,Nsamp,n_system) #This quantity is used for reweight measuremt
    quantity = zeros(Float64,Nsamp,Nmea,n_system)
    observables = zeros(Float64,Nblck,Nobs,n_system)
    #include("p_method.jl")
    #set_target!(target) # NEVER MISS THIS LINE!!!!!!!!!! 
    #set_target(target): this line will call function set_target!(target::Target) in MCModel.jl
    #So that the target will be set by configuration in target.jl

    #p_method
    # This jl file include an easy - optional way to measure p
    p = zeros(Int128,layers+1,layers+1,NV,NV,n_system)
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
    function measure_p!(model::Model,target::Target,p::Array{Int128,5},imodel::Int,weight::Array{Int,1},counts::Int)
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
                    p[j,k,temp[j],temp[k],imodel]+=counts
                end
            end
        end
    end

    
    #end settings
    #---Monte Carlo Simulation----
    #---Initialize energy and thermalization----
    energy = [cal_energy(model[imodel],target,vectors)*NC for imodel in 1:1:n_system]#use NC to normalize energy to [0,1]
    for imodel in 1:1:n_system
        # normal updates
        for itoss in 1:1:Ntoss
            for isamp in 1:1:Nsamp
                index = random_index(model[imodel])
                picker!(model[imodel],index)
                new_energy = cal_energy(model[imodel],target,vectors)*NC
                if rand()<exp(Jcp[imodel]*(energy[imodel]-new_energy))
                    energy[imodel] = new_energy
                else 
                    picker!(model[imodel],index)
                end
            end
        end
        #parallel tempering
        #we do not need parallel tempering in thermalization
    end
    # io = open("temp.dat","w")
    for iblck in 1:1:Nblck
        for isamp in 1:1:Nsamp
            for imodel in 1:1:n_system
                index = random_index(model[imodel])
                picker!(model[imodel],index)
                new_energy = cal_energy(model[imodel],target,vectors)*NC
                p_update = exp(Jcp[imodel]*(energy[imodel]-new_energy))
                picker!(model[imodel],index)
                # print(p_update,'\n')
                if p_update < 1.0
                    temp_weight = ceil(log(rand())/log(1-p_update-EPSILON))
                    sample_weight[isamp,imodel] = min(temp_weight,TOL_STEPS)
                    #---measure
                    quantity[isamp,1,imodel] = energy[imodel]      
                    quantity[isamp,2,imodel] = energy[imodel]^2
                    quantity[isamp,3,imodel] = sum(model[imodel].A[:,:,1])
                    quantity[isamp,4,imodel] = sum(model[imodel].A[:,:,2])
                    quantity[isamp,5,imodel] = sum(model[imodel].A[:,:,3])
                    quantity[isamp,6,imodel] = sum(model[imodel].A[:,:,1])^2
                    quantity[isamp,7,imodel] = sum(model[imodel].A[:,:,2])^2
                    quantity[isamp,8,imodel] = sum(model[imodel].A[:,:,3])^2
                    #---measure
                    measure_p!(model[imodel],target,p,imodel,weight,sample_weight[isamp,imodel]) #imodel modified
                     # if sample_weight> TOL_STEPS the move accept, else reject
                    if temp_weight <= TOL_STEPS
                        picker!(model[imodel],index)
                        energy[imodel] = new_energy
                    end
                else
                    sample_weight[isamp,imodel] = 1
                    #---measure---
                    quantity[isamp,1,imodel] = energy[imodel]      
                    quantity[isamp,2,imodel] = energy[imodel]^2
                    quantity[isamp,3,imodel] = sum(model[imodel].A[:,:,1])
                    quantity[isamp,4,imodel] = sum(model[imodel].A[:,:,2])
                    quantity[isamp,5,imodel] = sum(model[imodel].A[:,:,3])
                    quantity[isamp,6,imodel] = sum(model[imodel].A[:,:,1])^2
                    quantity[isamp,7,imodel] = sum(model[imodel].A[:,:,2])^2
                    quantity[isamp,8,imodel] = sum(model[imodel].A[:,:,3])^2
                    #---measure---
                    measure_p!(model[imodel],target,p,imodel,weight,sample_weight[isamp,imodel]) #imodel modified
                    picker!(model[imodel],index)       
                    energy[imodel] = new_energy             
                end
                # if rand()<exp(Jcp[imodel]*(energy[imodel] - new_energy)) #accept the move 
                #     energy[imodel] = new_energy
                # else 
                #     picker!(model[imodel],index)
                # end
            end
            #tempering
            for imodel in 1:1:(n_system-1)
                if rand()<exp((Jcp[imodel]- Jcp[imodel+1])*(energy[imodel]-energy[imodel+1]))
                    temp_model = model[imodel]
                    model[imodel] = model[imodel+1]
                    model[imodel+1] = temp_model
                    #----
                    temp_energy = energy[imodel] 
                    energy[imodel] = energy[imodel+1]
                    energy[imodel+1] = temp_energy
                end
            end
        end
        # normalize!(quantity, observables, iblck)
        # print(sample_weight[:,16])
        for imodel in 1:n_system
            for imea in 1:Nmea
                observables[iblck,imea,imodel] = sum(sample_weight[:,imodel].*quantity[:,imea,imodel])/sum(sample_weight[:,imodel])
            end
        end
    end
    # close(io)
    for i in 1:1:n_system
        print("beta=",Jcp[i],' ',"seed=",random_seed,'\n')
        statistics(observables[:,:,i])
        print_p_for_tempering(p[:,:,:,:,i])
        print("\n\n")
    end
end