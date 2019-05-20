struct Target
    A::Array{SpinValueType,3}
    B::Array{SpinValueType,2}
    C::Array{Float64,2}
    function Target(dim::Int,layers::Int)
        A = SpinValueType.(rand([-1,1],(dim, dim, layers)))
        B = SpinValueType.(rand([-1,1],(dim, layers)))
        C = Float64.(rand([-0.5,0.5],(dim, layers)))    
        new(A,B,C)
    end
end

struct Model
    A::Array{SpinValueType,3}
    B::Array{SpinValueType,2}
    function Model(dim::Int,layers::Int)
        A = SpinValueType.(rand([-1,1],(dim, dim, layers)))
        B = SpinValueType.(rand([-1,1],(dim, layers)))
        new(A,B)
    end
end

function picker!(model::Model,index::Int)
    A = model.A
    B = model.B
    if index > length(A)
        B[index-length(A)] = -B[index-length(A)]
    else 
        A[index] = -A[index]
    end
end

function random_index(model::Model)
    return Int(rand(UInt128)%(length(model.A)+length(model.B)))+1
end

function output(model::Model,target::Target,vectors::Array{SpinValueType,2})
    temp = vectors
    for i in 1:1:layers
        temp = model.A[:,:,i]*temp .+ model.B[:,i] .+ target.C[:,i]
        temp = sign.(temp)
    end
    return temp
end

function output(target::Target,vectors::Array{SpinValueType,2})
    temp = vectors
    for i in 1:1:layers
        temp = target.A[:,:,i]*temp .+ target.B[:,i] .+ target.C[:,i]
        temp = sign.(temp)
    end 
    return temp
end

function cal_energy(model::Model,target::Target,vectors::Array{SpinValueType,2})
    model_output = output(model,target,vectors)
    target_output = output(target,vectors)
    return sum(abs.(model_output-target_output))/2
end

function create_bit_verctor(dim::Int)
    temp_array = Array{SpinValueType}(UndefInitializer(),dim,2^dim)
    for i in 1:1:2^dim
        temp = i
        for j in 1:1:dim     
            temp_array[j,i] = (temp % 2)
            temp >>= 1
        end
    end
    return SpinValueType.(2*temp_array.-1)
end

function set_target!(target::Target)
    include(".\\target.jl")
    for i in 1:1:layers
        for j in 1:1:dim
            target.B[j,i] = B[j,i]
            target.C[j,i] = C[j,i]
        end
    end
    for i in 1:1:layers
        for j in 1:1:dim
            for k in 1:1:dim
                target.A[j,k,i] = A[j,k,i]
            end
        end
    end
end


function set_target!(target::Target,A::Array{SpinValueType,3},B::Array{SpinValueType,2},C::Array{Float64,2})
    for i in 1:1:layers
        for j in 1:1:dim
            target.B[j,i] = B[j,i]
            target.C[j,i] = C[j,i]
        end
    end
    for i in 1:1:layers
        for j in 1:1:dim
            for k in 1:1:dim
                target.A[j,k,i] = A[j,k,i]
            end
        end
    end
end
