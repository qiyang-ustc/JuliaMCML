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

function picker!(model::Model)
    A = model.A
    B = model.B
    p=rand(UInt16)%(length(A)+length(B))
    if p > length(A)
        B[p-length(A)] = -B[p-length(A)]
    else 
        A[p] = -A[p]
    end
end

function output(model::Model,target::Target,vectors::Array{SpinValueType,2})
    temp = vectors
    for i in 1:1:3
        temp = model.A[:,:,i]*temp + model.B[:,i] + target.C[:,i]
        temp = sign.(temp)
    end
    return temp
end

function output(target::Target,vectors::Array{SpinValueType,2})
    temp = vectors
    for i in 1:1:3
        temp = target.A[:,:,i]*temp + target.B[:,i] + target.C[:,i]
    end 
end

function energy(model::Model,target::Target,vectors::Array{SpinValueType,2})
    output = output(model,target,vectors)
    target_output = output(target,vectors)
    return sum(abs.(output-target_output))/2
end