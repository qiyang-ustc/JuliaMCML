# This jl file include an easy - optional way to measure p
p = zeros(Int128,layers+1,layers+1,NV,NV)
weight = [2^(i-1) for i = 1:1:dim]

if this_is_test # let test pass
    SpinValueType = Int8
end

function vector2number(v::Array{SpinValueType,1},weight::Array{Int,1})
    return sum(((v.+1)/2).*weight)+1
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

function measure_p!(model::Model,target::Target,p::Array{Int128,4})
    #TODO: This function has not been tested
    v = zeros(SpinValueType,dim)
    temp = zeros(Int,1+layers)
    for i in 1:1:NV
        v = number2vector!(v,i)
        temp[1] = i
        for j in 1:1:layers
            v = model.A[:,:,j]*v .+ model.B[:,j] .+ target.C[:,j]
            v = sign.(v)
            temp[1+j] = vector2number(v)
        end 
        for j in 1:1:layers+1
            for k in 1:1:layers+1
                p[j,k,temp[j],temp[k]]+=1
            end
        end
    end
end

function print_p(p::Array{Int128,4})
    #TODO: This function is only avaliable for n = 3 test!
        #TODO:This function need to be test
    #---- in order to save memory
    for i in 1:1:layers+1
        for j in 1:1:layers+1
            p[i,j,:,:] = p[i,j,:,:]/sum(p[i,j,:,:])
        end
    end

    s01 = sum(log.(p[1,2,:,:]).*p[1,2,:,:])
    s02 = sum(log.(p[1,3,:,:]).*p[1,3,:,:])
    s03 = sum(log.(p[1,4,:,:]).*p[1,4,:,:])
    s12 = sum(log.(p[2,3,:,:]).*p[2,3,:,:])
    s13 = sum(log.(p[2,4,:,:]).*p[2,4,:,:])
    s23 = sum(log.(p[3,4,:,:]).*p[3,4,:,:])
    s0 = sum(log.(sum(p[1,2,:,:],dims=2)).*sum(p[1,2,:,:],dims=1))
    s1 = sum(log.(sum(p[1,2,:,:],dims=1)).*sum(p[1,2,:,:],dims=2))
    s2 = sum(log.(sum(p[3,4,:,:],dims=2)).*sum(p[3,4,:,:],dims=1))
    s3 = sum(log.(sum(p[3,4,:,:],dims=1)).*sum(p[3,4,:,:],dims=2))

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