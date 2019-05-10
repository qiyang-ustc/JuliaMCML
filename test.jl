using Test
const SpinValueType = Int8
const layers = 3
const dim = 3
include(".\\MCModel.jl")

@testset "model" begin
target = Target(dim,layers)
model = Model(dim,layers)
temp = model.A[1]
picker!(model,1)
@test temp == -model.A[1]

temp = model.B[1]
picker!(model,dim^2*layers+1)
@test temp == -model.B[1]
for i in 1:dim
    for j in 1:dim
        for k in 1:layers
            model.A[i,j,k] = target.A[i,j,k]
        end 
    end
end
for i in 1:dim
    for k in 1:layers
        model.B[i,k] = target.B[i,k]
    end 
end
vectors = create_bit_verctor(dim)
@test cal_energy(model,target,vectors)==0
end

@testset "spectrum.jl" begin

end