using Test
const SpinValueType = Int8
const layers = 3
const dim = 3
const NA = layers*dim^2 
const NB = layers*dim
const NV = 2^dim
const NC = 1/(NV*dim)

this_is_test = true
include("./MCModel.jl")
include("./MCStatistics.jl")
include("p_method.jl")

@testset "information set: p_method.jl" begin
v = Int8.([-1,-1,-1])
@test vector2number(v,weight) == 1
v = Int8.([ 1,-1,-1])
@test vector2number(v,weight) == 2
v = Int8.([-1, 1,-1])
@test vector2number(v,weight) == 3
v = Int8.([ 1, 1,-1])
@test vector2number(v,weight) == 4
v = Int8.([-1,-1, 1])
@test vector2number(v,weight) == 5
v = Int8.([ 1,-1, 1])
@test vector2number(v,weight) == 6
v = Int8.([-1, 1, 1])
@test vector2number(v,weight) == 7
v = Int8.([ 1, 1, 1])
@test vector2number(v,weight) == 8

for i in 1:1:2^dim
    number2vector!(v,i)
    @test vector2number(v,weight) == i
end

end


@testset "spectrum.jl" begin
    test_spectrum = Spectrum(10)
    set_ground!(test_spectrum,-1)
    for i in 1:1:10
        add_factor!(test_spectrum,1,1.0)
    end
    @test get_factor(test_spectrum,1)≈10
    @test get_factor(test_spectrum,0)≈0
    add_factor!(test_spectrum,0,10000000.0)
    @test normalize_factor(test_spectrum)==1
    @test process_spectrum(test_spectrum,-1)==1
    @test maximum(test_spectrum.factor)≈0.0
    @test minimum(exp.(test_spectrum.factor))≈0.0
end

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
