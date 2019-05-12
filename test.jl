using Test
const SpinValueType = Int8
const layers = 3
const dim = 3
include(".\\MCModel.jl")
include(".\\MCStatistics.jl")

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
    @show exp.(test_spectrum.factor)
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
