struct Spectrum
    length::Int
    count ::Array{Int128,1}
    level ::Array{Float64,1}
    function Spectrum(length::Int)
        count = zeros(Int128,length)
        level = zeros(Float64,length)
        new(length,count,level)
    end
end