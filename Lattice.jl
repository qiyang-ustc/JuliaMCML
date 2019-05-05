module Lattices
export Lattice, SquareLattice
using Random

abstract type Lattice{D,BC} end 
struct SquareLattice{BC} <: Lattice{2,BC} 
    Nbr::Matrix{Int}

    function SquareLattice{:periodic}(Lx::Int, Ly::Int) 
        #Neighbor table
        Nbr = zeros(Int, 4, Lx * Ly)
        sites = reshape(1:Lx * Ly, Lx, Ly)
        Nbr[1, :] .= reshape(circshift(sites, (-1, 0)), Lx * Ly)
        Nbr[2, :] .= reshape(circshift(sites, (0, -1)), Lx * Ly)
        Nbr[3, :] .= reshape(circshift(sites, (1, 0)), Lx * Ly)
        Nbr[4, :] .= reshape(circshift(sites, (0, 1)), Lx * Ly)
        new{:periodic}(Nbr)
    end 

    function SquareLattice{:open}(Lx::Int, Ly::Int) 
        #Neighbor table
        Nbr = zeros(Int, 4, Lx * Ly)
        for k in 0:Lx * Ly - 1
            x = mod(k, Lx) 
            y = div(k, Lx)
                
            #right
            if x + 1 <= Lx - 1
                Nbr[1, k + 1] = (x + 1) + y * Lx + 1 
            end

            #up 
            if y + 1 <= Ly - 1
                Nbr[2, k + 1] = x + (y + 1) * Lx + 1
            end

            #left
            if x - 1 >= 0
                Nbr[3, k + 1] = (x - 1) + y * Lx + 1
            end

            #down
            if y - 1 >= 0
                Nbr[4, k + 1] = x + (y - 1) * Lx + 1
            end
        end
        new{:open}(Nbr)
    end 
    

end
#display(SquareLattice{:open}(4, 4))
end 
