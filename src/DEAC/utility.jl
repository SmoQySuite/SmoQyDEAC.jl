
# χ² fit 
function Χ²(observed::AbstractVector,calculated::AbstractMatrix,W::AbstractVector)
    Χ = zeros(Float64,(size(calculated,2),))
    for pop in 1:size(calculated,2)
        Χ[pop] = sum( ((observed .- calculated[:,pop]).^2) .* W )
    end
    return Χ
end # χ²

# return mutant indices
# links our genomes together for mutation
function get_mutant_indices(rng,pop_size)
    indices = zeros(Int64,(3,pop_size))
    for pop in 1:pop_size
        indices[:,pop] .= pop
        while (indices[1,pop] == pop)
            indices[1,pop] = 1 + mod(Random.rand(rng,Int64),pop_size)
        end
        while (indices[2,pop] == pop) || (indices[2,pop] == indices[1,pop])
            indices[2,pop] = 1 + mod(Random.rand(rng,Int64),pop_size)
        end
        while (indices[3,pop] == pop) || (indices[3,pop] == indices[1,pop]) || (indices[3,pop] == indices[2,pop])
            indices[3,pop] = 1 + mod(Random.rand(rng,Int64),pop_size)
        end
    end
    return indices
end # get_mutant_indices
