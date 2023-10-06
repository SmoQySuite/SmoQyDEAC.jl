
# χ² fit 
function Χ²(observed::AbstractVector{T},calculated::AbstractMatrix{T},W::AbstractVector{T}) where {T<:Real}
    Χ = zeros(T,(size(calculated,2),))
    for pop in 1:size(calculated,2)
        @inbounds Χ[pop] = sum( ((observed .- calculated[:,pop]).^2) .* W )
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

# Population update routine
#
# Below is the unoptimized verson which is easier to read
#
# for pop in 1:params.population_size
#     for ω in 1:size(params.out_ωs,1)
#         if mutate_indices[pop,ω]
#             population_new[ω,pop] = abs(population_old[ω,mutant_indices[1,pop]] + differential_weights_new[pop]*
#                                         (population_old[ω,mutant_indices[2,pop]]-population_old[ω,mutant_indices[3,pop]]))
#         else
#             population_new[ω,pop] = population_old[ω,pop]
#         end
#     end
# end
function update_populations!(params,population_new,population_old,mutate_indices,differential_weights_new,mutant_indices)
    @turbo for pop in 1:params.population_size, ω in 1:size(params.out_ωs,1)
        do_upd = convert(eltype(population_new),mutate_indices[ω,pop])
        population_new[ω,pop] = population_old[ω,pop]*(1.0-do_upd) + # If false keep old gene
                                #if true do update
                                do_upd * abs(population_old[ω,mutant_indices[1,pop]] + differential_weights_new[pop]*
                                (population_old[ω,mutant_indices[2,pop]]-population_old[ω,mutant_indices[3,pop]]))
                                
    end # pop
end

# Loop vectorized Matrix Multiply
function gemmavx!(C::Matrix{T}, A::Matrix{T}, B::Matrix{T}) where {T}
    @turbo for m ∈ axes(A,1), n ∈ axes(B,2)
        Cmn = zero(T)
        for k ∈ axes(A,2)
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end

