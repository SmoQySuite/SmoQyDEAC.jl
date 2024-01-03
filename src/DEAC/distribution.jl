function initialize_population(rng,params; user_init_dist=nothing)
    if user_init_dist != nothing

        population_old = use_user_init_dist(rng,params,user_init_dist)
    else
        population_old = default_dist(rng,params)
    end
    population_new = zeros(eltype(population_old),size(population_old))
    return population_old, population_new
end


function default_dist(rng,params)
    population_old  = reshape(Random.rand(rng,size(params.out_ωs,1)*params.population_size),(size(params.out_ωs,1),params.population_size))
    if occursin("bosonic",params.kernel_type)
        ωs_mod = copy(params.out_ωs)
        ωs_mod[ωs_mod .≈ 0.0] .= typemax(eltype(params.out_ωs))
        pop_mod = 1.0 ./ ωs_mod
        for pop ∈ 1:params.population_size
            population_old[:,pop] .* pop_mod
        end
    end

    return population_old
end










function use_user_init_dist(rng,params,user_init_dist)

end
