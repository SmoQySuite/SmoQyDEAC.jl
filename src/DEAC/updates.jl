function initial_fit!(pop,model,avg_p,Kp,W,params,use_SIMD)
    for p in axes(pop,2)
        pop[:,p] = pop[:,p] ./ sum(pop[:,p])
    end
    
    GEMM!(model,Kp,pop,use_SIMD)
    
    fitness = Χ²(avg_p,model,W) ./ size(params.input_grid,1)

    return fitness

end



function update_weights!(cp_new,dw_new,cp_old,dw_old,rng,params)
    for p in 1:params.population_size
        cp_new[p] = (Random.rand(rng,Float64)<params.self_adapting_crossover_probability) ? rand(rng,Float64) : cp_old[p]
        dw_new[p] = (Random.rand(rng,Float64)<params.self_adapting_differential_weight_probability) ? 2.0*rand(rng,Float64) : dw_old[p]
    end

end

function rand_mutate_array!(mutate_indices,crossover_probability_new,rng,params)
    mutate_indices_rnd = Random.rand(rng,Float64, (size(params.out_ωs,1),params.population_size)) 
    for pop in 1:params.population_size
        @. mutate_indices[:,pop] = Float64(mutate_indices_rnd[:,pop] < crossover_probability_new[pop])
    end
end


# Population update routine
#
# Below is the unoptimized verson which is easier to read
#   Note, mutate_indices[pop,ω] here is a bool, in the actual function it is 1.0 or 0.0
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
function propose_populations!(population_new,population_old,mutate_indices,differential_weights_new,mutant_indices,params)
    @turbo for pop in 1:params.population_size, ω in 1:size(params.out_ωs,1)
        population_new[ω,pop] = population_old[ω,pop]*(1.0-mutate_indices[ω,pop]) + # If false keep old gene
                                #if true do update
                                mutate_indices[ω,pop] * abs(population_old[ω,mutant_indices[1,pop]] + differential_weights_new[pop]*
                                (population_old[ω,mutant_indices[2,pop]]-population_old[ω,mutant_indices[3,pop]]))
                                
    end # pop
end


function update_populations!(fitness_old,crossover_probability_old,differential_weights_old,population_old,fitness_new,crossover_probability_new,differential_weights_new,population_new)
    for pop in axes(fitness_old,1)
        if fitness_new[pop] <= fitness_old[pop]
            fitness_old[pop] = fitness_new[pop]
            crossover_probability_old[pop] = crossover_probability_new[pop]
            differential_weights_old[pop] = differential_weights_new[pop]
            population_old[:,pop] = population_new[:,pop]
        end
    end
end


