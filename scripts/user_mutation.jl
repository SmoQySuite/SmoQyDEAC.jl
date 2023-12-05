using SmoQyDEAC
using FileIO
using Statistics
using LoopVectorization
using Random

loadfile = "greens.jld2";
input_dictionary = load(loadfile);
Gτ_bin =  input_dictionary["Gτ"];
τs = collect(input_dictionary["τs"]);
β = τs[end];

output_directory = "user_mutation_output/";
mkpath(output_directory);

number_of_bins = 2;
runs_per_bin = 5 ;
output_file = joinpath(output_directory, "user_mutation_out.jld2");
checkpoint_directory = output_directory;
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));

base_seed = 1000000;
keep_bin_data = false;

function user_mut!(pop_new,pop_old,rng)
    npop = size(pop_old,2)
    nω = size(pop_old,1)
    for ω ∈ 1:nω
        if (ω%2) == 1
            pop_new[ω,:] = pop_old[ω,:] .+ 0.05 .* (0.5 .- rand(rng,Float64,npop))
        else
            pop_new[ω,:] = pop_old[ω,:]
        end
    end
    return nothing
end

function user_mut_vec!(pop_new,pop_old,rng)
    rng_pop = rand(rng,Float64,(size(pop_new,1),size(pop_new,2)))
    @turbo for pop ∈ axes(pop_new,2), ω ∈ axes(pop_new,1)
        pop_new[ω,pop] = pop_old[ω,pop] + (0.05 * (ω%2) * (0.5 - rng_pop[ω,pop]))
    end
    return nothing
end

output_dictionary_τ = DEAC_Binned(Gτ_bin,β,τs,ωs,"time_fermionic",number_of_bins,runs_per_bin,output_file,
                                  checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data,
                                  user_mutation! =user_mut_vec!,stop_minimum_fitness=20.0,
                                  find_ideal_fitness=false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
