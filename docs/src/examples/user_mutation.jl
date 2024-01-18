# # Example 3: User Mutation
#
#  Usage:
#
#   `$ julia --threads=auto user_mutation.jl`
#
#   SmoQyDEAC uses multithreading for parallelizing runs. Multithreading it recommended.
#   --threads=auto will run a thread for each core available, or 2x for hyperthreading cores
#
# In this example we will show how to add an additional user-defined mutation to the DEAC base algorithm
# 

#md ## First we import all required packages
using SmoQyDEAC
using FileIO
using Statistics
using LoopVectorization
using Random

# We now load the data provided in our source file. This is the same as example 1, but we are just loading the binned τ data.
loadfile = "greens.jld2";
input_dictionary = load(loadfile);
Gτ_bin =  input_dictionary["Gτ"];
τs = collect(input_dictionary["τs"]);
β = τs[end];

# Make an output folder for checkpoint file and output file
output_directory = "user_mutation_output/";
mkpath(output_directory);

# Define necessary parameters for the DEAC run.
# Typically you will want at least 1,000 for number_of_bins * runs_per_bin. 
# For speed's sake we only do 2*5 in this example. 
number_of_bins = 2;
runs_per_bin = 5 ;
output_file = joinpath(output_directory, "user_mutation_out.jld2");
checkpoint_directory = output_directory;
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));


# Set optional parameters
base_seed = 1000000;
keep_bin_data = false;

# Define an additional mutation function which takes two populations of [1:nω,1:n_population] size and a random number generator.
# Per Julia convention, the first argument is the one which is updated. For those new to Julia, the "!" at the end of the function name means the function will mutate at least one of the arguments.
# This mutation just puts a slight blur on odd indexed ω values. This is not something you'd want to do, but, this is a simple example
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

# Vectorized version of the same function using the LoopVectorization.jl package.
function user_mut_vec!(pop_new,pop_old,rng)
    rng_pop = rand(rng,Float64,(size(pop_new,1),size(pop_new,2)))
    @turbo for pop ∈ axes(pop_new,2), ω ∈ axes(pop_new,1)
        pop_new[ω,pop] = pop_old[ω,pop] + (0.05 * (ω%2) * (0.5 - rng_pop[ω,pop]))
    end
    return nothing
end


# Now we run DEAC with the additional mutation
# Note, number_of_generations is set low for the example's run time
output_dictionary_τ = DEAC_Binned(
    Gτ_bin,
    β,
    τs,
    ωs,
    "time_fermionic",
    number_of_bins,
    runs_per_bin,
    output_file,
    checkpoint_directory;
    base_seed = base_seed,
    keep_bin_data = keep_bin_data,
    user_mutation! = user_mut_vec!,
    stop_minimum_fitness = 1.0,
    find_ideal_fitness = false,
    number_of_generations = 10000,
    verbose = true
)
