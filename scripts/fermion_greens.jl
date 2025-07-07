using SmoQyDEAC
using FileIO
using Statistics

loadfile = "greens.jld2";
input_dictionary = load(loadfile);

Gτ_bin =  input_dictionary["Gτ"];
Gτ_std = mean(Gτ_bin,dims=1)[1,:];
Gτ_err = std(Gτ_bin,dims=1)[1,:];
Gω_bin =  input_dictionary["Gω"];
Gω_std = mean(Gω_bin,dims=1)[1,:];
Gω_err = std(Gω_bin,dims=1)[1,:];

τs = collect(input_dictionary["τs"]); # must be evenly spaced.
β = τs[end];
ωₙ = collect(input_dictionary["ωns"]);

output_directory = "fermion_greens_output/";
mkpath(output_directory);

number_of_bins = 2;
runs_per_bin = 5 ;
output_file = joinpath(output_directory, "fermion_out.jld2");
checkpoint_file = joinpath(output_directory,"DEAC_checkpoint.jld2");
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));

base_seed = 1000000;
keep_bin_data = true;

output_dictionary_τ = DEAC_Binned(
    Gτ_bin,
    β,
    τs,
    ωs,
    "time_fermionic",
    number_of_bins,
    runs_per_bin,
    output_file,
    checkpoint_file;
    base_seed = base_seed,
    keep_bin_data = keep_bin_data,
    verbose = true
)
output_dictionary_τ_std = DEAC_Std(
    Gτ_std,
    Gτ_err,
    β,
    τs,
    ωs,
    "time_fermionic",
    number_of_bins,
    runs_per_bin,
    output_file,
    checkpoint_file;
    base_seed = base_seed,
    find_fitness_floor=false,
    keep_bin_data = keep_bin_data,
    number_of_generations = 10000,
    verbose = true
)
output_dictionary_ωₙ = DEAC_Binned(
    Gω_bin,
    β,
    ωₙ,
    ωs,
    "frequency_fermionic",
    number_of_bins,
    runs_per_bin,
    output_file,
    checkpoint_file;
    base_seed = base_seed,
    keep_bin_data = keep_bin_data,
    fitness=1.0,
    find_fitness_floor=false,
    number_of_generations = 20000,
    verbose = true
)
output_dictionary_ωₙ_std = DEAC_Std(
    Gω_std,
    Gω_err,
    β,
    ωₙ,
    ωs,
    "frequency_fermionic",
    number_of_bins,
    runs_per_bin,
    output_file,
    checkpoint_file;
    base_seed = base_seed,
    keep_bin_data = keep_bin_data,
    fitness = 1.0,
    find_fitness_floor = false,
    number_of_generations = 20000,
    verbose = true
)

A = output_dictionary_τ["A"];
fitnesses = output_dictionary_τ["fitness"];
A_σ = output_dictionary_τ["σ"];
ωs_out = output_dictionary_τ["ωs"];
zeroth_calc = output_dictionary_τ["zeroth_moment"];
zeroth_σ = output_dictionary_τ["zeroth_moment_σ"];
avg_generations = output_dictionary_τ["avg_generations"];

bin_data = output_dictionary_τ["bin_data"];
bin_zeroth = output_dictionary_τ["bin_zeroth_moment"];
full_eigenvalues = output_dictionary_τ["full_eigenvalues"];

test_dictionary = FileIO.load(output_file)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
