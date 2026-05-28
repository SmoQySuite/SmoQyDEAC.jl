# # Example 1: Fermion Greens function
#
#  Usage:
#
#   `$ julia --threads=auto fermion_greens.jl`
#
#   SmoQyDEAC uses multithreading for parallelizing runs. Multithreading it recommended.
#   `--threads=auto` will run a thread for each core available, or 2x for hyperthreading cores
#
# In this example we will take the up-spin electron Green's function output from a Determinant Quantum Monte Carlo run
# generated via [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl) using the Holstein Model. For information on how to load directly from SmoQyDQMC see example 2.
# 
# We note the convention that the correlation function reported as the Green's Function has ħ=1 and there is no leading negative sign
# ```math
# G(k,τ)=⟨T_{τ}c_k(τ)c_k^†(0)⟩
# ```
# Our relation uses the time fermionic kernel such that
# ```math
# G(k,τ)=\int_{-∞}^∞ dω K(ω,τ)A(ω)=\int_{-∞}^∞ dω \frac{e^{-τω}}{1+e^{-ωβ}}A(ω)
# ```
# Since $A(ω)=-ℑG(ω)/π$ both negative signs that would normally be in the expression and factors of π cancel.
#
# On to the example:

#md ## First we import all required packages
using SmoQyDEAC
using FileIO
using Statistics

# We now load the data provided in our source file.
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

# Make an output folder for checkpoint file and output file
output_directory = "fermion_greens_output/";
mkpath(output_directory);

# Define necessary parameters for the DEAC run
# Typically you will want at least 1,000 for `number_of_bins` * `runs_per_bin`
# For speed's sake we only do 2*5 in this example. 
number_of_bins = 2;
runs_per_bin = 5 ;
output_file = joinpath(output_directory, "fermion_out.jld2");
checkpoint_file = joinpath(output_directory,"DEAC_checkpoint.jld2");
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));


# Set optional parameters
base_seed = 1000000;
#md ## Note, the seed will incement for each run. 
#md ## Starting a new run at 1000020 will have output unique from this run
keep_bin_data = true;
#md ## If true, each bin will have it's data written to the output dictionary
#md ## Set to false to save disk space

# Run DEAC Algorithm for binned and unbinned data for τ and ωₙ spaces
#md ## Note, number_of_generations should typically be > 100,000 or even 1,000,000.
#md ## It is set lower to reduce run time in this example.
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

# Accessing output
#md ## Spectral function, 2D array size (nω,nFitness)
#md ## nFitness is the number of fitnesses you set. In this case, 1
A = output_dictionary_τ["A"];
fitnesses = output_dictionary_τ["fitness"];
#md ## Spectral function error, 1D array size (nω)
A_σ = output_dictionary_τ["σ"];
#md ## ω grid, 1D array size (nω)
ωs_out = output_dictionary_τ["ωs"];
#md ## zeroth moment: For fermions it is G(0) + G(β) which should = 1.0. Float64s
zeroth_calc = output_dictionary_τ["zeroth_moment"];
zeroth_σ = output_dictionary_τ["zeroth_moment_σ"];
#md ## Number of average generations to converge, Float64
avg_generations = output_dictionary_τ["avg_generations"];

# Binned information - not available if `keep_bin_data=false`
#md ## Bin data, 2D array size (nω,nbins)
bin_data = output_dictionary_τ["bin_data"];
#md ## zeroth moment, 1D array (nbins)
bin_zeroth = output_dictionary_τ["bin_zeroth_moment"];
#md ## full eigenvalues - only applies to DEAC_Binned
#md ## These are the eigenvalues of the covariance matrix.
#md ## Eigenvectors whose eigenvalues are a factor eigenvalue_ratio_min less than the highest eigenvalue
#md ## are dropped since not all correlation function inputs are linearly independent
full_eigenvalues = output_dictionary_τ["full_eigenvalues"];


# The dictionary will automatically be saved
#md ## Example of loading the data from the jld2
test_dictionary = FileIO.load(output_file)
#md ## This is identical to output_dictionary_ωₙ_std



