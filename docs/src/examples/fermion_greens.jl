# # Example 1: Fermion Greens function
#
#  Usage:
#
#   `$ julia --threads=auto fermion_greens.jl`
#
#   SmoQyDEAC uses multithreading for parallelizing runs. Multithreading it recommended.
#   --threads=auto will run a thread for each core available, or 2x for hyperthreading cores
#
# In this example we will take the up-spin electron Green's function output from a Determinant Quantum Monte Carlo run
# generated via [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl) using the Holstein Model. 
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

# We now load the data provided in our source file.
loadfile = joinpath(pkgdir(SmoQyDEAC), "docs/src/examples/greens.jld2")
input_dictionary = load(loadfile)

G_std = input_dictionary["G_std"];
G_error = input_dictionary["G_err"];
G_bin =  input_dictionary["G_bin"];
τs = input_dictionary["τs"]; # must be evenly spaced.
β = input_dictionary["β"];

# Make an output folder for checkpoint file and output file
output_directory = "fermion_greens_output/";
try
    mkdir(output_directory);
catch
end

# Define necessary parameters for the DEAC run
# Typically you will want at least 1,000 for number_of_bins * runs_per_bin
# For speed's sake we only do 2*1 in this example. 
number_of_bins = 2;
runs_per_bin = 1 ;
output_file = joinpath(output_directory, "fermion_out.jld2");
checkpoint_directory = output_directory;
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));


# Set optional parameters
base_seed = 1000000;
#md ## Note, the seed will incement for each run. 
#md ## Starting a new run at 1000002 will have output unique from this run
keep_bin_data = true;
#md ## If true, each bin will have it's data written to the output dictionary
#md ## Set to false to save disk space

# Run DEAC Algorithm for binned and unbinned
output_dictionary = DEAC_Binned(G_bin,β,τs,ωs,"time_fermionic",number_of_bins,runs_per_bin,output_file,
                         checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data)
output_dictionary_std = DEAC_Std(G_std,G_error,β,τs,ωs,"time_fermionic",number_of_bins,runs_per_bin,output_file,
                         checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data)

# Accessing output
#md ## Spectral function, 1D array size (nω)
A = output_dictionary["A"];
#md ## Spectral function error, 1D array size (nω)
A_σ = output_dictionary["σ"];
#md ## ω grid, 1D array size (nω)
ωs_out = output_dictionary["ωs"];
#md ## zeroth moment: For fermions it is G(0) + G(β) which should = 1.0. Float64s
zeroth_calc = output_dictionary["zeroth_moment"];
zeroth_σ = output_dictionary["zeroth_moment_σ"];
#md ## Number of average generations to converge, Float64
avg_generations = output_dictionary["avg_generations"];

# Binned information - not available if `keep_bin_data=false`
#md ## Bin data, 2D array size (nω,nbins)
bin_data = output_dictionary["bin_data"];
#md ## Standard error for each bin, 2D array size (nω,nbins)
bin_σ = output_dictionary["bin_σ"];
#md ## zeroth moment, 1D array (nbins)
bin_zeroth = output_dictionary["bin_zeroth_moment"];


# The dictionary will automatically be saved
#md ## Example of loading the data from the jld2
test_dictionary = FileIO.load(output_file)
# This is identical to output_dictionary



