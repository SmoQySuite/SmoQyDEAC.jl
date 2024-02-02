# # Example 2: Load from SmoQyDQMC and run SmoQyDEAC
#
#  Usage:
#
#   `$ julia --threads=auto SmoQyDQMC.jl`
#
#   SmoQyDEAC uses multithreading for parallelizing runs. Multithreading it recommended.
#   `--threads=auto` will run a thread for each core available, or 2x for hyperthreading cores
#
# In this example we explicitly load output from [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl).
# After running a simulation with [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl) you will need to use the 
# [`correlation_bins_to_csv`](https://smoqysuite.github.io/SmoQyDQMC.jl/dev/api/#SmoQyDQMC.correlation_bins_to_csv) method for the
# desired correlation function prior to using the loading script.
#
# This example utilizes the script [`SmoQyDQMCloader.jl`](https://github.com/SmoQySuite/SmoQyDEAC.jl/blob/main/scripts/SmoQyDQMCloader.jl) to parse
# [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl) csv files. This script may be moved to another repository at a future date.


#md ## First we import our SmoQyDQMC csv parser and needed packages for the example
include("SmoQyDQMCloader.jl")
using FileIO
using SmoQyDEAC


# Create directory for outputs 
output_directory = "SmoQyDQMC_DEAC_Out/";
mkpath(output_directory);


# Load data from fermion greens correlation functions
# This puts `real` in the format 
# `real["ORBITAL\_ID\_1","ORBITAL\_ID\_2","TAU","K\_1","K\_2","K\_3","BIN","PID"]``
#
# This example was from a 1D Holstein run with β = 20.0,
#   ORBITAL\_ID\_1 ∈ {1},
#   ORBITAL\_ID\_2 ∈ {1},
#   TAU          ∈ {0.0,Δτ,...,β-Δτ,β} for Nτ = 201,
#   Kx           ∈ {1,...,32},
#   Ky           ∈ {1},
#   Kz           ∈ {1},
#   Bin          ∈ {1,...,Nbin} for Nbin = 100,
#   PID          ∈ {1},
#
# Some dimensions are 1 deep. They are kept to ensure generality.
#
# See the [`SmoQyDQMCloader.jl`](https://github.com/SmoQySuite/SmoQyDEAC.jl/blob/main/scripts/SmoQyDQMCloader.jl) file for more information
input_directory = "SmoQyDQMC_sim-1/"

dims,real,image,sgnr,sgni,β = load_from_SmoQyDQMC(simulationfolder=input_directory,
                                                correlation="greens",
                                                space="momentum",
                                                type="time_displaced",bin=true);
#md ## dims is a dictionary which tells you what each dimension corresponds to.
println(dims) 



# Eliminate unnecessary dimensions
Gτ = real[1,1,:,:,1,1,:,1];
println(size(Gτ))

# Set parameters for DEAC
τs = collect(LinRange(0.0,β,size(Gτ,1)))
Nkx = size(Gτ,2)
number_of_bins = 2;
runs_per_bin = 10 ;
checkpoint_directory = output_directory;
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));


# Run DEAC over all k points in the x direction.
# For speed in this example I run 1:1 instead of 1:Nkx
for kx in 1:1 # 1:Nkx 
    output_file = joinpath(output_directory, string(kx) * ".jld2");
#md     ## put in [bins,τ] shape
    Gτ_temp = Matrix{Float64}(Gτ[:,kx,:]');
    deac_dict = DEAC_Binned(
        Gτ_temp,
        β,
        τs,
        ωs,
        "time_fermionic",
        number_of_bins,
        runs_per_bin,
        output_file,
        checkpoint_directory;
        stop_minimum_fitness = 1.0,
        find_ideal_fitness = false,
        number_of_generations = 20000,
        verbose = true
    )
end

# Note, these did not converge to a fitness of 1.0 within 20,000 generations. The number of generations is limited for speed when running this example.