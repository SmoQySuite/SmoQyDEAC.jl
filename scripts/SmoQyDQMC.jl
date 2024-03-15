include("SmoQyDQMCloader.jl")
using FileIO
using SmoQyDEAC

output_directory = "SmoQyDQMC_DEAC_Out/";
mkpath(output_directory);

input_directory = "SmoQyDQMC_sim-1/"

dims,real,image,sgnr,sgni,β = load_from_SmoQyDQMC(simulationfolder=input_directory,
                                                correlation="greens",
                                                space="momentum",
                                                type="time_displaced",bin=true);
println(dims)

Gτ = real[1,1,:,:,1,1,:,1];
println(size(Gτ))

τs = collect(LinRange(0.0,β,size(Gτ,1)))
Nkx = size(Gτ,2)
number_of_bins = 2;
runs_per_bin = 10 ;
checkpoint_file = joinpath(output_directory,"DEAC_checkpoint.jld2");
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));

for kx in 1:1 # 1:Nkx
    output_file = joinpath(output_directory, string(kx) * ".jld2");
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
        checkpoint_file;
        fitness = [1.5,1.0],
        find_fitness_floor = false,
        number_of_generations = 20000,
        verbose = true
    )
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

