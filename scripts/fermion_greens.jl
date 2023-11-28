using SmoQyDEAC
using FileIO
using Statistics

loadfile = joinpath("greens.jld2")
input_dictionary = load(loadfile)

Gτ_bin =  input_dictionary["Gτ"];
Gτ_std = mean(Gτ_bin,dims=1)[1,:];
Gτ_err = std(Gτ_bin,dims=1)[1,:];
Gω =  input_dictionary["Gω"];
Gω_std = mean(Gω,dims=1)[1,:];
Gω_err = std(Gω,dims=1)[1,:];

τs = collect(input_dictionary["τs"]); # must be evenly spaced.
β = τs[end];
ωₙ = collect(input_dictionary["ωns"]);

output_directory = "fermion_greens_output/";
try
    mkdir(output_directory);
catch
end

number_of_bins = 2;
runs_per_bin = 10 ;
output_file = joinpath(output_directory, "fermion_out.jld2");
checkpoint_directory = output_directory;
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));

base_seed = 1000000;
keep_bin_data = true;

output_dictionary_ωₙ = DEAC_Binned(Gτ_bin,β,ωₙ,ωs,"frequency_fermionic",number_of_bins,runs_per_bin,output_file,
                                  checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data)
output_dictionary_ωₙ_std = DEAC_Std(Gω_std,Gω_err,β,ωₙ,ωs,"frequency_fermionic",number_of_bins,runs_per_bin,output_file,
                                    checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data)

A = output_dictionary_τ["A"];
A_σ = output_dictionary_τ["σ"];
ωs_out = output_dictionary_τ["ωs"];
zeroth_calc = output_dictionary_τ["zeroth_moment"];
zeroth_σ = output_dictionary_τ["zeroth_moment_σ"];
avg_generations = output_dictionary_τ["avg_generations"];

bin_data = output_dictionary_τ["bin_data"];
bin_zeroth = output_dictionary_τ["bin_zeroth_moment"];

test_dictionary = FileIO.load(output_file)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
