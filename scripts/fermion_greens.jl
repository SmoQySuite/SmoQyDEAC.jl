using SmoQyDEAC
using FileIO

loadfile = joinpath(pkgdir(SmoQyDEAC), "docs/src/examples/greens.jld2")
input_dictionary = load(loadfile)

G_std = input_dictionary["G_std"];
G_error = input_dictionary["G_err"];
G_bin =  input_dictionary["G_bin"];
τs = input_dictionary["τs"]; # must be evenly spaced.
β = input_dictionary["β"];

output_directory = "fermion_greens_output/";
try
    mkdir(output_directory);
catch
end

number_of_bins = 2;
runs_per_bin = 1 ;
output_file = joinpath(output_directory, "fermion_out.jld2");
checkpoint_directory = output_directory;
nω = 401;
ωmin = -10.;
ωmax = 10.;
ωs = collect(LinRange(ωmin,ωmax,nω));

base_seed = 1000000;
keep_bin_data = true;

output_dictionary = DEAC_Binned(G_bin,β,τs,ωs,"time_fermionic",number_of_bins,runs_per_bin,output_file,
                         checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data)
output_dictionary_std = DEAC_Std(G_std,G_error,β,τs,ωs,"time_fermionic",number_of_bins,runs_per_bin,output_file,
                         checkpoint_directory,base_seed=base_seed,keep_bin_data=keep_bin_data)

A = output_dictionary["A"];
A_σ = output_dictionary["σ"];
ωs_out = output_dictionary["ωs"];
zeroth_calc = output_dictionary["zeroth_moment"];
zeroth_σ = output_dictionary["zeroth_moment_σ"];
avg_generations = output_dictionary["avg_generations"];

bin_data = output_dictionary["bin_data"];
bin_zeroth = output_dictionary["bin_zeroth_moment"];

test_dictionary = FileIO.load(output_file)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
