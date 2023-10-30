# List of allowed kernels
allowable_kernels = ["time_bosonic",
                     "time_fermionic",
                     "time_bosonic_symmetric",
                     "frequency_fermionic",
                     "frequency_bosonic",
                     "frequency_bosonic_symmetric"
                    ]



# h/t Przemyslaw Szufel
# https://stackoverflow.com/questions/62336686/struct-equality-with-arrays
abstract type Comparable end

import Base.==

function ==(a::T, b::T) where T <: Comparable
    f = fieldnames(T)
    getfield.(Ref(a),f) == getfield.(Ref(b),f)
end
                     

# DEACParameters
#
# A struct containing parameters passed in to the DEAC algorithm
#
mutable struct DEACParameters <: Comparable
    β::Float64
    input_grid::Vector{Float64}
    out_ωs::Vector{Float64}
    kernel_type::String
    output_file::String
    checkpoint_directory::String
    num_bins::Int64
    runs_per_bin::Int64
    population_size::Int64
    base_seed::Int64
    crossover_probability::Float64
    self_adapting_crossover_probability::Float64
    differential_weight::Float64
    self_adapting_differential_weight_probability::Float64
    self_adapting_differential_weight::Float64
    stop_minimum_fitness::Float64
    number_of_generations::Int64
end


# https://stackoverflow.com/questions/51956958/how-to-copy-a-struct-in-julia
Base.copy(x::DEACParameters) where DEACParameters = DEACParameters([getfield(x, k) for k ∈ fieldnames(DEACParameters)]...)