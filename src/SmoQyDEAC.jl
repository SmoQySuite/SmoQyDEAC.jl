module SmoQyDEAC

using Random
using Statistics
using FileIO
using Printf
using LinearAlgebra
using JLD2
using LoopVectorization

include("DEAC/types.jl")
include("DEAC/deac.jl")
include("DEAC/jackknife.jl")
include("DEAC/kernels.jl")
include("DEAC/utility.jl")
include("DEAC/checkpoint.jl")
include("DEAC/boostrap.jl")

export DEAC_Std
export DEAC_Binned

end # module SmoQyDEAC
