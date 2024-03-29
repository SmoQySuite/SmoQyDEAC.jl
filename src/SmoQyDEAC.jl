module SmoQyDEAC

using Random
using Statistics
using FileIO
using Printf
using LinearAlgebra
using JLD2
using LoopVectorization
using Chairmarks


include("DEAC/types.jl")
include("DEAC/deac.jl")
include("DEAC/jackknife.jl")
include("DEAC/kernels.jl")
include("DEAC/utility.jl")
include("DEAC/checkpoint.jl")
include("DEAC/boostrap.jl")
include("DEAC/updates.jl")

export DEAC_Std
export DEAC_Binned

function __init__()
    BLAS.set_num_threads(1)
end # __init__()

end # module SmoQyDEAC
