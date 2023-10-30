
# χ² fit 
function Χ²(observed::AbstractVector{T},calculated::AbstractMatrix{T},W::AbstractVector{T}) where {T<:Real}
    Χ = zeros(T,(size(calculated,2),))
    for pop in 1:size(calculated,2)
        @inbounds Χ[pop] = sum( ((observed .- calculated[:,pop]).^2) .* W )
    end
    return Χ
end # χ²

function Χ²(observed::AbstractVector{U},calculated::AbstractMatrix{U},W::AbstractVector{T}) where {T<:Real, U<:Complex}
    Χ = zeros(T,(size(calculated,2),))
    for pop in 1:size(calculated,2)
        @inbounds Χ[pop] = sum( norm.((observed .- calculated[:,pop]).^2) .* W )
    end
    return Χ
end # χ²



# return mutant indices
# links our genomes together for mutation
function get_mutant_indices(rng,pop_size)
    indices = zeros(Int64,(3,pop_size))
    for pop in 1:pop_size
        indices[:,pop] .= pop
        while (indices[1,pop] == pop)
            indices[1,pop] = 1 + mod(Random.rand(rng,Int64),pop_size)
        end
        while (indices[2,pop] == pop) || (indices[2,pop] == indices[1,pop])
            indices[2,pop] = 1 + mod(Random.rand(rng,Int64),pop_size)
        end
        while (indices[3,pop] == pop) || (indices[3,pop] == indices[1,pop]) || (indices[3,pop] == indices[2,pop])
            indices[3,pop] = 1 + mod(Random.rand(rng,Int64),pop_size)
        end
    end
    return indices
end # get_mutant_indices

# Loop vectorized Matrix Multiply
function gemmSIMD!(C::Matrix{T}, A::Matrix{T}, B::Matrix{T}) where {T}
    @turbo for m ∈ axes(A,1), n ∈ axes(B,2)
        Cmn = zero(T)
        for k ∈ axes(A,2)
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end

# Determine if to use gemmavx! or LinearAlgebra mul!
function useSIMD(K,M,N)

    # Matrices for GEMM
    C = Matrix{Float64}(undef, M, N)
    A = Random.rand(Float64,(M, K))
    B = Random.rand(Float64,(K, N))

    # Benchmark
    t_linAlg = @belapsed mul!($C,$A,$B)
    t_avx = @belapsed gemmSIMD!($C,$A,$B)
    use_SIMD = t_avx < t_linAlg
    if use_SIMD
        println("SIMD GEMM faster than BLAS, using internal gemmSIMD!()")
    else
        println("SIMD GEMM slower than BLAS, using BLAS mul!()")
    end
    
    return use_SIMD
end

# GEMM wrapper to select appropriate method
function GEMM!(C,A,B,use_SIMD)
    if use_SIMD && eltype(A) <: Real
        gemmSIMD!(C,A,B)
    else
        mul!(C,A,B)
    end
end

# Calculate binned data and save to a checkpoint
function bin_results!(bin_data,calculated_zeroth_moment,run_data,curbin,Δω,Greens_tuple,true_fitness,seed_vec,generations,params)
    bin_data[:,curbin] = run_data[:,curbin] / params.runs_per_bin
    calculated_zeroth_moment[1,curbin] = sum(bin_data[:,curbin]) .* Δω
        
    # Bosonic time kernels steal a factor of ω from the spectral function.
    # Multiply it back in if needed
    if  occursin("bosonic",params.kernel_type)
        bin_data[:,curbin] = bin_data[:,curbin] .* params.out_ωs
    end
    
    if curbin != params.num_bins
        save_checkpoint(bin_data,generations,curbin,params,Greens_tuple,calculated_zeroth_moment,true_fitness,seed_vec)
    end
    
    println("Finished bin ",curbin," of ",params.num_bins)
    
end

# Calculate matrices used to go from ω to τ space and χ² fit
function calculate_fit_matrices(Greens_tuple,K,W_ratio_max,use_SIMD)
    if Greens_tuple[2] == nothing
        # Covariance Methods
        
        # SVD on correlation bins
        corr_avg = Statistics.mean(Greens_tuple[1],dims=1)
        svd_corr = svd(Greens_tuple[1] .- corr_avg)
        sigma_corr = svd_corr.S

        # Unitary transformation matrix
        U_c = svd_corr.Vt
        
        # Inverse fit W array for χ^2
        # (2.0 * U_c1) factor generally gives ideal fit ~0.1-1.0
        U_c1 = size(U_c,1)
        W = (2.0 * U_c1) ./ (sigma_corr .* sigma_corr) 
        
        # Deal with nearly singular matrix
        W_cap = W_ratio_max * minimum(norm(W))
        clamp!(W,0.0,W_cap)
        
        Kp = similar(K)
        
        # rotate K and corr_avg
        GEMM!(Kp,U_c,K,use_SIMD)
        
        corr_avg_p = zeros(Float64,U_c1)
        for i in 1:U_c1
            corr_avg_p[i] = dot(view(U_c,i,:),corr_avg)
        end
    else
        # Diagonal error method
        W = norm.(1.0 ./ (Greens_tuple[2] .* Greens_tuple[2]))
        W_cap = W_ratio_max * minimum(norm(W))
        clamp!(W,0.0,W_cap)
        Kp = K
        corr_avg_p = Greens_tuple[1]
    end
    return W, Kp, corr_avg_p

end


# function clamp!(A::AbstractVector{T<:Complex},low::Real,hi::Real)
#     for k in axes(A,1)
#         if 

