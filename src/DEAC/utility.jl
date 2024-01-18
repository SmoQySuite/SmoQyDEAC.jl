
# χ² fit for real correlation functions
function Χ²(observed::AbstractVector{T},calculated::AbstractMatrix{T},W::AbstractVector{T}) where {T<:Real}
    Χ = zeros(T,(size(calculated,2),))
    for pop in 1:size(calculated,2)
        @inbounds Χ[pop] = sum( ((observed .- calculated[:,pop]).^2) .* W )
    end
    return Χ
end # χ²

# χ² fit for complex correlation functions
function Χ²(observed::AbstractVector{U},calculated::AbstractMatrix{U},W::AbstractVector{T}) where {T<:Real, U<:Complex}
    Χ = zeros(T,(size(calculated,2),))
    for pop in 1:size(calculated,2)
        @inbounds Χ[pop] = sum( (normsq.(observed .- calculated[:,pop])) .* W )
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
function gemmSIMD!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
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
    BenchmarkTools.DEFAULT_PARAMETERS.seconds=1.0
    t_linAlg = @belapsed mul!($C,$A,$B)
    t_avx = @belapsed gemmSIMD!($C,$A,$B)
    use_SIMD = t_avx < t_linAlg
    if use_SIMD
        println("SIMD GEMM faster than BLAS, using SmoQyDEAC's gemmSIMD!()")
    else
        println("SIMD GEMM slower than BLAS, using BLAS mul!()")
    end
    
    return use_SIMD
end

# GEMM wrapper to select appropriate method
function GEMM!(C,A,B,use_SIMD)
    if use_SIMD && eltype(C) <: Real
        gemmSIMD!(C,A,B)
    else
        mul!(C,A,B)
    end
end

# Calculate binned data and save to a checkpoint
function bin_results!(bin_data,calculated_zeroth_moment,run_data,weight_data,curbin,Δω,Greens_tuple,true_fitness,seed_vec,generations,params)
    bin_data[:,curbin] = run_data[:,curbin] ./  sum(weight_data[curbin])
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
function calculate_fit_matrices(Greens_tuple,K,use_SIMD,bootstrap,params)
    
    if Greens_tuple[2] == nothing
        # Covariance Methods
              
        corr_avg = Statistics.mean(Greens_tuple[1],dims=1)
        mask = get_covariance_mask(params)

        
        # Find eigenbasis for covariance matrix 
        cov_matrix = Statistics.cov(Greens_tuple[1],dims=1,corrected=true)
        F = eigen(cov_matrix)
        max_eig = maximum(F.values)

        # catch near-zero eigenvalues not found using get_covariance_mask
        mask = mask .&& (F.values .> max_eig * 1e-8)
        
        U = F.vectors[:,mask]
        corr_avg_p = Array{eltype(corr_avg)}(undef,size(corr_avg[:,mask]))
        
        # Rotate correlation functions
        GEMM!(corr_avg_p,corr_avg,U,use_SIMD)
        corr_avg_p = corr_avg_p[1,:]
        
        Kp = similar(K[mask,:])
        GEMM!(Kp,transpose(U),K,use_SIMD)
        
        Nsteps = size(corr_avg_p)
        Nbins =  (bootstrap) ? 1.0 : size(Greens_tuple[1],1)
        
        W = 0.5 * Nbins ./ abs.(F.values[mask] .* Nsteps)
        

    else
        # Diagonal error method
        Nsteps = size(Greens_tuple[1],1)
        W = 0.5 ./ real.(Greens_tuple[2] .* conj.(Greens_tuple[2]) .* Nsteps)
        Kp = K
        corr_avg_p = Greens_tuple[1]
    end
    return W, Kp, corr_avg_p

end

function normsq(val)
    return @. val * conj(val)
end

# eliminate ≈0.0 eigenvalues from diagonalized covariance matrix which arise due to symmetries
function get_covariance_mask(params)
    grid = params.input_grid
    mask = Array{Bool}(undef,size(grid,1))
    mask .= true

    # Fermions G(τ)+G(0) = 1.0
    if params.kernel_type == "time_fermionic"
        if findfirst(==(params.β),params.input_grid) != nothing
            mask = .!(params.input_grid .≈ 0.0)
        end
    # Bosons G(τ) = G(β-τ)
    elseif params.kernel_type == "time_bosonic_symmetric"
        β2 = params.β / 2.0
        for (τi,τ) in enumerate(params.input_grid)
            if τ < β2
                mask[τi] = !(findfirst(≈(τ),params.input_grid) != nothing)
            end
        end
    end
    return mask

end