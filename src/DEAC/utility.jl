
# χ² fit for real correlation functions
function Χ²(observed::AbstractVector{T},calculated::AbstractMatrix{T},W::AbstractVector{T}) where {T<:Real}
    Χ = zeros(T,(size(calculated,2),))

    for pop in 1:size(calculated,2)
        @inbounds Χ[pop] = sum( ((observed .- calculated[:,pop]).^2) .* W )
    end
# println(Χ);  exit()
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
    
    t_linAlg = @b mul!($C,$A,$B)
    t_avx = @b gemmSIMD!($C,$A,$B)
    use_SIMD = t_avx.time < t_linAlg.time
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
function bin_results!(bin_data,calculated_zeroth_moment,run_data,weight_data,curbin,Δω,Greens_tuple,fitness,seed_vec,generations,params)
    for fit_idx ∈ 1:size(fitness,1)
        bin_data[:,curbin,fit_idx] = run_data[:,curbin,fit_idx] ./  sum(weight_data[curbin,fit_idx])
        calculated_zeroth_moment[1,curbin,fit_idx] = sum(bin_data[:,curbin,fit_idx]) .* Δω 
    
        # Bosonic time kernels steal a factor of ω from the spectral function.
        # Multiply it back in if needed
    
        if  occursin("bosonic",params.kernel_type)
            bin_data[:,curbin,fit_idx] = @. 0.5 * bin_data[:,curbin,fit_idx] * (1- exp(-params.β * params.out_ωs))
        end
    end    
    
    if curbin != params.num_bins
        save_checkpoint(bin_data,generations,curbin,params,Greens_tuple,calculated_zeroth_moment,fitness,seed_vec)
    end
    
    println("Finished bin ",curbin," of ",params.num_bins)
    
end

# Calculate matrices used to go from ω to τ space and χ² fit
function calculate_fit_matrices(Greens_tuple,K,use_SIMD,bootstrap,params,eigenvalue_ratio_min)

    if Greens_tuple[2] == nothing
        # Covariance Methods
              
        
        mask = get_covariance_mask(params)


        # datas = Greens_tuple[1]# * U1
        corr_avg = Statistics.mean(Greens_tuple[1],dims=1)
        # Find eigenbasis for covariance matrix 
        cov_matrix = Statistics.cov(Greens_tuple[1],dims=1,corrected=true)
        F = eigen(cov_matrix)
        max_eig = maximum(F.values)
        # println(F.values)
        # catch near-zero eigenvalues not found using get_covariance_mask
        mask = (F.values .> max_eig * eigenvalue_ratio_min) .&& mask 
        # println(mask)
        println(size(F.vectors))
        
        U = F.vectors#[:,mask]
        corr_avg_p = corr_avg * U
        
        # Rotate correlation functions
        GEMM!(corr_avg_p,corr_avg,U,use_SIMD)
        # corr_avg_p = corr_avg_p * U1
        K = transpose(U) * K
        corr_avg_p = corr_avg_p[1,:]
        # exit()
        

        ################################
        S = svd(K)
        println(size(K), " K")
        svd_ratio = 1e-14
        println(S.S)
        # exit()
        
        minval = maximum(S.S) * svd_ratio
        S2 = S.S[S.S .> minval]
        nS2 = size(S2,1)
        U1 = S.U[:,1:nS2] 
        
        println(nS2, " nS2")
        Kp =  Diagonal(S2) * S.Vt[1:nS2,:]
        # println("\n",K2[1,:])
        # println("\n",K[1,:])
        corr_avg_p = (U1' *corr_avg_p)[:,1]
        println(size(corr_avg_p))
        # println(size(U1' * U1))
        # println(size((S.Vt[1:nS2,:])'))
        # println(size(Greens_tuple[1]))
        # println(size(Greens_tuple[1] * U1))
        # println(diag(U1' * U1),"\n")
        # println((U1' * U1)[2,:],"\n")
        # exit()
        err = U1' * F.values#[mask] 
        
        #############################
        println(err)

        Nsteps = size(corr_avg_p)
        Nbins =  (bootstrap) ? 1.0 : size(err,1)
        
        W = 0.5 .* Nbins ./ abs.(err .* Nsteps)
        full_eigen = F.values
        # println(F.values)
        # exit()

    else
        # Diagonal error method
        S = svd(K)
        println(size(K))
        svd_ratio = 1e-14
        println(S.S)
        # exit()
        println(size(S.U))
        println(size(S.S))
        println(size(S.Vt))
        minval = maximum(S.S) * svd_ratio
        S2 = S.S[S.S .> minval]
        nS2 = size(S2,1)
        U1 = S.U[:,1:nS2] 
        
        println(nS2, " nS2")
        K2 =  Diagonal(S2) * S.Vt[1:nS2,:]
        # println("\n",K2[1,:])
        # println("\n",K[1,:])
        
        # println(size(U1' * U1))
        # println(size((S.Vt[1:nS2,:])'))
        # println(size(Greens_tuple[1]))
        # println(size(Greens_tuple[1] * U1))
        # println(diag(U1' * U1),"\n")
        # println((U1' * U1)[2,:],"\n")
        # exit()
        err = (reshape(Greens_tuple[2],(1,:))*U1)[1,:]
        # println(size(err))
        # exit()

        Nsteps = size(Greens_tuple[1],1)
        W = 0.5 ./ real.(err .* conj.(err) .* Nsteps)
        Kp = K2
        corr_avg_p = (reshape(Greens_tuple[1],(1,:))*U1)[1,:]
        full_eigen = nothing
    end
    return W, Kp, corr_avg_p, full_eigen

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