
# Integration weights for a grid whose points represent the centers of energy bins.
# The first and last bins are extrapolated by half of the adjacent grid spacing.
# Consequently, an evenly spaced grid retains the historical weight Δω at
# every point, while an uneven grid uses the local bin width.
function omega_weights(out_ωs::AbstractVector{T}) where {T<:Real}
    nω = length(out_ωs)
    nω >= 2 || throw(ArgumentError("out_ωs must contain at least two points"))

    Δω = diff(out_ωs)
    all(>(zero(T)), Δω) || throw(ArgumentError("out_ωs must be strictly increasing"))

    weights = similar(out_ωs, float(T), nω)
    weights[1] = Δω[1]
    @inbounds for i in 2:nω-1
        weights[i] = (Δω[i-1] + Δω[i]) / 2
    end
    weights[end] = Δω[end]
    return weights
end

# Generate Kernel times the integration weights for the energy grid.
# Notably, for bosonic kernels we multiply in a factor of ω in the n_b routine.
# This makes the kernel and the spectral function positive and analytic for all ω
# Before we return data, though, we will multiply bosonic functions by ω
function generate_K(params::DEACParameters)
    nω = size(params.out_ωs,1)
    ngrid = size(params.input_grid,1)
    if occursin("frequency_fermionic",params.kernel_type)  || occursin("frequency_bosonic", params.kernel_type)
        K = zeros(ComplexF64,(ngrid,nω))
    else
        K = zeros(Float64,(ngrid,nω))
    end
    

    ω_weights = omega_weights(params.out_ωs)


    if params.kernel_type == "time_bosonic"
        
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] = ω_weights[ω]*exp(-params.out_ωs[ω]*params.input_grid[τ])
            end
        end
    elseif params.kernel_type == "time_bosonic_symmetric" 
        
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] =   ω_weights[ω]*(exp(-params.out_ωs[ω]*params.input_grid[τ]) + exp(-params.out_ωs[ω]*(params.β - params.input_grid[τ])))
            end
        end
    elseif params.kernel_type == "time_bosonic_symmetric_w"
        nb = n_b(params)
        if params.input_grid[end] > params.β /2
            coeff = 0.5
        else
            coeff = 0.25
        end
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] =    coeff * ω_weights[ω]*(exp(-params.out_ωs[ω]*params.input_grid[τ]) + exp(-params.out_ωs[ω]*(params.β - params.input_grid[τ]))) * nb[ω]
            end
        end
    elseif params.kernel_type == "time_fermionic"
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] = ω_weights[ω] / (exp(params.out_ωs[ω] * params.input_grid[τ]) + exp(-params.out_ωs[ω] * (params.β - params.input_grid[τ])))
            end
        end
    elseif params.kernel_type == "frequency_fermionic"
        for (iω, ω) in enumerate(params.out_ωs)
            for (iωn, ωn) in enumerate(params.input_grid)
                K[iωn,iω] = - ω_weights[iω] / (1im*ωn - ω)  #(1im * ωn + ω )/ ( ωn^2 + ω^2 )
            end
        end
    elseif params.kernel_type == "frequency_bosonic"
        for (iω, ω) in enumerate(params.out_ωs)
            for (iωn, ωn) in enumerate(params.input_grid)
                if ωn == 0.0 && ω ≈ 0.0
                    K[iωn,iω] = -1.0 * ω_weights[iω]
                else
                    K[iωn,iω] =  ω_weights[iω] * ω / (1im * ωn - ω )
                end
            end
        end
    elseif params.kernel_type == "frequency_bosonic_symmetric"
        for (iω, ω) in enumerate(params.out_ωs)
            for (iωn, ωn) in enumerate(params.input_grid)
                if ωn == 0.0 && ω ≈ 0.0
                    K[iωn,iω] = 2.0 * ω_weights[iω]
                else
                    K[iωn,iω] =  2.0 * ω_weights[iω] * ω * ω / ( ωn^2 + ω^2 )
                end
            end
        end
    elseif params.kernel_type == "time_bosonic_no_bose_einstein"
        if params.input_grid[end] > params.β /2
            coeff = 1.0
        else
            coeff = 0.5
        end
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] =     coeff*ω_weights[ω]*(exp(-params.out_ωs[ω]*params.input_grid[τ]) + exp(-params.out_ωs[ω]*(params.β - params.input_grid[τ])))
            end
        end
    end # kernel_type

    

    return K
end # generate_K()

# calculate Bose factor * ω, 
function n_b(params::DEACParameters)
    close = 1.0e-6
    nω = size(params.out_ωs,1)
    arr = zeros(Float64,nω)
    for ω in 1:nω
        # L'hopital
        if abs(params.out_ωs[ω]) < close
            arr[ω] = 1.0 /params.β
        else
            arr[ω] = params.out_ωs[ω] / (1.0 - exp(-params.β * params.out_ωs[ω]))
        end
    end 
    return arr
end # n_b()
