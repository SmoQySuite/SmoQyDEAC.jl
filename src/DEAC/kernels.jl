
# Generate Kernel * Δω
# Notably, for bosonic kernels we multiply in a factor of ω in the n_b routine.
# This makes the kernel and the spectral function positive and analytic for all ω
# Before we return data, though, we will multiply bosonic functions by ω
function generate_K(params::DEACParameters)
    nω = size(params.out_ωs,1)
    ngrid = size(params.input_grid,1)
    if occursin("frequency",params.kernel_type)
        K = zeros(ComplexF64,(ngrid,nω))
    else
        K = zeros(Float64,(ngrid,nω))
    end

    Δω = (params.out_ωs[end]-params.out_ωs[1])/(size(params.out_ωs,1)-1)


    if params.kernel_type == "time_bosonic"
        nb =  n_b(params)
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] = Δω*exp(-params.out_ωs[ω]*params.input_grid[τ]) * nb[ω]
            end
        end
    elseif params.kernel_type == "time_bosonic_symmetric"
        nb = n_b(params)
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] = Δω*(exp(-params.out_ωs[ω]*params.input_grid[τ]) + exp(-params.out_ωs[ω]*(params.β - params.input_grid[τ]))) * nb[ω]
            end
        end
    elseif params.kernel_type == "time_fermionic"
        for ω in 1:nω
            for τ in 1:ngrid
                K[τ,ω] = Δω / (exp(params.out_ωs[ω] * params.input_grid[τ]) + exp(-params.out_ωs[ω] * (params.β - params.input_grid[τ])))
            end
        end
    elseif params.kernel_type == "frequency_fermionic"
        for (iω, ω) in enumerate(params.out_ωs)
            for (iωn, ωn) in enumerate(params.input_grid)
                K[iωn,iω] = Δω / (1im * ωn - ω )#(exp(params.out_ωs[ω] * params.input_grid[τ]) + exp(-params.out_ωs[ω] * (params.β - params.input_grid[τ])))
            end
        end
    elseif params.kernel_type == "frequency_bosonic"
        for (iω, ω) in enumerate(params.out_ωs)
            for (iωn, ωn) in enumerate(params.input_grid)
                if ωn == 0.0 && ω ≈ 0.0
                    K[iωn,iω] = -1.0 * Δω
                else
                    K[iωn,iω] = Δω * ω / (1im * ωn - ω )#(exp(params.out_ωs[ω] * params.input_grid[τ]) + exp(-params.out_ωs[ω] * (params.β - params.input_grid[τ])))
                end
            end
        end
    elseif params.kernel_type == "frequency_bosonic_symmetric"
        for (iω, ω) in enumerate(params.out_ωs)
            for (iωn, ωn) in enumerate(params.input_grid)
                if ωn == 0.0 && ω ≈ 0.0
                    K[iωn,iω] = -2.0*Δω 
                else
                    K[iωn,iω] = -2.0 * Δω * ω * ω / ( ωn^2 + ω^2 )#(exp(params.out_ωs[ω] * params.input_grid[τ]) + exp(-params.out_ωs[ω] * (params.β - params.input_grid[τ])))
                end
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
