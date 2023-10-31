@doc raw"""
    DEAC_Std(correlation_function::AbstractVector,
         correlation_function_error::AbstractVector,
         β::Float64,
         input_grid::Vector{Float64},
         out_ωs::Vector{Float64},
         kernel_type::String,
         num_bins::Int64,
         runs_per_bin::Int64,
         output_file::String,
         checkpoint_directory::String;

         find_ideal_fitness::Bool=true
         population_size::Int64=8,
         base_seed::Integer=8675309,
         number_of_generations::Int64=1000000,
         keep_bin_data=true,
         autoresume_from_checkpoint=false,
         verbose::Bool=false,
         
         stop_minimum_fitness::Float64=1.0,
         crossover_probability::Float64=0.9,
         self_adapting_crossover_probability::Float64=0.1,
         differential_weight::Float64=0.9,
         self_adapting_differential_weight_probability::Float64=0.1,
         self_adapting_differential_weight::Float64=0.9)
    
Runs the DEAC algorithm on data passed in `correlation_function` using $\Chi^2$ fitting from the error passed in by `correlation_function_error`.
# Arguments
- `correlation_function::AbstractVector`: Input data in τ space 
- `correlation_function_error::AbstractVector`: Error associated with input data
- `β::Float64`: Inverse temperature
- `input_grid::Vector{Float64}`: Evenly spaced values in τ from 0 to β, including end points
- `out_ωs::Vector{Float64}`: Energies for AC output
- `kernel_type::String`: See below for allowable kernels
- `num_bins::Int64`: Bins to generate
- `runs_per_bin::Int64`: Number of runs per bin for statistics
- `output_file::String`: File to store output dictionary. jld2 format recommended
- `checkpoint_directory::String`: Directory to store checkpoint data. 

# Optional Arguments
- `find_ideal_fitness::Bool=true`: Use ideal fitness finder
- `population_size::Int64=8`: DEAC population size. Must be ≥ 6
- `base_seed::Int64=8675309`: Seed
- `number_of_generations::Int64=100000`: Maximum number of mutation loops
- `keep_bin_data::Bool=true`: Save binned data or not
- `autoresume_from_checkpoint::Bool=false`: Resume from checkpoint if possible
- `verbose::Bool=false`: Print stats per run
- `stop_minimum_fitness::Float64=1.0`: Value below which fit is considered good, only applies if `find_ideal_fitness=false`

# Optional algorithm arguments
- `crossover_probability::Float64=0.9`: Starting likelihood of crossover
- `self_adapting_crossover_probability::Float64=0.1`: Chance of crossover probability changing
- `differential_weight::Float64=0.9`: Weight for second and third mutable indices
- `self_adapting_differential_weight_probability::Float64=0.1`: Likelihood of SAD changing
- `self_adapting_differential_weight::Float64=0.9`: SAD

Each run will use its own seed. E.g. if you run 10 bins with 100 runs per bin, you will use seeds `base_seed:base_seed+999`. 
You may increment your base seed by 1000, use another output file name, and generate more statistics later.

# Output
SmoQyDEAC will save a dictionary to the file name passed via the `output_file` argument. The same dictionary will be returned by the function.

# Checkpointing
SmoQyDEAC will place a file named `DEAC_checkpoint.jld2` in the directory passed in `checkpoint_directory`. After completing every bin this file will be saved.
After the last bin is finished the file will be deleted. When `autoresume_from_checkpoint=true` SmoQyDEAC will attempt to resume from the checkpoint. 
If the arguments passed do not match those in the checkpoint the code will exit.
"""
function DEAC_Std(correlation_function::AbstractVector,
                  correlation_function_error::AbstractVector,
                  β::Float64,
                  input_grid::Vector{Float64},
                  out_ωs::Vector{Float64},
                  kernel_type::String,
                  num_bins::Int64,
                  runs_per_bin::Int64,
                  output_file::String,
                  checkpoint_directory::String;

                  population_size::Int64=8,
                  base_seed::Integer=8675309,
                  crossover_probability::Float64=0.9,
                  self_adapting_crossover_probability::Float64=0.1,
                  differential_weight::Float64=0.9,
                  self_adapting_differential_weight_probability::Float64=0.1,
                  self_adapting_differential_weight::Float64=0.9,
                  stop_minimum_fitness::Float64=1.0,
                  number_of_generations::Int64=1000000,
                  autoresume_from_checkpoint=false,
                  keep_bin_data=true,
                  W_ratio_max = 1.0e6,
                  find_ideal_fitness::Bool=true,
                  verbose::Bool=false,
                )
    #
    println("\n*** It is highly recommended to use binned data and the covariant matrix method instead (DEAC_Binned) if possible ***\n")
    params = DEACParameters(β,input_grid,out_ωs,kernel_type,output_file,checkpoint_directory,
                            num_bins,runs_per_bin,population_size,base_seed,
                            crossover_probability,self_adapting_crossover_probability,
                            differential_weight,self_adapting_differential_weight_probability,
                            self_adapting_differential_weight,stop_minimum_fitness,number_of_generations)
    #
    return run_DEAC((correlation_function,correlation_function_error),params,autoresume_from_checkpoint,keep_bin_data,W_ratio_max,find_ideal_fitness,verbose)
end # DEAC_Std


@doc raw"""
    DEAC_Binned(correlation_function::AbstractMatrix,
         β::Float64,
         input_grid::Vector{Float64},
         out_ωs::Vector{Float64},
         kernel_type::String,
         num_bins::Int64,
         runs_per_bin::Int64,
         output_file::String,
         checkpoint_directory::String;

         find_ideal_fitness::Bool=true
         population_size::Int64=8,
         base_seed::Integer=8675309,
         number_of_generations::Int64=1000000,
         keep_bin_data=true,
         autoresume_from_checkpoint=false,
         verbose::Bool=false,
         stop_minimum_fitness::Float64=1.0,
         
         crossover_probability::Float64=0.9,
         self_adapting_crossover_probability::Float64=0.1,
         differential_weight::Float64=0.9,
         self_adapting_differential_weight_probability::Float64=0.1,
         self_adapting_differential_weight::Float64=0.9,
         )
    
Runs the DEAC algorithm on data passed in `correlation_function` using $\Chi^2$ fitting from the error passed in by `correlation_function_error`.
# Arguments
- `correlation_function::AbstractVector`: Input data in τ space 
- `correlation_function_error::AbstractVector`: Error associated with input data
- `β::Float64`: Inverse temperature
- `input_grid::Vector{Float64}`: Evenly spaced values in τ from 0 to β, including end points
- `out_ωs::Vector{Float64}`: Energies for AC output
- `kernel_type::String`: See below for allowable kernels
- `num_bins::Int64`: Bins to generate
- `runs_per_bin::Int64`: Number of runs per bin for statistics
- `output_file::String`: File to store output dictionary. jld2 format recommended
- `checkpoint_directory::String`: Directory to store checkpoint data. 

# Optional Arguments
- `find_ideal_fitness::Bool=true`: Use ideal fitness finder
- `population_size::Int64=8`: DEAC population size. Must be ≥ 6
- `base_seed::Int64=8675309`: Seed
- `number_of_generations::Int64=100000`: Maximum number of mutation loops
- `keep_bin_data::Bool=true`: Save binned data or not
- `autoresume_from_checkpoint::Bool=false`: Resume from checkpoint if possible
- `verbose::Bool=false`: Print stats per run
- `stop_minimum_fitness::Float64=1.0`: Value below which fit is considered good, only applies if `find_ideal_fitness=false`

# Optional algorithm arguments
- `crossover_probability::Float64=0.9`: Starting likelihood of crossover
- `self_adapting_crossover_probability::Float64=0.1`: Chance of crossover probability changing
- `differential_weight::Float64=0.9`: Weight for second and third mutable indices
- `self_adapting_differential_weight_probability::Float64=0.1`: Likelihood of SAD changing
- `self_adapting_differential_weight::Float64=0.9`: SAD
- `W_ratio_max::Float64=1.0e6`: Χ² ~ 1.0/σ², this parameter prevents [near] singularities for very small σ 
- `bootstrap_bins::Int=0`: The algorithm requires more bins than τ steps. We use bootstrapping to get 5 * nτ bins by default. User may set this higher 
                  

Each run will use its own seed. E.g. if you run 10 bins with 100 runs per bin, you will use seeds `base_seed:base_seed+999`. 
You may increment your base seed by 1000, use another output file name, and generate more statistics later.

# Output
SmoQyDEAC will save a dictionary to the file name passed via the `output_file` argument. The same dictionary will be returned by the function.

# Checkpointing
SmoQyDEAC will place a file named `DEAC_checkpoint.jld2` in the directory passed in `checkpoint_directory`. After completing every bin this file will be saved.
After the last bin is finished the file will be deleted. When `autoresume_from_checkpoint=true` SmoQyDEAC will attempt to resume from the checkpoint. 
If the arguments passed do not match those in the checkpoint the code will exit. 
"""
function DEAC_Binned(correlation_function::AbstractMatrix,
                  β::Float64,
                  input_grid::Vector{Float64},
                  out_ωs::Vector{Float64},
                  kernel_type::String,
                  num_bins::Int64,
                  runs_per_bin::Int64,
                  output_file::String,
                  checkpoint_directory::String;

                  population_size::Int64=8,
                  base_seed::Integer=8675309,
                  crossover_probability::Float64=0.9,
                  self_adapting_crossover_probability::Float64=0.1,
                  differential_weight::Float64=0.9,
                  self_adapting_differential_weight_probability::Float64=0.1,
                  self_adapting_differential_weight::Float64=0.9,
                  stop_minimum_fitness::Float64=1.0,
                  number_of_generations::Int64=1000000,
                  autoresume_from_checkpoint::Bool=false,
                  keep_bin_data::Bool=true,
                  W_ratio_max::Float64 = 1.0e6,
                  bootstrap_bins::Int = 0,
                  find_ideal_fitness::Bool = true,
                  verbose::Bool=false,
                )
    #
    # Bootstrap bins to ensure sufficient bin size
    if bootstrap_bins ≤ 0 || (size(correlation_function,1) < 5 * size(correlation_function,2))
        bootstrap_bins = max(bootstrap_bins,5*size(correlation_function,2))
        correlation_function = bootstrap_samples(correlation_function,bootstrap_bins,base_seed )
    end


    params = DEACParameters(β,input_grid,out_ωs,kernel_type,output_file,checkpoint_directory,
                            num_bins,runs_per_bin,population_size,base_seed,
                            crossover_probability,self_adapting_crossover_probability,
                            differential_weight,self_adapting_differential_weight_probability,
                            self_adapting_differential_weight,stop_minimum_fitness,number_of_generations)
    #
    return run_DEAC((correlation_function,nothing),params,autoresume_from_checkpoint,keep_bin_data,W_ratio_max,find_ideal_fitness,verbose)
end # DEAC_Binned()

# Run the DEAC algorithm
function run_DEAC(Greens_tuple,
                  params::DEACParameters,
                  autoresume_from_checkpoint::Bool,
                  keep_bin_data::Bool,
                  W_ratio_max::Float64,
                  find_ideal_fitness::Bool,
                  verbose::Bool)
    
    # Assert parameters are within allowable/realistic ranges
    @assert params.population_size >= 6 # DEAC can be run with as few as 4, but it gives garbage results
    @assert params.β > 0.0
    @assert params.num_bins > 1
    @assert params.runs_per_bin >= 1
    @assert params.kernel_type in allowable_kernels
    @assert params.crossover_probability > 0.0 && params.crossover_probability < 1.0
    @assert params.self_adapting_crossover_probability > 0.0 && params.self_adapting_crossover_probability < 1.0
    @assert params.differential_weight > 0.0 && params.differential_weight < 1.0
    @assert params.self_adapting_differential_weight > 0.0 && params.self_adapting_differential_weight < 1.0
    @assert params.self_adapting_differential_weight_probability > 0.0 && params.self_adapting_differential_weight_probability < 1.0
    @assert params.stop_minimum_fitness > 0.0
    @assert params.number_of_generations >= 1
    @assert params.base_seed >= 1

    # Declare variables, set defaults, pre-calculations
    bin_data = zeros(Float64,(size(params.out_ωs,1),params.num_bins))
    generations = zeros(UInt64,params.num_bins)
    run_data = zeros(Float64,(size(params.out_ωs,1),params.num_bins))
    calculated_zeroth_moment = zeros(Float64,(1,params.num_bins))
    
    
    correlation_function = Greens_tuple[1]
    start_bin = 1
    true_fitness = params.stop_minimum_fitness
    total_runs = params.runs_per_bin*params.num_bins
    seed_vec = collect(params.base_seed:params.base_seed +(total_runs -1))
    Δω = (params.out_ωs[end]-params.out_ωs[1])/(size(params.out_ωs,1)-1)
    nthreads = Threads.nthreads()
        
    # Utilize the correct kernel
    K = generate_K(params)
    
    # prevent race conditions when threads finish
    thread_lock = ReentrantLock()

    # Test to see if SIMD instructions via LoopVectorization are working
    if eltype(Greens_tuple[1]) <: Real
        use_SIMD = useSIMD(size(params.out_ωs,1),size(params.input_grid,1),params.population_size)
    else
        use_SIMD = false
    end
    
    # Load from checkpoint if applicable
    if autoresume_from_checkpoint
        chk_exists, chk_dict = find_checkpoint(params)

        if chk_exists
            if compare_checkpoint(chk_dict,params,Greens_tuple)
                println("Checkpoint found at ",joinpath(params.checkpoint_directory,"DEAC_checkpoint.jld2"))
                start_bin = chk_dict["bin_num"] + 1
                println("Parameters match. Resuming at bin ",start_bin,"\n")
                bin_data = chk_dict["bin_data"]
                generations = chk_dict["generations"]
                seed_vec = chk_dict["seeds"]
                calculated_zeroth_moment = chk_dict["zeroth"]
                true_fitness = chk_dict["true_fitness"]
            else
                println("Checkpoint found at ",joinpath(params.checkpoint_directory,"DEAC_checkpoint.jld2"))
                println("Mismatched parameters. Exiting")
                exit()
            end
        end # if file exists
    end # if checkpoint

    start_thread = (start_bin-1) * params.runs_per_bin +1

    # Matrices for calculating fit
    W, Kp, corr_avg_p = calculate_fit_matrices(Greens_tuple,K,W_ratio_max,use_SIMD)
    

    
    
    ### Find Ideal Fitness
    # If, for a thread, over the course of two consecutive fit_check_frequency generations
    # there is ≤ 10% improvement that thread's current fitness goes into fitness array
    # Ideal fit is fit_mod * minimum(fitness)
    if find_ideal_fitness && start_bin == 1
        println("Finding Ideal Fitness Parameter")
        fit_check_frequency = 10000
        fit_check_difference = 0.1
        fit_mod = 1.025
        nfinder = max(nthreads,10)
        fitness = zeros(Float64,nfinder)

        Δt = @elapsed begin
            Threads.@threads for thd in 1:nfinder
                
                # Setup RNG
                seed = params.base_seed + thd
                rng = Random.Xoshiro(seed)

                # Allocate arrays, set random initial state
                population_old  = reshape(Random.rand(rng,size(params.out_ωs,1)*params.population_size),(size(params.out_ωs,1),params.population_size))
                population_new = zeros(Float64,(size(params.out_ωs,1),params.population_size))
                model = zeros(eltype(Greens_tuple[1]),(size(Kp,1),size(population_old,2)))

                crossover_probability_new = zeros(Float64,params.population_size)
                crossover_probability_old = zeros(Float64,params.population_size)
                crossover_probability_old .= params.crossover_probability

                differential_weights_new = zeros(Float64,params.population_size)
                differential_weights_old = zeros(Float64,params.population_size)
                differential_weights_old .= params.differential_weight

                mutate_indices = Array{Float64}(undef,(size(params.out_ωs,1),params.population_size))

                fitness_old = initial_fit!(population_old,model,corr_avg_p,Kp,W,params,use_SIMD)
                
                # IFF tracker variables
                lastDelta = typemax(Float64)
                lastFitness = minimum(fitness_old)
                last2Fitness = typemax(Float64)
                Delta = typemax(Float64)

                # main loop
                for gen in 1:params.number_of_generations
                    
                    update_weights!(crossover_probability_new,differential_weights_new,crossover_probability_old,differential_weights_old,rng,params)

                    # Randomly set some ω points to 'mutate'
                    rand_mutate_array!(mutate_indices,crossover_probability_new,rng,params)

                    # Set triplet of other populations for mutations
                    mutant_indices = get_mutant_indices(rng,params.population_size)
                    
                    # if mutate_indices, do mutation, else keep same
                    propose_populations!(population_new,population_old,mutate_indices,differential_weights_new,mutant_indices,params)
                    
                    # calculate new fitness
                    GEMM!(model,Kp,population_new,use_SIMD)
                    fitness_new = Χ²(corr_avg_p,model,W) ./ size(params.input_grid,1)

                    # update populations if fitness improved
                    update_populations!(fitness_old,crossover_probability_old,differential_weights_old,population_old,fitness_new,crossover_probability_new,differential_weights_new,population_new)

                    # check for low improvement for two consecutive slices, break if so
                    if (gen % fit_check_frequency) == 0
                        
                        last2Fitness = lastFitness
                        lastFitness = minimum(fitness_new)
                        lastDelta = Delta
                        Delta = (last2Fitness-lastFitness)/lastFitness
                        fitness[thd] = lastFitness
                        if (Delta < fit_check_difference) && (lastDelta < fit_check_difference)
                            break
                        end
                    end # if fit_check_frequency
                end # generations
            end #threads
        end  # Δt   
        # Calculate ideal Fitness, print 
        true_fitness = minimum(fitness) * fit_mod
        
        println(@sprintf("\nFitness found in %01.5fs",Δt))
        println(@sprintf("Using Ideal Fitness:  %01.5f\n",true_fitness))
    end # IFF

    finished_runs = (start_bin -1) * params.runs_per_bin
    Δt = @elapsed begin
        # loop over bins*runs_per_bin
        Threads.@threads for thd in start_thread:total_runs

            # Each run utilizes its own RNG with a unique seed
            # println("seed_vec ",1+thd-start_thread)
            seed = seed_vec[1+thd-start_thread]
            rng = Random.Xoshiro(seed)
            # Set initial parameters for algo
            crossover_probability_new = zeros(Float64,params.population_size)
            crossover_probability_old = zeros(Float64,params.population_size)
            crossover_probability_old .= params.crossover_probability

            differential_weights_new = zeros(Float64,params.population_size)
            differential_weights_old = zeros(Float64,params.population_size)
            differential_weights_old .= params.differential_weight

            # Randomly set initial populations, initialize arrays
            population_old  = reshape(Random.rand(rng,size(params.out_ωs,1)*params.population_size),(size(params.out_ωs,1),params.population_size))
            population_new = zeros(Float64,(size(params.out_ωs,1),params.population_size))
            model = zeros(eltype(Greens_tuple[1]),(size(Kp,1),size(population_old,2)))
            mutate_indices = Array{Float64}(undef,(size(params.out_ωs,1),params.population_size))

            fitness_old = initial_fit!(population_old,model,corr_avg_p,Kp,W,params,use_SIMD)
          
            # track number of generations to get to fitness
            numgen = 0


            # Loop over generations until number_of_generations or fitness is achieved
            for gen in 1:params.number_of_generations
            
                # If fitness achieved, exit loop
                if (minimum(fitness_old) <= true_fitness) break; end
                    
                # Modify DEAC parameters stochastically
                update_weights!(crossover_probability_new,differential_weights_new,crossover_probability_old,differential_weights_old,rng,params)


                # Randomly set some ω points to 'mutate'
                rand_mutate_array!(mutate_indices,crossover_probability_new,rng,params)

                # Set triplet of other populations for mutations
                mutant_indices = get_mutant_indices(rng,params.population_size)
                
                # if mutate_indices, do mutation, else keep same
                propose_populations!(population_new,population_old,mutate_indices,differential_weights_new,mutant_indices,params)
                
                # get new fitness
                GEMM!(model,Kp,population_new,use_SIMD)
                fitness_new = Χ²(corr_avg_p,model,W) ./ size(params.input_grid,1)

                # if improved do updates
                update_populations!(fitness_old,crossover_probability_old,differential_weights_old,population_old,fitness_new,crossover_probability_new,differential_weights_new,population_new)

                
                numgen = numgen + 1
            end # generations
            fit, fit_idx = findmin(fitness_old)
            
            # lock to prevent race conditions
            lock(thread_lock) do 
                
                thisrun = 1 + finished_runs % params.runs_per_bin
                curbin = (1 + finished_runs ÷ params.runs_per_bin)
                
                finished_runs += 1
                
                run_data[:,curbin] += population_old[:,fit_idx]
                generations[curbin] += numgen

                # setting seed to 0 allows culling of used seeds in the checkpoint
                # ensuring each seed is used once
                seed_vec[1+thd-start_thread] = 0

                if verbose
                    println(@sprintf("  Bin %3u | Run %4u | Fitness %8.6f | Generations %u",curbin,thisrun,fit,numgen))
                end

                # Calculate bin data if enough to finish a bin, then checkpoint
                if finished_runs % params.runs_per_bin == 0

                    bin_results!(bin_data,calculated_zeroth_moment,run_data,curbin,Δω,Greens_tuple,true_fitness,seed_vec,generations,params)

                end # bin completing
            end # lock(thread_lock)
            

        end # threads
    end
    # time stats
    
    t_per_run = Δt/(total_runs - (start_bin-1)*params.runs_per_bin)
    t_per_bin = Δt/(params.num_bins + 1 - start_bin)

    # Merge data
    zero_avg, zero_err = jackknife(calculated_zeroth_moment)
    gen_per_run = sum(generations ./ params.runs_per_bin)[1]/params.num_bins
    differential = 100.0*abs(1.0-zero_avg[1])
    data, err = jackknife(bin_data)
    
    # Print statistics
    println("\nSaving data to ",params.output_file," and deleting checkpoint file\n")
    
    println("Run Statistics")
    if occursin("fermionic",params.kernel_type)
        println(@sprintf(" Expected 0th moment:   1.00") )
        println(@sprintf(" DEAC 0th moment:       %01.3f ± %01.3f",zero_avg[1],zero_err[1]))
        println(@sprintf(" 0th moment difference: %01.3f%%",differential))
    end
    println(@sprintf(" Mean generations/run:  %01.3f",gen_per_run))
    println(@sprintf(" Total Run time:        %01.3fs",t_per_bin * params.num_bins))
    println(@sprintf(" Run time per bin:      %01.3fs",t_per_bin))
    println(@sprintf(" Run time per genome:   %01.3fs\n",t_per_run*nthreads))
    
    # Create data dictionaries, save to file, and return via function call
    if keep_bin_data
        bin_dict = Dict{String,Any}(
            "A" => data,
            "σ" => err,
            "ωs" => params.out_ωs,
            "zeroth_moment" => zero_avg[1],
            "zeroth_moment_σ" => zero_err[1],
            "avg_generations" => gen_per_run,
            "bin_data" => bin_data,
            "bin_zeroth_moment" => calculated_zeroth_moment,
            "fitness_target" =>  true_fitness,
            "runtime" => Δt
        )
    else
        bin_dict = Dict{String,Any}(
            "A" => data,
            "σ" => err,
            "zeroth_moment" => zero_avg[1],
            "zeroth_moment_σ" => zero_err[1],
            "avg_generations" => gen_per_run,
            "ωs" => params.out_ωs,
            "fitness_target" =>  true_fitness,
            "runtime" => Δt
        )
    end
    FileIO.save(params.output_file,bin_dict)
    delete_checkpoint(params)
    return bin_dict
end # run_DEAC()



