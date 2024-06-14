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
         checkpoint_file::String;

         find_fitness_floor::Bool=true
         population_size::Int64=8,
         base_seed::Integer=8675309,
         number_of_generations::Int64=1000000,
         keep_bin_data=true,
         autoresume_from_checkpoint=false,
         verbose::Bool=false,
         
         fitness = 1.0,
         crossover_probability::Float64 = 0.9,
         self_adapting_crossover_probability::Float64 = 0.1,
         differential_weight::Float64 = 0.9,
         self_adapting_differential_weight_probability::Float64 = 0.1,
         self_adapting_differential_weight::Float64 = 0.9,
         user_mutation! = nothing 
         )
    
Runs the DEAC algorithm on data passed in `correlation_function` using $\Chi^2$ fitting from the error passed in by `correlation_function_error`.
# Arguments
- `correlation_function::AbstractVector`: Input data in τ/ωₙ space 
- `correlation_function_error::AbstractVector`: Error associated with input data
- `β::Float64`: Inverse temperature
- `input_grid::Vector{Float64}`: Evenly spaced values in τ from 0 to β, including end points or ωₙ
- `out_ωs::Vector{Float64}`: Evenly spaced energies for AC output
- `kernel_type::String`: See below for allowable kernels
- `num_bins::Int64`: Bins to generate
- `runs_per_bin::Int64`: Number of runs per bin for statistics
- `output_file::String`: File to store output dictionary. jld2 format recommended
- `checkpoint_file::String`: File to store checkpoint data. jld2 format recommended

# Optional Arguments
- `find_fitness_floor::Bool = false`: Find a reasonable floor for fitness and append that value to target fitnesses.
- `population_size::Int64 = 8`: DEAC population size. Must be ≥ 6
- `base_seed::Int64 = 8675309`: Seed
- `number_of_generations::Int64 = 100000`: Maximum number of mutation loops
- `keep_bin_data::Bool = true`: Save binned data or not
- `autoresume_from_checkpoint::Bool = false`: Resume from checkpoint if possible
- `verbose::Bool = false`: Print stats per run
- `fitness = 1.0`: single value or an array of values for χ² fitting

# Optional algorithm arguments
- `crossover_probability::Float64 = 0.9`: Starting likelihood of crossover
- `self_adapting_crossover_probability::Float64 = 0.1`: Chance of crossover probability changing
- `differential_weight::Float64 = 0.9`: Weight for second and third mutable indices
- `self_adapting_differential_weight_probability::Float64 = 0.1`: Likelihood of SAD changing
- `self_adapting_differential_weight::Float64 = 0.9`: SAD
- `user_mutation! = nothing`: User passed function to add additional mutation to each iteration. See below for more information

Each run will use its own seed. E.g. if you run 10 bins with 100 runs per bin, you will use seeds `base_seed:base_seed+999`. 
You may increment your base seed by 1000, use another output file name, and generate more statistics later.

# Output
SmoQyDEAC will save a dictionary to the file name passed via the `output_file` argument. The same dictionary will be returned by the function.

# Checkpointing
SmoQyDEAC will periodically save a file named according to the parameter passed in`checkpoint_file`. JLD2 format is recommended, e.g. a file with the .jld2 extension. After completing every bin this file will be saved.
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
                  checkpoint_file::String;

                  population_size::Int64=8,
                  base_seed::Integer=8675309,
                  crossover_probability::Float64=0.9,
                  self_adapting_crossover_probability::Float64=0.1,
                  differential_weight::Float64=0.9,
                  self_adapting_differential_weight_probability::Float64=0.1,
                  self_adapting_differential_weight::Float64=0.9,
                  fitness=1.0,
                  number_of_generations::Int64=1000000,
                  autoresume_from_checkpoint=false,
                  keep_bin_data=true,
                  find_fitness_floor::Bool=false,
                  verbose::Bool=false,
                  user_mutation! =nothing
                )
    #
    println("\n*** It is highly recommended to use binned data and the covariant matrix method instead (DEAC_Binned) if possible ***\n")

    if isa(fitness,Number)
        fitness = [fitness,]
    end
    sort!(fitness,rev=true)
    
    params = DEACParameters(β,input_grid,out_ωs,kernel_type,output_file,checkpoint_file,
                            num_bins,runs_per_bin,population_size,base_seed,
                            crossover_probability,self_adapting_crossover_probability,
                            differential_weight,self_adapting_differential_weight_probability,
                            self_adapting_differential_weight,number_of_generations)
    #
    return run_DEAC(
        (correlation_function,correlation_function_error),
        params,
        autoresume_from_checkpoint,
        keep_bin_data,
        find_fitness_floor,
        verbose,
        user_mutation!,
        fitness
    )
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
         checkpoint_file::String;

         find_fitness_floor::Bool = false,
         population_size::Int64 = 8,
         base_seed::Integer = 8675309,
         number_of_generations::Int64 = 1000000,
         keep_bin_data = true,
         autoresume_from_checkpoint = false,
         verbose::Bool = false,
         fitness = 1.0,
         
         crossover_probability::Float64 = 0.9,
         self_adapting_crossover_probability::Float64 = 0.1,
         differential_weight::Float64 = 0.9,
         self_adapting_differential_weight_probability::Float64 = 0.1,
         self_adapting_differential_weight::Float64 = 0.9,
         user_mutation! = nothing,
         eigenvalue_ratio_min::Float64 = 1e-8
         )
    
Runs the DEAC algorithm on data passed in `correlation_function` using $\Chi^2$ fitting using the eigenvalues of the covariance matrix
# Arguments
- `correlation_function::AbstractMatrix`: Input data in τ/ωₙ space, shape [Bins,τ/ωₙ] 
- `β::Float64`: Inverse temperature
- `input_grid::Vector{Float64}`: Evenly spaced values in τ from 0 to β, including end points
- `out_ωs::Vector{Float64}`: Evenly spaced energies for AC output
- `kernel_type::String`: See below for allowable kernels
- `num_bins::Int64`: Bins to generate
- `runs_per_bin::Int64`: Number of runs per bin for statistics
- `output_file::String`: File to store output dictionary. jld2 format recommended
- `checkpoint_file::String`: File to store checkpoint data. jld2 format recommended

# Optional Arguments
- `find_fitness_floor::Bool = false`: Find the floor for fitness and append value to target fitnesses
- `population_size::Int64 = 8`: DEAC population size. Must be ≥ 6
- `base_seed::Int64 = 8675309`: Seed
- `number_of_generations::Int64 = 100000`: Maximum number of mutation loops
- `keep_bin_data::Bool = true`: Save binned data or not
- `autoresume_from_checkpoint::Bool = false`: Resume from checkpoint if possible
- `verbose::Bool = false`: Print stats per run
- `fitness = 1.0`: single value or an array of values for χ² fitting

# Optional algorithm arguments
- `crossover_probability::Float64 = 0.9`: Starting likelihood of crossover
- `self_adapting_crossover_probability::Float64 = 0.1`: Chance of crossover probability changing
- `differential_weight::Float64 = 0.9`: Weight for second and third mutable indices
- `self_adapting_differential_weight_probability::Float64 = 0.1`: Likelihood of SAD changing
- `self_adapting_differential_weight::Float64 = 0.9`: SAD
- `user_mutation! = nothing`: User passed function to add additional mutation to each iteration. See below for more information
- `eigenvalue_ratio_min::Float64 = 1e-8`: Cutoff to mask out eigenvectors with eigenvalues ≈ 0.0 that arise due to symmetries in input data

Each run will use its own seed. E.g. if you run 10 bins with 100 runs per bin, you will use seeds `base_seed:base_seed+999`. 
You may increment your base seed by 1000, use another output file name, and generate more statistics later.

# Output
SmoQyDEAC will save a dictionary to the file name passed via the `output_file` argument. The same dictionary will be returned by the function.

# Checkpointing
SmoQyDEAC will periodically save a file named according to the parameter passed in`checkpoint_file`. JLD2 format is recommended, e.g. a file with the .jld2 extension. After completing every bin this file will be saved.
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
                  checkpoint_file::String;

                  population_size::Int64=20,
                  base_seed::Integer=8675309,
                  crossover_probability::Float64=0.9,
                  self_adapting_crossover_probability::Float64=0.1,
                  differential_weight::Float64=0.9,
                  self_adapting_differential_weight_probability::Float64=0.1,
                  self_adapting_differential_weight::Float64=0.9,
                  fitness=1.0,
                  number_of_generations::Int64=1000000,
                  autoresume_from_checkpoint::Bool = false,
                  keep_bin_data::Bool = true,
                  bootstrap_bins::Int = 0,
                  find_fitness_floor::Bool = false,
                  verbose::Bool=false,
                  user_mutation! =nothing,
                  eigenvalue_ratio_min::Float64 = 1e-8,
                )
    #
    #- `bootstrap_bins::Int = 0`: The algorithm requires more bins than τ steps. We use bootstrapping to get 5 * nτ bins by default. User may set this higher 
    if isa(fitness,Number)
        fitness = [fitness,]
    end
    sort!(fitness,rev=true)


    params = DEACParameters(β,input_grid,out_ωs,kernel_type,output_file,checkpoint_file,
                            num_bins,runs_per_bin,population_size,base_seed,
                            crossover_probability,self_adapting_crossover_probability,
                            differential_weight,self_adapting_differential_weight_probability,
                            self_adapting_differential_weight,number_of_generations)


    # Bootstrap bins to ensure sufficient bin size
    do_bootstrap = false
    if bootstrap_bins > 0 || (size(correlation_function,1) < 5 * size(correlation_function,2)) 
        do_bootstrap = true
        bootstrap_bins = max(bootstrap_bins,5*size(correlation_function,2))
        println("Bootstrapping to ", bootstrap_bins, " samples")
        correlation_function = bootstrap_samples(correlation_function,bootstrap_bins,base_seed )
    end

    
    #
    return run_DEAC(
        (correlation_function,nothing),
        params,autoresume_from_checkpoint,
        keep_bin_data,
        find_fitness_floor,
        verbose,
        user_mutation!,
        fitness
        ;bootstrap=do_bootstrap,
        eigenvalue_ratio_min=eigenvalue_ratio_min


    )
end # DEAC_Binned()

# Run the DEAC algorithm
function run_DEAC(Greens_tuple,
                  params::DEACParameters,
                  autoresume_from_checkpoint::Bool,
                  keep_bin_data::Bool,
                  find_fitness_floor::Bool,
                  verbose::Bool,
                  user_mutation!,
                  fitness;
                  bootstrap=false,
                  eigenvalue_ratio_min=nothing
    )

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
    @assert minimum(fitness) > 0.0 || find_fitness_floor
    @assert params.number_of_generations >= 1
    @assert params.base_seed >= 1

    
    # Declare variables, set defaults, pre-calculations
    bin_data = nothing # zeros(Float64,(size(params.out_ωs,1),params.num_bins))
    generations = nothing # zeros(UInt64,params.num_bins)
    run_data = nothing # zeros(Float64,(size(params.out_ωs,1),params.num_bins))
    weight_data = nothing #zeros(Float64,params.num_bins)
    calculated_zeroth_moment = nothing #zeros(Float64,(1,params.num_bins))
    
    start_bin = 1
    total_runs = params.runs_per_bin*params.num_bins
    seed_vec = collect(params.base_seed:params.base_seed +(total_runs -1))
    Δω = (params.out_ωs[end]-params.out_ωs[1])/(size(params.out_ωs,1)-1)
    nthreads = Threads.nthreads()
        
    # Utilize the correct kernel
    K = generate_K(params)
    normalize = occursin("time",params.kernel_type)
    normK = ones(Float64,(1,size(params.out_ωs,1))) .* Δω
    target_zeroth = 1.0
    if normalize && occursin("bosonic",params.kernel_type)
        if Greens_tuple[2] == nothing
            target_zeroth =  mean(Greens_tuple[1][:,1])
        else
            target_zeroth =  Greens_tuple[1][1]
        end
        if (params.kernel_type == "time_bosonic_symmetric") && normalize
            normK[1,:] =   K[1,:] 
            
        end
    end
    

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
                println("Checkpoint found at ",params.checkpoint_file)
                start_bin = chk_dict["bin_num"] + 1
                println("Parameters match. Resuming at bin ",start_bin,"\n")
                bin_data = chk_dict["bin_data"]
                generations = chk_dict["generations"]
                seed_vec = chk_dict["seeds"]
                calculated_zeroth_moment = chk_dict["zeroth"]
                fitness = chk_dict["fitness"]
            else
                println("Checkpoint found at ",params.checkpoint_file)
                println("Mismatched parameters. Exiting")
                exit()
            end
        end # if file exists
    end # if checkpoint
    
    start_thread = (start_bin-1) * params.runs_per_bin +1

    # Matrices for calculating fit
    W, Kp, corr_avg_p, full_eigen = calculate_fit_matrices(Greens_tuple,K,use_SIMD,bootstrap,params,eigenvalue_ratio_min)
    
    
    
    ### Fitness Floor Finder
    # If, for a thread, over the course of two consecutive fit_check_frequency generations
    # there is ≤ 10% improvement that thread's current fitness goes into fitness array
    # Ideal fit is fit_mod * minimum(fitness)
    if find_fitness_floor && start_bin == 1
        println("\nFinding Fitness Floor")
        fit_check_frequency = 10000
        fit_check_difference = 0.1
        fit_mod = 1.025
        nfinder = max(nthreads,10)
        fitness_floor = zeros(Float64,nfinder)

        Δt = @elapsed begin
            # Threads.@threads 
            for thd in 1:nfinder
                
                # Setup RNG
                seed = params.base_seed + thd
                rng = Random.Xoshiro(seed)

                # Allocate arrays, set random initial state
                population_old  = reshape(Random.rand(rng,size(params.out_ωs,1)*params.population_size),(size(params.out_ωs,1),params.population_size))
                population_new = zeros(Float64,(size(params.out_ωs,1),params.population_size))
                norm_array = zeros(Float64,(1,params.population_size))
                
                model = zeros(eltype(Greens_tuple[1]),(size(Kp,1),size(population_old,2)))

                crossover_probability_new = zeros(Float64,params.population_size)
                crossover_probability_old = zeros(Float64,params.population_size)
                crossover_probability_old .= params.crossover_probability

                differential_weights_new = zeros(Float64,params.population_size)
                differential_weights_old = zeros(Float64,params.population_size)
                differential_weights_old .= params.differential_weight

                mutate_indices = Array{Float64}(undef,(size(params.out_ωs,1),params.population_size))

                fitness_old = initial_fit!(population_old,model,corr_avg_p,Kp,W,params,use_SIMD,normalize,normK,target_zeroth,norm_array)
                
                # FFF tracker variables
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
                    propose_populations!(population_new,population_old,mutate_indices,differential_weights_new,mutant_indices,params,normalize,normK,target_zeroth,norm_array,use_SIMD)
                    
                    # calculate new fitness
                    GEMM!(model,Kp,population_new,use_SIMD)
                    fitness_new = Χ²(corr_avg_p,model,W) #./ size(params.input_grid,1)
                    
                    # update populations if fitness improved
                    update_populations!(fitness_old,crossover_probability_old,differential_weights_old,population_old,fitness_new,crossover_probability_new,differential_weights_new,population_new)

                    # check for low improvement for two consecutive slices, break if so
                    if (gen % fit_check_frequency) == 0
                        
                        last2Fitness = lastFitness
                        lastFitness = minimum(fitness_new)
                        lastDelta = Delta
                        Delta = (last2Fitness-lastFitness)/lastFitness
                        fitness_floor[thd] = lastFitness
                        if (Delta < fit_check_difference) && (lastDelta < fit_check_difference)
                            break
                        end
                    end # if fit_check_frequency
                end # generations
            end #threads
        end  # Δt   
        # Calculate ideal Fitness, print 
        true_fitness = minimum(fitness_floor) * fit_mod
        
        println(@sprintf("Fitness floor found in %01.5fs",Δt))
        println(@sprintf("Adding Fitness Floor to output arrays:  %01.5f",true_fitness))

        # Remove fitness targets below IFF from array and append IFF value
        idx_low = findall(fitness .< true_fitness)
        if idx_low != []
            println("Fitnesses ", fitness[idx_low], " below floor. Removing.")
        end
        idx_high = findall(fitness .> true_fitness)
        old_params = copy(fitness[idx_high])
        veclen = size(old_params,1)
        fitness = zeros(veclen+1)
        fitness[1:veclen] = old_params
        fitness[end] = true_fitness
        
        
    end # IFF
    println("Target fitnesses are: ",fitness,"\n")
        
    if bin_data == nothing
        bin_data = zeros(Float64,(size(params.out_ωs,1),params.num_bins,size(fitness,1)))
        generations = zeros(UInt64,(params.num_bins,size(fitness,1)))
        calculated_zeroth_moment = zeros(Float64,(1,params.num_bins,size(fitness,1)))
        
    end
    
    run_data = zeros(Float64,(size(params.out_ωs,1),params.num_bins,size(fitness,1)))
    weight_data = zeros(Float64,(params.num_bins,size(fitness,1)))
    
    finished_runs = (start_bin -1) * params.runs_per_bin
    Δt = @elapsed begin
        # loop over bins*runs_per_bin
        Threads.@threads for thd in start_thread:total_runs

            # Track fits in case there are multiple fits.
            current_fit_idx = 1
            tmp_generations = zeros(UInt64,size(fitness,1))
            tmp_run_data = zeros(Float64,(size(params.out_ωs,1),size(fitness,1)))
            tmp_fit_data = zeros(Float64,(size(fitness,1)))
        
            # Each run utilizes its own RNG with a unique seed
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
            norm_array = zeros(Float64,(1,params.population_size))
            
            model = zeros(eltype(Greens_tuple[1]),(size(Kp,1),size(population_old,2)))
            mutate_indices = Array{Float64}(undef,(size(params.out_ωs,1),params.population_size))

            fitness_old = initial_fit!(population_old,model,corr_avg_p,Kp,W,params,use_SIMD,normalize,normK,target_zeroth,norm_array)
          
            # track number of generations to get to fitness
            numgen = 0


            # Loop over generations until number_of_generations or fitness is achieved
            for gen in 1:params.number_of_generations
            
                # If fitness achieved, exit loop
                if (minimum(fitness_old) <= fitness[end]) break; end
                    
                # Modify DEAC parameters stochastically
                update_weights!(crossover_probability_new,differential_weights_new,crossover_probability_old,differential_weights_old,rng,params)


                # Randomly set some ω points to 'mutate'
                rand_mutate_array!(mutate_indices,crossover_probability_new,rng,params)

                # Set triplet of other populations for mutations
                mutant_indices = get_mutant_indices(rng,params.population_size)
                
                # if mutate_indices, do mutation, else keep same
                propose_populations!(population_new,population_old,mutate_indices,differential_weights_new,mutant_indices,params,normalize,normK,target_zeroth,norm_array,use_SIMD)
                
                # get new fitness
                GEMM!(model,Kp,population_new,use_SIMD)
                fitness_new = Χ²(corr_avg_p,model,W) #./ size(params.input_grid,1)

                # if improved do updates
                update_populations!(fitness_old,crossover_probability_old,differential_weights_old,population_old,fitness_new,crossover_probability_new,differential_weights_new,population_new)

                # do user mutation if applicable
                if (user_mutation! != nothing)
                    user_mutation!(population_new,population_old,rng)
                    GEMM!(model,Kp,population_new,use_SIMD)
                    fitness_new = Χ²(corr_avg_p,model,W) #./ size(params.input_grid,1)
                    
                    for pop in axes(fitness_old,1)
                        if fitness_new[pop] <= fitness_old[pop]
                            fitness_old[pop] = fitness_new[pop]
                            population_old[:,pop] = population_new[:,pop]
                        end
                    end
                end
                
                # Populate any values in fitness array which are met
                while (minimum(fitness_old) < fitness[current_fit_idx] )
                    tmp_fit_data[current_fit_idx], fit_idx = findmin(fitness_old)
                    tmp_run_data[:,current_fit_idx] = population_old[:,fit_idx]
                    tmp_generations[current_fit_idx] = numgen
                    current_fit_idx += 1
                    if current_fit_idx > size(fitness,1); break; end
                end

                numgen = numgen + 1
            end # generations

            if numgen ≥ params.number_of_generations
                fit_data, fit_idx = findmin(fitness_old)
                cur_best_pop = population_old[:,fit_idx]
                for current_fit ∈ current_fit_idx:size(tmp_fit_data,1)
                    tmp_fit_data[current_fit] = fit_data
                    tmp_run_data[:,current_fit] = cur_best_pop
                    tmp_generations[current_fit] = params.number_of_generations
                end
            end
            # lock to prevent race conditions
            lock(thread_lock) do 
                
                thisrun = 1 + finished_runs % params.runs_per_bin
                curbin = (1 + finished_runs ÷ params.runs_per_bin)
                
                finished_runs += 1
                
                for fit_idx ∈ 1:size(fitness,1)
                    run_data[:,curbin,fit_idx] += tmp_run_data[:,fit_idx] ./ (tmp_fit_data[fit_idx])
                    weight_data[curbin,fit_idx] += 1.0 / (tmp_fit_data[fit_idx])
                    generations[curbin,fit_idx] += tmp_generations[fit_idx] 
                end
                
                # setting seed to 0 allows culling of used seeds in the checkpoint
                # ensuring each seed is used once
                seed_vec[1+thd-start_thread] = 0

                if verbose
                    println(@sprintf("  Bin %3u | Run %4u | Fitness %8.6f | Generations %u",curbin,thisrun,tmp_fit_data[end],numgen))
                end
                
                # Calculate bin data if enough to finish a bin, then checkpoint
                if finished_runs % params.runs_per_bin == 0
                    
                    bin_results!(bin_data,calculated_zeroth_moment,run_data,weight_data,curbin,Δω,Greens_tuple,fitness,seed_vec,generations,params)

                end # bin completing
            end # lock(thread_lock)
            

        end # threads
    end
    # time stats
    
    t_per_run = Δt/(total_runs - (start_bin-1)*params.runs_per_bin)
    t_per_bin = Δt/(params.num_bins + 1 - start_bin)

    # Merge data
    zero_tmp = mean(calculated_zeroth_moment,dims=2)[:,1,:]
    gen_per_run = sum(generations ./ params.runs_per_bin)[1]/params.num_bins
    zero_avg, zero_err = jackknife(zero_tmp)
    differential = 100.0*abs(1.0-zero_avg[1])
    data = zeros(Float64,(size(params.out_ωs,1),size(fitness,1)))
    err = zeros(Float64,(size(params.out_ωs,1),size(fitness,1)))
    
    for fit_idx ∈ 1:size(fitness,1)
        data[:,fit_idx], err[:,fit_idx] = jackknife(bin_data[:,:,fit_idx])
    end
    
    # Print statistics
    println("\nSaving data to ",params.output_file," and deleting checkpoint file\n")
    
    println("Run Statistics")
    if occursin("fermionic",params.kernel_type)
        println(@sprintf(" Expected 0th moment:   1.00") )
        println(@sprintf(" DEAC 0th moment:       %01.3f",zero_avg[1]))
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
            "fitness" =>  fitness,
            "runtime" => Δt,
            "full_eigenvalues" => full_eigen
        )
    else
        bin_dict = Dict{String,Any}(
            "A" => data,
            "σ" => err,
            "zeroth_moment" => zero_avg[1],
            "zeroth_moment_σ" => zero_err[1],
            "avg_generations" => gen_per_run,
            "ωs" => params.out_ωs,
            "fitness" =>  fitness,
            "runtime" => Δt,
            "full_eigenvalues" => full_eigen
        )
    end
    FileIO.save(params.output_file,bin_dict)
    delete_checkpoint(params)
    return bin_dict
end # run_DEAC()



