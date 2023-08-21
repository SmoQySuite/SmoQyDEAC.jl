
# check if checkpoint file exists, if so return dictionary saved in
# checkpoint file
function find_checkpoint(params::DEACParameters)
    file = params.checkpoint_directory*"/DEAC_checkpoint.jld2"
    check_exists = isfile(file)
    if check_exists
        check_dict = FileIO.load(file)
        return true, check_dict
    else
        return false, nothing
    end
end # find_checkpoint()


# Check if parameters that must be identical between current run and checkpoint 
# are identical
function compare_checkpoint(checkpoint_dict,params::DEACParameters,G_tuple)
    check_params = checkpoint_dict["params"]
    cor_dat = checkpoint_dict["G_tuple"]
    return check_params == params && cor_dat == G_tuple
end # compare_checkpoint()

# Deletes checkpoint file
function delete_checkpoint(params::DEACParameters)
    file = params.checkpoint_directory*"/DEAC_checkpoint.jld2"
    if isfile(file)
        rm(file)
    end
end # delete_checkpoint()

# Save a checkpoint file
function save_checkpoint(bin_data, generations, bin_num, params::DEACParameters,G_tuple,zeroth_momentum::AbstractArray,true_fitness,seed_vec)
    file = params.checkpoint_directory*"/DEAC_checkpoint.jld2"
    seeds = filter(x->x≠0,seed_vec)
    chk_data = Dict{String,Any}(
        "bin_data" => bin_data,
        "generations" => generations,
        "bin_num" => bin_num,
        "params" => params,
        "G_tuple" => G_tuple,
        "zeroth" => zeroth_momentum,
        "true_fitness" => true_fitness,
        "seeds" => seeds
    )
    FileIO.save(file,chk_data)
end # save_checkpoint()
