
# Create bootstrapped bins from current number of bins. 
function bootstrap_samples(corr_bin_in,N_boots,seed)
    NBins_in = size(corr_bin_in,1)
    nτ = size(corr_bin_in,2)

    corr_bin_out = Array{eltype(corr_bin_in)}(undef,(N_boots,nτ))

    Threads.@threads for boot in 1:N_boots
        rng = Random.Xoshiro(seed + boot)
        indices = (rand(rng,UInt32,NBins_in) .% NBins_in) .+ 1
        sum = zeros(eltype(corr_bin_in),nτ)
        for i in 1:NBins_in
            sum .+= corr_bin_in[indices[i],:]
        end
        corr_bin_out[boot,:] = sum ./ NBins_in 
    end # bootstrap bins
    return corr_bin_out
end # bootstrap_samples()
