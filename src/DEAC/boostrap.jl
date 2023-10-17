
# Create bootstrapped bins from current number of bins. 
function bootstrap_samples(corr_bin_in,N_boots,seed)
    NBins_in = size(corr_bin_in,1)
    nτ = size(corr_bin_in,2)

    corr_bin_out = Array{Float64}(undef,(N_boots,nτ))

    rng = Random.Xoshiro(seed)

    for boot in 1:N_boots
        indices = (rand(rng,UInt32,NBins_in) .% NBins_in) .+ 1
        sum = zeros(Float64,nτ)
        for i in 1:NBins_in
            sum .+= corr_bin_in[indices[i],:]
        end
        corr_bin_out[boot,:] = sum ./ NBins_in 
    end # bootstrap bins
    return corr_bin_out
end # bootstrap_samples()
