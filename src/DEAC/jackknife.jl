
# Calculate average and error using jackknife method
function jackknife(datums::AbstractArray{T}) where {T}
    nbin = size(datums,2)
    avg = Statistics.mean(datums, dims=2)[:,1]
    
    datums = datums .- avg
    err = sum(datums .* datums,dims = 2) ./ (nbin * ( nbin -1 ))
    err = sqrt.(err)[:,1]
    return avg, err
end # jacknife