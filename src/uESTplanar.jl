module uESTplanar

using Scratch
using Logging
cache = ""

include("kernels.jl")
include("convolutions.jl")


function __init__()
    global cache = @get_scratch!("cache")

    # Load fast kernel versions into runtime
    cache_path = joinpath(cache, "fast_verions.jldat")
    if !ispath(cache_path)
        @info "Precomputing fast kernels..."
        compute_fast_versions!(cache_path)
        @info "Fast kernels computed."
    end
    load_fast_verions!(cache_path)
end
end
