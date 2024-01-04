"""
FFT convolution methods for computing the electrostatic potentials
"""

using FFTW
import Base: size

export KernelConvolution, round_pow2

"""
    KernelConvolution{T, M}

Convolution kernel computed over a fixed 2-dimensional grid `xs`, `ys`. Can be
used to convolve 2-dimensional functions.

The convolution is zero-padded.

# Properties
- `xs::StepRangeLen` grid of `x` values in cartesian grid
- `ys::StepRangeLen` grid of `y` values in cartesian grid
- `kermatt` 

"""
struct KernelConvolution{T<:AbstractFloat, M}
    xs::StepRangeLen{T}
    ys::StepRangeLen{T}
    kermat::Array{Complex{T}}
    buffer::Array{Complex{T}}
end

# TODO: Docstring!!
"""

"""
function KernelConvolution{T}(
            xs::StepRangeLen{T}, ys::StepRangeLen{T}, kernel, M=1) where {T}
    # Build kernel matrix
    kxs = kernel_samples(xs)
    kys = kernel_samples(ys)

    kermat = zeros(Complex{T}, 2 * length(xs), 2 * length(ys), M)
    for i ∈ eachindex(kxs), j ∈ eachindex(kys)
        kermat[i, j, :] .= kernel(kxs[i], kys[j]) * step(xs) * step(ys)
    end
    fft!(kermat, 1:2)

    # Allocate buffer
    buffer = Array{Complex{T}}(undef, 2 * length(xs), 2 * length(ys), M)

    return KernelConvolution{T, M}(
        xs,
        ys,
        kermat,
        buffer
    )
end

function KernelConvolution(xs::StepRangeLen{T}, ys::StepRangeLen{T}, kernel, M=1) where {T}
    return KernelConvolution{T}(xs, ys, kernel, M)
end

Base.size(kconv::KernelConvolution{T, M}) where {T, M} = (length(kconv.xs), length(kconv.ys))


function (kconv::KernelConvolution{T, M})(out::Array{T, 3}, mat::Matrix{T}) where {T, M}
    if size(kconv) != size(mat)
        throw(DomainError("`mat` must have same size as KernelConvolution"))
    end

    # Prepare buffer
    kconv.buffer .= 0.0

    for i ∈ axes(mat, 1), j ∈ axes(mat, 2)
        kconv.buffer[i, j, 1] = mat[i, j]
    end

    fft!(view(kconv.buffer, :, :, 1), 1:2)
    @. kconv.buffer = kconv.kermat * kconv.buffer[:, :, 1]
    ifft!(kconv.buffer, 1:2)

    Nx, Ny = size(kconv)
    i0 = Nx ÷ 2 - 1
    j0 = Ny ÷ 2 - 1

    for i ∈ axes(out, 1), j ∈ axes(out, 2)
        @. out[i, j, :] = real(kconv.buffer[i0 + i, j0 + j, :])
    end
end

function (kconv::KernelConvolution{T, M})(mat::Matrix{T}) where {T, M}
    out = Array{T, 3}(undef, size(kconv)..., M)
    kconv(out, mat)
    return out
end

function (kconv::KernelConvolution{T, 1})(mat::Matrix{T}) where {T}
    out = Array{T, 2}(undef, size(kconv)...)
    kconv(reshape(out, size(kconv)..., 1), mat)
    return out
end

function kernel_samples(xs)
    N = length(xs)
    m = N ÷ 2
    return step(xs) * (iseven(N) ?  range(-m, m-1) : range(-m, m)) 
end

"""
    round_pow2(n::Integer)

Round an integer to the nearest power of two.

# Examples
```julia-repl
julia> round_pow2(4)
4
julia> round_pow2(5)
4
```
"""
function round_pow2(n::Integer)
    return 1 << round(Int, log2(n))
end

"""
    round_pow2(xs::AbstractRange)

Modify the step size so the range length is the nearest power of two.

# Examples
```julia-repl
julia> xs = 1:7;
julia> round_pow2(xs)
1:1:8
```
"""
function round_pow2(xs::AbstractRange)
    n = round_pow2(length(xs))
    return range(
        start=xs[begin],
        stop=xs[end],
        length=n
    )
end