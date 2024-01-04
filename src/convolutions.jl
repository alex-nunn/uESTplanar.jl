"""
FFT convolution methods for computing the electrostatic potentials
"""

using FFTW
import Base: size

export KernelConvolution, round_pow2

"""
    KernelConvolution{T, M}

Convolution kernel type contains a grid information and the kernel convolution
matrix. Instances of KernelConvolution can be used to convolve 2-dimensional
functions.

The convolution is zero-padded and the type also allocates a buffer allowing the
kernel transform to be applied many times with minimal memory usage.

# Types
- `T<:AbstactFloat` floating point precision type
- `M::Int` dimension of kernel function. For scalar kernals `M =1`, for vector 
    kernels `M = 3`

# Properties
- `xs::StepRangeLen{T}` grid of `x` values in cartesian grid
- `ys::StepRangeLen{T}` grid of `y` values in cartesian grid
- `kermat::Array{Complex{T}}` M-dimensional Fourier space kernel matrix
- `buffer::Array{Complex{T}}` M-dimensional buffer

"""
struct KernelConvolution{T<:AbstractFloat, M}
    xs::StepRangeLen{T}
    ys::StepRangeLen{T}
    kermat::Array{Complex{T}}
    buffer::Array{Complex{T}}
end


"""
    KernelConvolution{T}(xs, ys, kernel, M=1)

Construct a kernel convolution over grid `xs`, `ys` using `kernel(x, y)` function.
The `M` specifies the dimension of the kernel function output.

For example, if the kernel is scalar function then `kernel(x, y)` should return a
scalar value and `M = 1`. If the kernel is a vector function where `kernel(x,y)`
outputs a vector of length `n` then `M = n`.

# Examples
Defining a 2-dimensional kernel function
```julia-repl
julia> kernel = (x, y) -> [cos(x), sin(y)];
julia> xs = range(0, 1, 128)
julia> KernelConvolution{Float64}(xs, xs, kernel, 2);
```
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

"""
    KernelConvolution(xs, ys, kernel, M=1)

Behaves identically to `KernelConvolution{T}` except the numerical type is
inferred from the arguments.
"""
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