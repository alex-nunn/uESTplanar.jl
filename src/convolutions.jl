"""
FFT convolution methods for computing the electrostatic potentials
"""

using FFTW

export KernelConvolution


Base.@kwdef struct KernelConvolution{T}
    xs::StepRangeLen{T}
    ys::StepRangeLen{T}
    kermat::Matrix{Complex{T}}
    buffer::Matrix{Complex{T}}
end

size(kconv::KernelConvolution) = (length(xs), length(ys))

KernelConvolution{T}(xs::T, ys::T, kernel) where {T} = begin
    # Build kernel matrix
    kxs = kernel_samples(xs)
    kys = kernel_samples(ys)

    kermat = zeros(Complex{T}, 2 * length(xs), 2 * length(ys))
    for i ∈ eachindex(kxs), j ∈ eachindex(kys)
        kermat[i, j] = kernel(kxs[i], kys[j]) * step(xs) * step(ys)
    end
    ftt!(kermat)

    # Allocate buffer
    buffer = Matrix{Complex{T}}(undef, 2 * length(xs), 2 * length(ys))

    return KernelConvolution{T}(
        xs,
        ys,
        kermat,
        buffer
    )
end


function (kconv::KernelConvolution{T})(out::Matrix{T}, mat::Matrix{T}) where {T}
    if size(kconv) != size(mat)
        throw(DomainError("`mat` must have same size as KernelConvolution"))
    end

    # Prepare buffer
    kconv.buffer .= 0.0

    for i ∈ axes(mat, 1), j ∈ axes(mat, 2)
        kconv.buffer[i, j] = mat[i, j]
    end

    fft!(kconv.buffer)
    @. kconv.buffer *= kconv.kermat
    ifft!(kconv.buffer)

    Nx, Ny = size(kconv)
    i0 = Nx ÷ 2 - 1
    j0 = Ny ÷ 2 - 1

    for i ∈ 1:Nx, j ∈ 1:Ny
        out[i, j] = real(result[i0 + i, j0 + j])
    end
end


function (kconv::KernelConvolution{T})(mat::Matrix{T}) where {T}
    out = Matrix{T}(undef, size(kconv)...)
    kconv(out, mat)
    return out
end


function kernel_samples(xs)
    N = length(xs)
    m = N ÷ 2
    return step(xs) * (iseven(N) ?  range(-m, m-1) : range(-m, m)) 
end