# 1. Convolutions in 2D dimensions

The [`KernelConvolution`](@ref) type allows us to rapidly convolve 2-dimensional functions with a kernel function. 

For example, we can convolve the 2-dimensional box function ``f(x, y)`` with a Gaussian kernel using the following code.

```@example
using Plots
using uESTplanar

# 1. Define convolution kernel
σ = 0.3
kconv = KernelConvolution(
    range(-2, 2; step=0.005),
    range(-2, 2; step=0.005),
    (x, y)-> exp(-(x^2 + y^2) / (2σ^2)) / (sqrt(2π) * σ)
)

# 2. Define box function & compute over grid
f = (x, y) ->  -1 < x < 1 && -1 < y < 1 ? 1.0 : 0.0
(; xs, ys) = kconv
boundary_values = f.(xs, ys')

# 3. Apply kernel transform
result = kconv(boundary_values)

# 4. Plot potential
heatmap(
    xs, ys, transpose(result);
    xlabel="x",
    ylabel="y",
    size=(450, 450),
    aspect_ratio=:equal,
    dpi=100
)
savefig("01_convolution_example.png"); nothing # hide
```

![](01_convolution_example.png)


The Fast Fourier Transform (FFT) used to perform the 2-dimensional convolution is faster on arrays where the length is a power of two. We can sometimes increase evaluation speed by rounding our grid size to the nearest power of two using the [`round_pow2`](@ref) method.

For example,
```@example
using uESTplanar

# Kernel
σ = 0.3
kernel = (x, y)-> exp(-(x^2 + y^2) / (2σ^2)) / (sqrt(2π) * σ)

# Grid
xs = range(-2.5, 2.5; step=0.005)
xs_pow2 = round_pow2(xs)

# Convolutions
kc1 = KernelConvolution(xs, xs, kernel)
kc2 = KernelConvolution(xs_pow2, xs_pow2, kernel)

# Random boundary potentials
fs1 = randn(size(kc1)...)
fs2 = randn(size(kc2)...)

# Compare evaluation speeds
t1 = @timed kc1(fs1)
t2 = @timed kc2(fs2)

println("normal grid : $(t1.time) sec")
println("power 2 grid : $(t2.time) sec")
```