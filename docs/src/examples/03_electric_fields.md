# Electric fields

We can also compute the electric field using convolutions. If the potential on
the bottom plate, ``z = 0``, is given by ``f(x,y)`` then the electric field has
the form,

```math
\begin{aligned}
E &= - \nabla \phi(x, y, z),\\
&= -\nabla\iint K(x - x', y - y', z) f(x', y') dx'dy'. 
\end{aligned}
```

The `uESTplanar` module provides the partial derivatives of the kernel function
``K(x, y, z)`` with respect to each cartesian coordinate as [`∂xK`](@ref), 
[`∂yK`](@ref), [`∂zK`](@ref). The example below shows how to use these functions
to compute the electric field.

```@example
using Plots
using LaTeXStrings
using uESTplanar

# 1. Define kernel for each electric field component
z = 0.25
kconv = KernelConvolution{Float32}(
    range(-5, 10; step=0.01f0),
    range(-5, 5; step=0.01f0),
    (x, y)-> [
        ∂xK(x, y, z) + ∂xK(x, y, 1 - z),
        ∂yK(x, y, z) + ∂yK(x, y, 1 - z),
        ∂zK(x, y, z) - ∂zK(x, y, 1 - z),
    ],
    3
)

# 2. Define a potential function with a wavy boundary
x0 = 0.0  # position of interface
λ = 2.0  # wavelength of variation
A = 2.0   # amplitude of variation
wavy_boundary = (x, y) -> (A * cos(2π / λ * y) < x - x0) ? 1.0f0 : 0.0f0

# 3. Compute potential over grid of convolution
boundary_matrix = wavy_boundary.(kconv.xs, kconv.ys')

# 4. Apply convolution operator
electric_field = kconv(boundary_matrix)

# 5. Plot potential
(; xs, ys) = kconv
pkws = (
    xlabel="x",
    ylabel="y",
    xlims=(-5, 10),
    ylims=(-5, 5),
    size=(600, 400),
    colorbar=false,
    margin=3mm,
)

p1 = heatmap(
    xs, ys, transpose(boundary_matrix);
    title="boundary potential",
    pkws...
)

ps = [
    heatmap(
        xs, ys, transpose(electric_field[:, :, i]);
        title=L"E_%$label",
        titlefontsize=16,
        pkws...
    )
    for (i, label) ∈ enumerate("xyz")
]


p = plot(p1, ps...; layout=(2, 2), size=(800, 800))
savefig("03_electric_example.png"); nothing # hide
```

![](03_electric_example.png)