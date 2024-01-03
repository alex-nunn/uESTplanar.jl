"""
Kernels for computing the electrostatic field in a planar electrostatic trap via
convolution.
"""

using QuadGK
using SpecialFunctions
using Serialization
using Interpolations

export K, ∂xK, ∂yK, ∂zK

# Components for computing convolution kernel
I = (ρ, z) -> begin 
    result, error = quadgk(
        θ -> real(polygamma(1, 0.5 * (z - 1.0im * ρ * sin(θ)))), 
        0, π/2
    )
    return result
end

I_asym = (ρ, z) -> begin
    D = sqrt(ρ^2 + z^2)
    return 4real(atan((z -1.0im * ρ) / D)) / D
end

∂ρI = (ρ, z) -> begin 
    result, error = quadgk(
        θ -> 0.5sin(θ) * imag(polygamma(2, 0.5 * (z - 1.0im * ρ * sin(θ)))), 
        0, π/2
    )
    return result
end

∂ρI_asym = (ρ, z) -> begin 
    D = sqrt(ρ^2 + z^2)
    return -4ρ * real(atan((z-1.0im * ρ) / D)) / D^3
end

∂zI = (ρ, z) -> begin 
    result, error = quadgk(
        θ -> 0.5real(polygamma(2, 0.5 * (z - 1.0im * ρ * sin(θ)))), 
        0, π/2
    )
    return result
end

∂zI_asym = (ρ, z) -> begin 
    D = sqrt(ρ^2 + z^2)
    return 4z * real(atan((z-1.0im * ρ) / D)) / D^3
end


# Convolution kernels for evaluating the potential and associated derivatives
"""
    K(x, y, z)

Planar electrostatic trap kernel in cartesian coordinates where `0 < z < 1` and
`x`, `y` are real valued horizontal coordinates. The bottom and top plates of
the trap are located at `z = 0` and `z = 1`, respectively.

Useful for computing the electrostatic field via convolution.

See also [`∂xK`](@ref), [`∂yK`](@ref), [`∂zK`](@ref).
"""
function K(x, y, z)
    ρ = sqrt(x^2 + y^2)
    return (I_fast(ρ, z) - I_fast(ρ, 2 - z)) / 4π^2
end

"""
    K(ρ, z)

Planar electrostatic trap kernel in cylindrical coordinates where `0 < ρ`,
`0 < z < 1`.
"""
function K(ρ, z)
    return (I_fast(ρ, z) - I_fast(ρ, 2 - z)) / 4π^2
end

"""
    ∂xK(x, y, z)

Partial derivative wrt. `x` of planar electrostatic trap kernel in cartersian 
coordinates where `0 < z < 1`.

See also [`K`](@ref)
"""
∂xK = (x, y, z) -> begin
    ρ = sqrt(x^2 + y^2)
    if ρ < 1e-12 return 0.0 end
    return (x / ρ) * (∂ρI_fast(ρ, z) - ∂ρI_fast(ρ, 2 - z)) / 4π^2
end

"""
    ∂yK(x, y, z)

Partial derivative wrt. `y` of planar electrostatic trap kernel in cartersian 
coordinates where `0 < z < 1`.

See also [`K`](@ref)
"""
∂yK = (x, y, z) -> begin
    ρ = sqrt(x^2 + y^2)
    if ρ < 1e-12 return 0.0 end
    return (y / ρ) * (∂ρI_fast(ρ, z) - ∂ρI_fast(ρ, 2 - z)) / 4π^2
end

"""
    ∂zK(x, y, z)

Partial derivative wrt. `z` of planar electrostatic trap kernel in cartersian 
coordinates where `0 < z < 1`.

See also [`K`](@ref)
"""
∂zK = (x, y, z) -> begin
    ρ = sqrt(x^2 + y^2)
    return (∂zI_fast(ρ, z) + ∂zI_fast(ρ, 2 - z)) / 4π^2
end


# Use combination of interpolation and asymptotic forms to increase evaluation
# performance
function fast_factory(F, F_asmp, ρ_max=8)
    ρs = range(0, ρ_max, 2000)
    zs = range(0, 2, 300)[2:end]
    values = F.(ρs, zs')
    F_interp = cubic_spline_interpolation((ρs, zs), values)
    #F_interp = scale(
    #    interpolate(values, BSpline(Cubic(Flat(OnGrid())))),
    #    ρs, zs
    #)

    F_fast = (ρ, z) -> begin
        if ρ > ρ_max
            return F_asmp(ρ, z)
        else
            return F_interp(ρ, z)
        end
    end
    return F_fast
end

function compute_fast_versions!(cache_path)
    I_fast = fast_factory(I, I_asym)
    ∂ρI_fast = fast_factory(∂ρI, ∂ρI_asym)
    ∂zI_fast = fast_factory(∂zI, ∂zI_asym)

    serialize(cache_path, (;I_fast, ∂ρI_fast, ∂zI_fast))
end

function load_fast_verions!(cache_path)
    global (;I_fast, ∂ρI_fast, ∂zI_fast) = deserialize(cache_path)
end