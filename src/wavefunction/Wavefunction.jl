"""
    Wavefunction

Module for quantum wavefunction representation in spherical coordinates.

Represents wavefunctions in the radial-angular basis with expansion:
    ψ(r,θ,φ) = Σₗ Σₘ g(r,m,l) Yₗᵐ(θ,φ)

where g(r,m,l) are radial coefficients on the GPS grid.

# Implementation Notes

## Fortran Blueprint Fidelity
Matches the array structure from D_inner_out_volkov_3d_with_prob.f90:
- Array g(nrmax, -lmax:lmax, 0:lmax) stores radial coefficients
- Complex*16 representation for time evolution
- Spherical harmonic quantum numbers: l ∈ [0,lmax], m ∈ [-l,l]

## Memory Layout
Julia uses column-major order (like Fortran), so direct translation is efficient.
"""
module Wavefunction

using Printf
using LinearAlgebra

export WavefunctionData, create_wavefunction, initialize_ground_state!
export compute_norm, compute_radial_density, copy_wavefunction
export get_component, set_component!, print_wavefunction_info

"""
    WavefunctionData

Container for quantum wavefunction in radial-angular basis.

# Fields
- `g::Array{ComplexF64,3}`: Radial coefficients [nrmax, 2lmax+1, lmax+1]
  - First index: radial grid point
  - Second index: m quantum number (shifted to 1-based: m+lmax+1)
  - Third index: l quantum number (shifted to 1-based: l+1)
- `nrmax::Int`: Number of radial grid points
- `lmax::Int`: Maximum angular momentum
- `grid_weights::Vector{Float64}`: GPS quadrature weights for integration

# Indexing Convention
To access g(r,m,l) from physics notation:
    data.g[ir, m+lmax+1, l+1]
where ir ∈ [1,nrmax], m ∈ [-l,l], l ∈ [0,lmax]
"""
struct WavefunctionData
    g::Array{ComplexF64,3}  # Radial coefficients: g[ir, m_index, l_index]
    nrmax::Int
    lmax::Int
    grid_weights::Vector{Float64}  # For norm computation
end

"""
    create_wavefunction(nrmax::Int, lmax::Int, grid_weights::Vector{Float64}) -> WavefunctionData

Create empty wavefunction container with proper array dimensions.

# Arguments
- `nrmax::Int`: Number of radial grid points
- `lmax::Int`: Maximum angular momentum quantum number
- `grid_weights::Vector{Float64}`: GPS quadrature weights from grid

# Returns
- `WavefunctionData`: Initialized with zeros

# Example
```julia
using GPSGrid
grid = create_gps_grid(400, 150.0)
wfn = create_wavefunction(grid.nrmax, 50, grid.quadrature_weights)
```
"""
function create_wavefunction(nrmax::Int, lmax::Int, grid_weights::Vector{Float64})
    if nrmax < 10
        error("nrmax too small: $nrmax (need at least 10)")
    end
    if lmax < 0
        error("lmax must be non-negative: $lmax")
    end
    if length(grid_weights) != nrmax
        error("grid_weights length $(length(grid_weights)) does not match nrmax $nrmax")
    end

    # Allocate array with proper dimensions
    # Dimension 2: m quantum number ranges from -lmax to +lmax (2*lmax+1 values)
    # Dimension 3: l quantum number ranges from 0 to lmax (lmax+1 values)
    g = zeros(ComplexF64, nrmax, 2*lmax+1, lmax+1)

    return WavefunctionData(g, nrmax, lmax, grid_weights)
end

"""
    initialize_ground_state!(wfn::WavefunctionData,
                            ground_state_radial::Vector{Float64},
                            n_ground::Int=1, l_ground::Int=0)

Initialize wavefunction with ground state eigenvector.

# Arguments
- `wfn::WavefunctionData`: Wavefunction to initialize (modified in-place)
- `ground_state_radial::Vector{Float64}`: Ground state radial wavefunction φ(r)
- `n_ground::Int=1`: Principal quantum number (for documentation)
- `l_ground::Int=0`: Angular momentum quantum number

# Physics
For hydrogen-like atoms, ground state is typically (n=1, l=0, m=0):
    ψ₁ₛ(r,θ,φ) = R₁₀(r) Y₀⁰(θ,φ)

Only the (l=0, m=0) component is populated with the radial function.

# Example
```julia
# After solving eigenvalue problem for Hamiltonian
ground_state = eigenvectors[:, 1]  # Lowest energy eigenstate
initialize_ground_state!(wfn, ground_state, 1, 0)
```

# Implementation Notes
Follows Fortran code initialization (Lines 400-420 in reference program):
- Sets g(ir, m=0, l=0) = ground_state_radial[ir]
- All other components remain zero
- Normalizes the wavefunction
"""
function initialize_ground_state!(wfn::WavefunctionData,
                                  ground_state_radial::Vector{Float64},
                                  n_ground::Int=1, l_ground::Int=0)
    if length(ground_state_radial) != wfn.nrmax
        error("Ground state radial function length $(length(ground_state_radial)) " *
              "does not match wfn.nrmax $(wfn.nrmax)")
    end

    if l_ground > wfn.lmax
        error("l_ground ($l_ground) exceeds wfn.lmax ($(wfn.lmax))")
    end

    # Zero out all components
    fill!(wfn.g, 0.0 + 0.0im)

    # Set (l=l_ground, m=0) component
    # Indexing: g[ir, m+lmax+1, l+1]
    # For m=0, l=l_ground: g[ir, lmax+1, l_ground+1]
    m_index = wfn.lmax + 1  # m=0
    l_index = l_ground + 1

    for ir in 1:wfn.nrmax
        wfn.g[ir, m_index, l_index] = ground_state_radial[ir] + 0.0im
    end

    # Normalize
    norm = compute_norm(wfn)
    wfn.g ./= norm

    @info "Initialized ground state" n=n_ground l=l_ground norm_before=norm norm_after=compute_norm(wfn)
end

"""
    compute_norm(wfn::WavefunctionData) -> Float64

Compute L² norm of wavefunction using GPS quadrature.

# Formula
    ⟨ψ|ψ⟩ = Σᵣ Σₗ Σₘ |g(r,m,l)|² wᵣ

where wᵣ are the GPS quadrature weights (includes r² and dr/dx Jacobian).

# Returns
- `Float64`: Norm ⟨ψ|ψ⟩ (should be ≈ 1.0 for normalized states)

# Physics
GPS weights already include:
- Gauss-Legendre weights
- Jacobian dr/dx from algebraic mapping
- Properly accounts for volume element in spherical coordinates
"""
function compute_norm(wfn::WavefunctionData)
    norm_sq = 0.0

    # Sum over all radial points, m values, and l values
    for l_index in 1:(wfn.lmax+1)
        l = l_index - 1
        for m_index in 1:(2*wfn.lmax+1)
            m = m_index - (wfn.lmax + 1)
            # Check if m is valid for this l
            if abs(m) <= l
                for ir in 1:wfn.nrmax
                    norm_sq += abs2(wfn.g[ir, m_index, l_index]) * wfn.grid_weights[ir]
                end
            end
        end
    end

    return sqrt(norm_sq)
end

"""
    compute_radial_density(wfn::WavefunctionData) -> Vector{Float64}

Compute radial probability density ρ(r) = Σₗ Σₘ |g(r,m,l)|².

# Returns
- `Vector{Float64}`: Radial density at each grid point [nrmax]

# Physics
Integrates over angular coordinates:
    ρ(r) = ∫ |ψ(r,θ,φ)|² sin(θ) dθ dφ = Σₗ Σₘ |g(r,m,l)|²

Useful for visualization and ionization analysis.
"""
function compute_radial_density(wfn::WavefunctionData)
    density = zeros(Float64, wfn.nrmax)

    for l_index in 1:(wfn.lmax+1)
        l = l_index - 1
        for m_index in 1:(2*wfn.lmax+1)
            m = m_index - (wfn.lmax + 1)
            if abs(m) <= l
                for ir in 1:wfn.nrmax
                    density[ir] += abs2(wfn.g[ir, m_index, l_index])
                end
            end
        end
    end

    return density
end

"""
    copy_wavefunction(wfn::WavefunctionData) -> WavefunctionData

Create a deep copy of wavefunction data.

# Returns
- `WavefunctionData`: Independent copy with new array allocation

# Use Cases
- Checkpointing during time evolution
- Region splitting (inner/outer regions)
- Backup before applying absorbers
"""
function copy_wavefunction(wfn::WavefunctionData)
    g_copy = copy(wfn.g)
    weights_copy = copy(wfn.grid_weights)
    return WavefunctionData(g_copy, wfn.nrmax, wfn.lmax, weights_copy)
end

"""
    get_component(wfn::WavefunctionData, l::Int, m::Int) -> Vector{ComplexF64}

Extract radial function for specific (l,m) component.

# Arguments
- `wfn::WavefunctionData`: Wavefunction
- `l::Int`: Angular momentum quantum number
- `m::Int`: Magnetic quantum number

# Returns
- `Vector{ComplexF64}`: Radial function g(r,m,l) for all r

# Example
```julia
# Get ground state radial function (l=0, m=0)
ground_radial = get_component(wfn, 0, 0)
```
"""
function get_component(wfn::WavefunctionData, l::Int, m::Int)
    if l > wfn.lmax || l < 0
        error("l=$l out of range [0, $(wfn.lmax)]")
    end
    if abs(m) > l
        error("m=$m out of range [-$l, $l] for l=$l")
    end

    l_index = l + 1
    m_index = m + wfn.lmax + 1

    return wfn.g[:, m_index, l_index]
end

"""
    set_component!(wfn::WavefunctionData, l::Int, m::Int, radial::Vector{ComplexF64})

Set radial function for specific (l,m) component.

# Arguments
- `wfn::WavefunctionData`: Wavefunction (modified in-place)
- `l::Int`: Angular momentum quantum number
- `m::Int`: Magnetic quantum number
- `radial::Vector{ComplexF64}`: Radial function values [nrmax]
"""
function set_component!(wfn::WavefunctionData, l::Int, m::Int, radial::Vector{ComplexF64})
    if l > wfn.lmax || l < 0
        error("l=$l out of range [0, $(wfn.lmax)]")
    end
    if abs(m) > l
        error("m=$m out of range [-$l, $l] for l=$l")
    end
    if length(radial) != wfn.nrmax
        error("radial length $(length(radial)) does not match nrmax $(wfn.nrmax)")
    end

    l_index = l + 1
    m_index = m + wfn.lmax + 1

    wfn.g[:, m_index, l_index] .= radial
end

"""
    print_wavefunction_info(wfn::WavefunctionData)

Print diagnostic information about wavefunction.
"""
function print_wavefunction_info(wfn::WavefunctionData)
    println("Wavefunction Information:")
    println("  Grid points (nrmax): $(wfn.nrmax)")
    println("  Max angular momentum (lmax): $(wfn.lmax)")
    println("  Array dimensions: $(size(wfn.g))")
    println("  Memory usage: $(@sprintf("%.2f", sizeof(wfn.g) / 1024^2)) MB")
    println("  Norm: $(@sprintf("%.6f", compute_norm(wfn)))")

    # Find dominant components
    max_amplitude = 0.0
    max_l = 0
    max_m = 0
    for l_index in 1:(wfn.lmax+1)
        l = l_index - 1
        for m_index in 1:(2*wfn.lmax+1)
            m = m_index - (wfn.lmax + 1)
            if abs(m) <= l
                for ir in 1:wfn.nrmax
                    amp = abs(wfn.g[ir, m_index, l_index])
                    if amp > max_amplitude
                        max_amplitude = amp
                        max_l = l
                        max_m = m
                    end
                end
            end
        end
    end
    println("  Dominant component: (l=$max_l, m=$max_m) with max amplitude $(@sprintf("%.6f", max_amplitude))")
end

end  # module Wavefunction
