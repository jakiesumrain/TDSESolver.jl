"""
    AngularGrid

Module for creating angular grids (θ, φ) for 3D TDSE calculations.

Provides uniform and Gauss-Legendre grids for angular coordinates
in spherical polar representation of wavefunctions.
"""
module AngularGrid

using FastGaussQuadrature
using GSL: sf_legendre_sphPlm
using Printf  # For @sprintf

export AngularGridData, create_angular_grid, compute_solid_angle_element
export compute_direction_vectors

"""
    AngularGridData

Angular grid specification for 3D calculations with precomputed transformation arrays.

# Fields
## Basic Grid
- `nthmax::Int`: Number of θ points (corresponds to nthmax_sp in Fortran)
- `nphimax::Int`: Number of φ points (corresponds to nfimax_sp in Fortran)
- `theta_grid::Vector{Float64}`: θ values in [0, π]
- `phi_grid::Vector{Float64}`: φ values in [0, 2π]
- `theta_weights::Vector{Float64}`: Quadrature weights for θ
- `phi_weights::Vector{Float64}`: Quadrature weights for φ

## Transformation Arrays (matching Fortran)
- `plgd_sp::Array{Float64,3}`: Associated Legendre polynomials [l, m, nth]
  Corresponds to plgd_sp_lm(l,m,nth) in Fortran (lines 644, 702)
- `expfi_sp::Matrix{ComplexF64}`: exp(imφ) phase factors [m, nfi]
  Corresponds to expfi_sp(m,nfi) in Fortran (line 650)
- `expfi_sp_conj::Matrix{ComplexF64}`: exp(-imφ) conjugate [nfi, m]
  Corresponds to expfi_sp_conj(nfi,m) in Fortran (line 694)
- `direction::Array{Float64,3}`: Cartesian unit vectors [component, nfi, nth]
  Corresponds to direction(1:3,nfi,nth) in Fortran (line 663)

# Fortran Reference
Lines 638-708: Transformation between radial and angular representations
"""
struct AngularGridData
    # Basic grid
    nthmax::Int
    nphimax::Int
    theta_grid::Vector{Float64}
    phi_grid::Vector{Float64}
    theta_weights::Vector{Float64}
    phi_weights::Vector{Float64}

    # Transformation arrays (matching Fortran)
    plgd_sp::Array{Float64,3}           # [l, m, nth] - Associated Legendre
    expfi_sp::Matrix{ComplexF64}        # [m, nfi] - exp(imφ)
    expfi_sp_conj::Matrix{ComplexF64}   # [nfi, m] - exp(-imφ)
    direction::Array{Float64,3}          # [1:3, nfi, nth] - Cartesian vectors
end

"""
    compute_associated_legendre(l::Int, m::Int, x::Float64) -> Float64

Compute **spherical harmonic normalized** associated Legendre polynomial at x.

# Arguments
- `l::Int`: Degree (l ≥ 0)
- `m::Int`: Order (|m| ≤ l)
- `x::Float64`: Argument in [-1, 1]

# Returns
- Spherical harmonic normalized P_l^m(x) with normalization:
  sqrt((2l+1)/(4π) * (l-|m|)!/(l+|m|)!) * P_l^|m|(x)

# Note
Uses GSL sf_legendre_sphPlm which gives the correct normalization for
spherical harmonics Y_l^m = P_l^m(cos θ) * exp(imφ) (without additional factors).

This ensures that coordinate transformations preserve wavefunction norm.
"""
function compute_associated_legendre(l::Int, m::Int, x::Float64)
    m_abs = abs(m)
    if m_abs > l
        return 0.0
    end

    # GSL function: sf_legendre_sphPlm(l, m_abs, x)
    # Returns spherical harmonic normalized P_l^|m|(x)
    # Normalization: sqrt((2l+1)/(4π) * (l-m)!/(l+m)!)
    return sf_legendre_sphPlm(l, m_abs, x)
end

"""
    create_angular_grid(nthmax::Int, nphimax::Int, lmax::Int; method::Symbol=:gauss) -> AngularGridData

Create angular grid for spherical coordinates with precomputed transformation arrays.

# Arguments
- `nthmax::Int`: Number of θ (polar angle) points
- `nphimax::Int`: Number of φ (azimuthal angle) points
- `lmax::Int`: Maximum angular momentum (for precomputing Legendre polynomials)
- `method::Symbol=:gauss`: Grid type (:uniform or :gauss)

# Returns
- `AngularGridData`: Complete angular grid specification with transformation arrays

# Grid types:
- `:uniform`: Uniformly spaced points (simple but less efficient)
- `:gauss`: Gauss-Legendre points for θ, uniform for φ (optimal for integrals)

# Precomputed Arrays
Following Fortran implementation (lines 638-708):
- Associated Legendre polynomials P_l^|m|(cos θ) at all grid points
- Phase factors exp(imφ) and exp(-imφ)
- Cartesian direction vectors (sin θ cos φ, sin θ sin φ, cos θ)

# Example
```julia
grid = create_angular_grid(180, 60, 50, method=:gauss)
println("Angular grid: \$(grid.nthmax) × \$(grid.nphimax) points")
```

# Fortran Reference
Lines 295-310: Angular grid initialization
Lines 638-708: Usage in transformations
"""
function create_angular_grid(nthmax::Int, nphimax::Int, lmax::Int; method::Symbol=:gauss)
    # Validate inputs
    if nthmax < 2
        error("nthmax must be at least 2, got $nthmax")
    end
    if nphimax < 2
        error("nphimax must be at least 2, got $nphimax")
    end
    if lmax < 0
        error("lmax must be non-negative, got $lmax")
    end

    @info "Creating angular grid" nthmax=nthmax nphimax=nphimax lmax=lmax method=method

    # Create basic grid points and weights
    if method == :uniform
        # Uniform spacing
        theta_grid = range(0.0, π, length=nthmax) |> collect
        phi_grid = range(0.0, 2π, length=nphimax+1)[1:end-1] |> collect

        # Trapezoidal rule weights
        Δθ = π / (nthmax - 1)
        theta_weights = fill(Δθ, nthmax)
        theta_weights[1] /= 2.0
        theta_weights[end] /= 2.0

        Δφ = 2π / nphimax
        phi_weights = fill(Δφ, nphimax)

    elseif method == :gauss
        # Gauss-Legendre for θ using x = cos(θ) transformation
        # This is the correct approach for spherical harmonics with P_l^m(cos θ)
        #
        # For ∫₀^π f(θ) sin(θ) dθ:
        # Substitute x = cos(θ), then sin(θ) dθ = -dx
        # ∫₀^π f(θ) sin(θ) dθ = ∫₁^{-1} f(arccos(x)) (-dx) = ∫₋₁¹ f(arccos(x)) dx
        # Use Gauss-Legendre quadrature: ∫₋₁¹ g(x) dx ≈ Σ g(x_i) w_i

        x_gl, w_gl = gausslegendre(nthmax)

        # θ = arccos(x), NOT linear mapping!
        theta_grid = @. acos(x_gl)

        # Weights are w_gl directly (sin(θ) factor already in transformation)
        theta_weights = w_gl

        # Uniform for φ
        phi_grid = range(0.0, 2π, length=nphimax+1)[1:end-1] |> collect
        Δφ = 2π / nphimax
        phi_weights = fill(Δφ, nphimax)

    else
        error("Unknown grid method: $method. Use :uniform or :gauss")
    end

    # Precompute transformation arrays matching Fortran
    @info "Precomputing transformation arrays..."

    # 1. Associated Legendre polynomials plgd_sp[l, m, nth]
    # Fortran: plgd_sp_lm(l,m,nth) used in lines 644, 702
    plgd_sp = zeros(Float64, lmax+1, 2*lmax+1, nthmax)
    for nth in 1:nthmax
        cos_theta = cos(theta_grid[nth])
        for l in 0:lmax
            for m in -l:l
                plgd_sp[l+1, m+lmax+1, nth] = compute_associated_legendre(l, m, cos_theta)
            end
        end
    end

    # 2. Phase factors exp(imφ) - expfi_sp[m, nfi]
    # Fortran: expfi_sp(m,nfi) used in line 650
    expfi_sp = zeros(ComplexF64, 2*lmax+1, nphimax)
    for nfi in 1:nphimax
        φ = phi_grid[nfi]
        for m in -lmax:lmax
            expfi_sp[m+lmax+1, nfi] = exp(im * m * φ)
        end
    end

    # 3. Conjugate phase factors exp(-imφ) - expfi_sp_conj[nfi, m]
    # Fortran: expfi_sp_conj(nfi,m) used in line 694
    expfi_sp_conj = zeros(ComplexF64, nphimax, 2*lmax+1)
    for nfi in 1:nphimax
        φ = phi_grid[nfi]
        for m in -lmax:lmax
            expfi_sp_conj[nfi, m+lmax+1] = exp(-im * m * φ)
        end
    end

    # 4. Cartesian direction vectors - direction[1:3, nfi, nth]
    # Fortran: direction(1:3,nfi,nth) used in line 663
    # direction(1) = sin(θ)cos(φ) - x component
    # direction(2) = sin(θ)sin(φ) - y component
    # direction(3) = cos(θ)       - z component
    direction = zeros(Float64, 3, nphimax, nthmax)
    for nth in 1:nthmax
        θ = theta_grid[nth]
        sin_theta = sin(θ)
        cos_theta = cos(θ)

        for nfi in 1:nphimax
            φ = phi_grid[nfi]
            direction[1, nfi, nth] = sin_theta * cos(φ)  # x
            direction[2, nfi, nth] = sin_theta * sin(φ)  # y
            direction[3, nfi, nth] = cos_theta           # z
        end
    end

    @info "  Precomputed plgd_sp: $(size(plgd_sp))"
    @info "  Precomputed expfi_sp: $(size(expfi_sp))"
    @info "  Precomputed direction: $(size(direction))"

    return AngularGridData(nthmax, nphimax, theta_grid, phi_grid,
                          theta_weights, phi_weights,
                          plgd_sp, expfi_sp, expfi_sp_conj, direction)
end

"""
    compute_solid_angle_element(grid::AngularGridData, i_theta::Int, i_phi::Int) -> Float64

Compute solid angle element for integration.

# Arguments
- `grid::AngularGridData`: Angular grid
- `i_theta::Int`: θ grid index
- `i_phi::Int`: φ grid index

# Returns
- Solid angle element (steradian) for quadrature

# Implementation
Detects grid type by checking θ monotonicity:
- Gauss grids (θ = arccos(x)): decreasing θ → return w_θ * w_φ (sin(θ) implicit)
- Uniform grids: increasing θ → return sin(θ) * w_θ * w_φ (sin(θ) explicit)

# Usage
For surface integrals: ∫∫ f(θ,φ) dΩ ≈ Σᵢⱼ f(θᵢ,φⱼ) * compute_solid_angle_element(grid, i, j)
"""
function compute_solid_angle_element(grid::AngularGridData, i_theta::Int, i_phi::Int)
    θ = grid.theta_grid[i_theta]
    w_θ = grid.theta_weights[i_theta]
    w_φ = grid.phi_weights[i_phi]

    # Detect grid type: Gauss grids have decreasing θ (arccos is decreasing function)
    # Uniform grids have increasing θ
    is_gauss_grid = (grid.nthmax > 1) && (grid.theta_grid[2] < grid.theta_grid[1])

    if is_gauss_grid
        # Gauss-Legendre grid: sin(θ) already in weights
        return w_θ * w_φ
    else
        # Uniform grid: sin(θ) explicit
        return sin(θ) * w_θ * w_φ
    end
end

"""
    compute_direction_vectors(grid::AngularGridData) -> Array{Float64, 3}

Compute Cartesian unit vectors at all angular grid points.

# Returns
- `direction[1:3, nfi, nth]`: Cartesian components at (φ, θ) grid points

Already precomputed in grid.direction during initialization.
"""
function compute_direction_vectors(grid::AngularGridData)
    return grid.direction
end

"""
    print_grid_info(grid::AngularGridData)

Print diagnostic information about angular grid.
"""
function print_grid_info(grid::AngularGridData)
    println("Angular Grid Information:")
    println("  θ points: $(grid.nthmax) in [0, π]")
    println("  φ points: $(grid.nphimax) in [0, 2π)")
    println("  Total angular points: $(grid.nthmax * grid.nphimax)")
    println("  θ range: [$(@sprintf("%.4f", minimum(grid.theta_grid))), $(@sprintf("%.4f", maximum(grid.theta_grid)))]")
    println("  φ range: [$(@sprintf("%.4f", minimum(grid.phi_grid))), $(@sprintf("%.4f", maximum(grid.phi_grid)))]")

    # Compute total solid angle (should be 4π)
    total_solid_angle = 0.0
    for i in 1:grid.nthmax
        for j in 1:grid.nphimax
            total_solid_angle += compute_solid_angle_element(grid, i, j)
        end
    end
    println("  Total solid angle: $(@sprintf("%.4f", total_solid_angle)) (exact: $(4π))")
end

end  # module AngularGrid
