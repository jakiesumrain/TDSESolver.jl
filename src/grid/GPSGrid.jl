"""
    GPSGrid

Module for creating Generalized Pseudospectral (GPS) radial grids with algebraic mapping.

Implements the GPS method from Tong & Chu (1997) using Gauss-Lobatto-Legendre
interior points and algebraic coordinate mapping for efficient TDSE propagation.

# Key Design Decision
Uses Gauss-Lobatto interior points (excluding endpoints ±1) to ensure:
1. Dirichlet boundary conditions ψ(-1) = ψ(1) = 0 are naturally enforced
2. Spectral accuracy for the Hamiltonian eigenvalue problem
3. Eigenstates computed on this grid have exact orthonormality (no interpolation needed)
"""
module GPSGrid

using FastGaussQuadrature
using LinearAlgebra
using Printf  # For @sprintf

export GPSGridData, create_gps_grid, get_grid_spacing

"""
    GPSGridData

GPS radial grid with algebraic mapping using Gauss-Lobatto interior points.

# Fields
- `nrmax::Int`: Number of radial grid points (interior points only)
- `rmax::Float64`: Maximum radius (a.u.)
- `L::Float64`: Mapping parameter
- `α::Float64`: Smoothness parameter (α = 2L/rmax for exact mapping)
- `legendre_points::Vector{Float64}`: Gauss-Lobatto interior points x_i ∈ (-1, 1)
- `radial_grid::Vector{Float64}`: Mapped radial coordinates r_i
- `radial_derivative::Vector{Float64}`: Jacobian dr/dx for coordinate transform
- `quadrature_weights::Vector{Float64}`: Weights for radial integration (includes Jacobian)
- `gauss_weights::Vector{Float64}`: Bare Gauss-Lobatto weights (interior only)

# Mapping formula:
r(x) = L(1+x)/(1-x+α) where x ∈ (-1, 1)

This maps x=-1 → r=0 and x=1 → r=L(2)/(α) ≈ rmax
"""
struct GPSGridData
    nrmax::Int
    rmax::Float64
    L::Float64
    α::Float64
    legendre_points::Vector{Float64}
    radial_grid::Vector{Float64}
    radial_derivative::Vector{Float64}
    quadrature_weights::Vector{Float64}
    gauss_weights::Vector{Float64}  # Bare weights (interior only, no Jacobian)
end

"""
    create_gps_grid(nrmax::Int, rmax::Float64; L::Float64=30.0) -> GPSGridData

Create GPS radial grid with algebraic coordinate mapping.

Implements the Generalized Pseudospectral method following Tong & Chu (1997).
Uses Gauss-Lobatto-Legendre interior points (excluding boundary ±1) mapped to
radial coordinates [0, rmax] via algebraic transformation.

# Arguments
- `nrmax::Int`: Number of interior radial grid points
- `rmax::Float64`: Maximum radius in atomic units
- `L::Float64=30.0`: Mapping parameter (controls grid density near origin)

Note: α is computed automatically as α = 2L/rmax to ensure r(x=1) = rmax exactly.

# Returns
- `GPSGridData`: Complete grid specification

# Algorithm (Tong & Chu 1997):
1. Get N+2 Gauss-Lobatto points on [-1, 1] (includes endpoints)
2. Extract N interior points (excludes ±1) - this enforces ψ(0)=ψ(∞)=0
3. Apply algebraic mapping: r_i = L(1+x_i)/(1-x_i+α)
4. Compute Jacobian: dr/dx = L(2+α)/(1-x+α)²
5. Adjust quadrature weights: w_r,i = w_i * (dr/dx)_i

# Example
```julia
grid = create_gps_grid(200, 100.0, L=30.0)
println("Grid spans [0, \$(grid.rmax)] a.u. with \$(grid.nrmax) points")
```

# References
- Tong & Chu, Chem. Phys. 217, 119-130 (1997)
"""
function create_gps_grid(nrmax::Int, rmax::Float64; L::Float64=30.0)
    # Validate inputs
    if nrmax < 10
        error("nrmax must be at least 10, got $nrmax")
    end
    if rmax <= 0
        error("rmax must be positive, got $rmax")
    end
    if L <= 0
        error("Mapping parameter L must be positive, got $L")
    end

    # Compute α from rmax constraint: r(x=1) = L(2)/(0+α) = 2L/α = rmax
    # Therefore α = 2L/rmax
    α = 2.0 * L / rmax

    # Get Gauss-Lobatto points: need nrmax interior points
    # gausslobatto(N) returns N points including both endpoints
    # So we need N = nrmax + 2 total points to get nrmax interior points
    N_total = nrmax + 2
    x_all, w_all = gausslobatto(N_total)

    # Extract interior points (exclude endpoints at indices 1 and N_total)
    # This enforces Dirichlet boundary conditions ψ(-1) = ψ(1) = 0
    x_interior = x_all[2:end-1]
    w_interior = w_all[2:end-1]

    # Apply algebraic mapping: r(x) = L(1+x)/(1-x+α)
    radial_grid = zeros(Float64, nrmax)
    radial_derivative = zeros(Float64, nrmax)

    for i in 1:nrmax
        x = x_interior[i]
        denominator = 1.0 - x + α

        if denominator < 1e-10
            error("Mapping denominator too small at x=$x")
        end

        # Radial coordinate
        radial_grid[i] = L * (1.0 + x) / denominator

        # Jacobian dr/dx = L(2+α)/(1-x+α)²
        # Note: derivative of r = L(1+x)/(1-x+α) is:
        # dr/dx = L * [(1-x+α) + (1+x)] / (1-x+α)² = L(2+α)/(1-x+α)²
        radial_derivative[i] = L * (2.0 + α) / (denominator * denominator)
    end

    # Quadrature weights for radial integration
    # ∫f(r)dr = ∫f(r(x))(dr/dx)dx ≈ Σ w_i * f(r_i) * (dr/dx)_i
    quadrature_weights = w_interior .* radial_derivative

    # Verify grid properties
    r_actual_max = maximum(radial_grid)
    if abs(r_actual_max - rmax) > 0.05 * rmax
        @warn "Actual grid maximum r=$(@sprintf("%.2f", r_actual_max)) differs from requested rmax=$rmax by >5%"
    end

    # Check monotonicity
    for i in 2:nrmax
        if radial_grid[i] <= radial_grid[i-1]
            error("GPS grid not monotonic at index $i")
        end
    end

    return GPSGridData(nrmax, rmax, L, α, x_interior, radial_grid, radial_derivative, quadrature_weights, w_interior)
end

"""
    get_grid_spacing(grid::GPSGridData) -> Vector{Float64}

Compute local grid spacing Δr_i = r_{i+1} - r_i for diagnostics.

# Arguments
- `grid::GPSGridData`: GPS grid data

# Returns
- Vector of grid spacings (length nrmax-1)
"""
function get_grid_spacing(grid::GPSGridData)
    return diff(grid.radial_grid)
end

"""
    print_grid_info(grid::GPSGridData)

Print diagnostic information about GPS grid.
"""
function print_grid_info(grid::GPSGridData)
    spacing = get_grid_spacing(grid)
    println("GPS Grid Information:")
    println("  Points: $(grid.nrmax)")
    println("  Range: [0, $(@sprintf("%.2f", grid.rmax))] a.u.")
    println("  Mapping: L=$(grid.L), α=$(grid.α)")
    println("  Spacing: min=$(@sprintf("%.4f", minimum(spacing))) a.u., max=$(@sprintf("%.4f", maximum(spacing))) a.u.")
    println("  r[1]=$(@sprintf("%.6f", grid.radial_grid[1])) a.u. (near origin)")
    println("  r[end]=$(@sprintf("%.2f", grid.radial_grid[end])) a.u. (boundary)")
end

end  # module GPSGrid
