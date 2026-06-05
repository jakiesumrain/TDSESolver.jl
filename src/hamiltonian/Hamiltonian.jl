"""
    Hamiltonian

Module for field-free Hamiltonian construction and eigenstate calculation.

Constructs radial Hamiltonian operator in GPS discretization:
    H_l = -1/2 d²/dr² + l(l+1)/(2r²) + V(r)

For each angular momentum l, solves the eigenvalue problem:
    H_l φₙₗ(r) = Eₙₗ φₙₗ(r)

# Implementation Notes

## GPS Spectral Method
Uses GPS spectral eigenstate solver following Tong & Chu (1997) exactly.
Key features:
- Gauss-Lobatto-Legendre interior points (boundary ±1 excluded)
- Spectral differentiation via Canuto formula
- Boundary extraction enforces Dirichlet BCs: ψ(0)=ψ(∞)=0
- GPS coordinate transformation for radial grid

## NO INTERPOLATION
Eigenstates are computed DIRECTLY on the same grid used for propagation.
This ensures eigenvectors are exactly orthonormal in the Euclidean inner product,
which is REQUIRED for S-matrix unitarity.

## Eigenvector Orthonormality
CRITICAL: Eigenvectors from eigen(Symmetric(H)) are orthonormal in the
STANDARD EUCLIDEAN inner product: Σᵢ φₙ(xᵢ)φₘ(xᵢ) = δₙₘ

Do NOT renormalize with quadrature weights - this would break orthonormality
and cause the S-matrix to become non-unitary.

## References
- Tong & Chu, Chemical Physics 217 (1997) 119-130
- Canuto, Hussaini, Quarteroni, Zang, "Spectral Methods in Fluid Dynamics" (1988)
"""
module Hamiltonian

using LinearAlgebra
using Printf
using FastGaussQuadrature
using HDF5
using Dates

# Required types from other modules in parent scope
using ..GPSGrid: GPSGridData
using ..Potential: PotentialFunction

export HamiltonianData, solve_eigenstates, load_eigenstates_from_hdf5
export get_ground_state, get_eigenstate, print_hamiltonian_info
export save_eigenstates_to_hdf5

"""
    HamiltonianData

Container for Hamiltonian matrices and eigenstates.

# Fields
- `grid::GPSGridData`: Radial grid
- `potential::PotentialFunction`: Atomic potential V(r)
- `lmax::Int`: Maximum angular momentum
- `eigenvalues::Matrix{Float64}`: Energy eigenvalues [n, l+1]
- `eigenvectors::Array{Float64,3}`: Radial eigenfunctions [nrmax, n, l+1]
- `n_states::Vector{Int}`: Number of states per l channel
- `E_cutoff::Float64`: Energy cutoff for bound states (Ha)

# Indexing
- eigenvalues[n, l+1]: Energy of n-th state with angular momentum l
- eigenvectors[:, n, l+1]: Radial wavefunction φₙₗ(r)

# Important Note on Orthonormality
The eigenvectors are stored WITHOUT quadrature weight normalization.
They are orthonormal in the standard Euclidean inner product:
    Σᵢ φₙ(rᵢ) φₘ(rᵢ) = δₙₘ

This is REQUIRED for S-matrix unitarity in the propagator.
"""
struct HamiltonianData
    grid::GPSGridData
    potential::PotentialFunction
    lmax::Int
    eigenvalues::Matrix{Float64}
    eigenvectors::Array{Float64,3}
    n_states::Vector{Int}
    E_cutoff::Float64
end

"""
    build_D2_matrix_lobatto_interior(x_interior::Vector{Float64}, nrmax::Int) -> Matrix{Float64}

Build second derivative matrix D² at Gauss-Lobatto interior points.

The matrix is constructed by:
1. Building the full (N+2)×(N+2) D² matrix at all Gauss-Lobatto points (including ±1)
2. Extracting the interior N×N block (rows and columns 2 to N+1)

This automatically enforces Dirichlet boundary conditions ψ(-1) = ψ(1) = 0.

# Arguments
- `x_interior::Vector{Float64}`: The N interior Gauss-Lobatto points (from grid.legendre_points)
- `nrmax::Int`: Number of interior points (should equal length(x_interior))

# Returns
- `Matrix{Float64}`: N×N second derivative matrix at interior points
"""
function build_D2_matrix_lobatto_interior(x_interior::Vector{Float64}, nrmax::Int)
    # Reconstruct full Gauss-Lobatto points by adding endpoints
    N_total = nrmax + 2
    x_all = vcat(-1.0, x_interior, 1.0)

    # Build full first derivative matrix using Canuto formula
    D1_full = zeros(Float64, N_total, N_total)

    # c_i = 2 at endpoints, 1 at interior
    c = ones(N_total)
    c[1] = 2.0
    c[end] = 2.0

    # Polynomial degree N = N_total - 1
    N = N_total - 1

    # Off-diagonal elements
    for i in 1:N_total
        for j in 1:N_total
            if i != j
                D1_full[i,j] = (c[i]/c[j]) * (-1)^(i+j) / (x_all[i] - x_all[j])
            end
        end
    end

    # Diagonal elements
    D1_full[1,1] = -N*(N+1)/4.0
    D1_full[end,end] = N*(N+1)/4.0
    # Interior diagonals are zero for Gauss-Lobatto
    for i in 2:N_total-1
        D1_full[i,i] = 0.0
    end

    # Second derivative: D2 = D1 × D1
    D2_full = D1_full * D1_full

    # Extract interior-interior block (indices 2:N_total-1)
    D2_interior = D2_full[2:end-1, 2:end-1]

    return D2_interior
end

"""
    build_kinetic_matrix_on_gps_grid(grid::GPSGridData) -> Matrix{Float64}

Build kinetic energy matrix T = -½d²/dr² directly on the GPS grid points.

Uses the GPS transformation (Tong & Chu 1997) to convert the second derivative
from computational domain x ∈ (-1,1) to physical domain r(x).

# GPS Transformation
The kinetic operator transforms as:
    T = -½ d²/dr²

In spectral collocation, the second derivative matrix in r-space is:
    (D²_r)ᵢⱼ = (D²_x)ᵢⱼ / (r'(xᵢ) × r'(xⱼ))

where r'(x) = dr/dx is the Jacobian stored in grid.radial_derivative.
"""
function build_kinetic_matrix_on_gps_grid(grid::GPSGridData)
    N = grid.nrmax
    x = grid.legendre_points
    r_prime = grid.radial_derivative

    # Build D² in computational (x) space using Gauss-Lobatto formula
    D2_x = build_D2_matrix_lobatto_interior(x, N)

    # Transform to physical (r) space using GPS formula
    # (D²_r)ᵢⱼ = (D²_x)ᵢⱼ / (r'(xᵢ) × r'(xⱼ))
    D2_r = zeros(Float64, N, N)
    for i in 1:N
        for j in 1:N
            D2_r[i,j] = D2_x[i,j] / (r_prime[i] * r_prime[j])
        end
    end

    # Kinetic energy: T = -½ D²_r
    T = -0.5 * D2_r

    # Symmetrize to ensure numerical symmetry
    T = (T + T') / 2

    return T
end

"""
    build_hamiltonian_on_gps_grid(grid::GPSGridData, T::Matrix{Float64}, l::Int, V) -> Symmetric{Float64}

Build Hamiltonian matrix directly on GPS grid.

    [Ĥₗ]ᵢⱼ = Tᵢⱼ + [l(l+1)/(2rᵢ²) + V(rᵢ)]δᵢⱼ
"""
function build_hamiltonian_on_gps_grid(grid::GPSGridData, T::Matrix{Float64}, l::Int, V)
    N = grid.nrmax
    H = copy(T)

    for i in 1:N
        ri = grid.radial_grid[i]
        centrifugal = l * (l + 1) / (2.0 * ri^2)
        H[i, i] += centrifugal + V(ri)
    end

    return Symmetric(H)
end

"""
    solve_eigenstates(grid::GPSGridData,
                     potential::PotentialFunction,
                     lmax::Int;
                     n_max::Int=50,
                     E_cutoff::Float64=5.0) -> HamiltonianData

Solve eigenvalue problem for all angular momentum channels up to lmax.

# Arguments
- `grid::GPSGridData`: Radial GPS grid (eigenstates computed directly on this grid)
- `potential::PotentialFunction`: Atomic potential
- `lmax::Int`: Maximum angular momentum to compute
- `n_max::Int=50`: Maximum number of eigenstates per l channel
- `E_cutoff::Float64=5.0`: Energy cutoff for states (Ha)

# Returns
- `HamiltonianData`: Container with eigenvalues and eigenvectors

# Method
Computes eigenstates DIRECTLY on the input GPS grid (NO interpolation).
This ensures eigenvectors are orthonormal in the Euclidean inner product
on the same grid used for propagation, which is REQUIRED for S-matrix unitarity.

# Example
```julia
grid = create_gps_grid(200, 100.0)
pot = get_potential(:hydrogen)
ham = solve_eigenstates(grid, pot, 5, n_max=50, E_cutoff=5.0)

# Get ground state (n=1, l=0)
E_ground = ham.eigenvalues[1, 1]
φ_ground = ham.eigenvectors[:, 1, 1]
```
"""
function solve_eigenstates(grid::GPSGridData,
                          potential::PotentialFunction,
                          lmax::Int;
                          n_max::Int=50,
                          E_cutoff::Float64=5.0)
    nrmax = grid.nrmax

    @info "Computing eigenstates directly on GPS grid (no interpolation)" nrmax=nrmax lmax=lmax n_max=n_max E_cutoff=E_cutoff

    # Build kinetic energy matrix once (same for all l)
    T = build_kinetic_matrix_on_gps_grid(grid)

    @info "Kinetic matrix built on GPS grid" r_min=grid.radial_grid[1] r_max=grid.radial_grid[end]

    # Allocate storage for eigenvalues and eigenvectors
    eigenvalues = fill(Inf, n_max, lmax+1)
    eigenvectors = zeros(Float64, nrmax, n_max, lmax+1)
    n_states = zeros(Int, lmax+1)

    # Get potential function
    V = potential.V

    # Solve for each angular momentum channel
    for l in 0:lmax
        H = build_hamiltonian_on_gps_grid(grid, T, l, V)
        eigs = eigen(H)

        # Sort by eigenvalue
        perm = sortperm(eigs.values)
        E_sorted = eigs.values[perm]
        V_sorted = eigs.vectors[:, perm]

        # Keep states below E_cutoff
        n_kept = 0
        for n in 1:min(length(E_sorted), n_max)
            if E_sorted[n] < E_cutoff
                n_kept += 1
                eigenvalues[n_kept, l + 1] = E_sorted[n]
            else
                break
            end
        end

        n_states[l + 1] = n_kept

        if n_kept > 0
            # Store eigenvectors directly - they are already orthonormal
            # from eigen(Symmetric(H)) in Euclidean inner product
            # NO renormalization needed!
            for n in 1:n_kept
                eigenvectors[:, n, l + 1] .= V_sorted[:, n]
            end

            n_bound = count(eigenvalues[1:n_kept, l + 1] .< 0)
            @info "Solved l=$l" n_total=n_kept n_bound=n_bound E_min=round(eigenvalues[1,l+1], digits=6)
        end
    end

    return HamiltonianData(grid, potential, lmax, eigenvalues, eigenvectors,
                          n_states, E_cutoff)
end

"""
    get_ground_state(ham::HamiltonianData) -> Tuple{Float64, Vector{Float64}, Int, Int}

Extract ground state energy and wavefunction.

# Returns
- `E_ground::Float64`: Ground state energy (Ha)
- `φ_ground::Vector{Float64}`: Radial wavefunction
- `n_ground::Int`: Principal quantum number (usually 1)
- `l_ground::Int`: Angular momentum (usually 0)

# Example
```julia
ham = solve_eigenstates(grid, pot, 50)
E, φ, n, l = get_ground_state(ham)
println("Ground state: E = \$E Ha, (n=\$n, l=\$l)")
```
"""
function get_ground_state(ham::HamiltonianData)
    # Find the lowest energy state across all l channels
    E_min = Inf
    n_min = 1
    l_min = 0

    for l in 0:ham.lmax
        if ham.n_states[l+1] > 0
            E = ham.eigenvalues[1, l+1]
            if E < E_min
                E_min = E
                n_min = 1
                l_min = l
            end
        end
    end

    if isinf(E_min)
        error("No bound states found in Hamiltonian")
    end

    φ_ground = ham.eigenvectors[:, n_min, l_min+1]

    return E_min, φ_ground, n_min, l_min
end

"""
    get_eigenstate(ham::HamiltonianData, n::Int, l::Int) -> Tuple{Float64, Vector{Float64}}

Extract specific eigenstate by quantum numbers.

# Arguments
- `ham::HamiltonianData`: Hamiltonian with solved eigenstates
- `n::Int`: Principal quantum number (1-indexed)
- `l::Int`: Angular momentum quantum number

# Returns
- `E::Float64`: Eigenvalue (Ha)
- `φ::Vector{Float64}`: Radial eigenfunction

# Example
```julia
# Get 2p state (n=2, l=1)
E_2p, φ_2p = get_eigenstate(ham, 2, 1)
```
"""
function get_eigenstate(ham::HamiltonianData, n::Int, l::Int)
    if l > ham.lmax || l < 0
        error("l=$l out of range [0, $(ham.lmax)]")
    end

    if n < 1 || n > ham.n_states[l+1]
        error("n=$n out of range [1, $(ham.n_states[l+1])] for l=$l")
    end

    E = ham.eigenvalues[n, l+1]
    φ = ham.eigenvectors[:, n, l+1]

    return E, φ
end

"""
    print_hamiltonian_info(ham::HamiltonianData)

Print summary of Hamiltonian eigenstates.
"""
function print_hamiltonian_info(ham::HamiltonianData)
    println("Hamiltonian Information:")
    println("  Atom type: $(ham.potential.atom_type)")
    println("  Grid points: $(ham.grid.nrmax)")
    println("  Max angular momentum: $(ham.lmax)")
    println("  Energy cutoff: $(ham.E_cutoff) Ha")
    println()
    println("  States per l:")

    total_states = 0
    total_bound = 0

    for l in 0:ham.lmax
        n_total = ham.n_states[l+1]
        if n_total > 0
            n_bound = count(ham.eigenvalues[1:n_total, l+1] .< 0)
            n_cont = n_total - n_bound
            total_states += n_total
            total_bound += n_bound
            E_min = ham.eigenvalues[1, l+1]
            E_max = ham.eigenvalues[n_total, l+1]
            @printf("    l=%d: %d states (%d bound, %d continuum), E ∈ [%.6f, %.4f] Ha\n",
                    l, n_total, n_bound, n_cont, E_min, E_max)
        end
    end

    # Print ground state
    E_ground, _, n_ground, l_ground = get_ground_state(ham)
    println()
    println("  Ground state: E = $(@sprintf("%.6f", E_ground)) Ha (n=$n_ground, l=$l_ground)")
    println("  Ionization potential: $(@sprintf("%.6f", -E_ground)) Ha ($(@sprintf("%.2f", -E_ground * 27.2114)) eV)")
    println()
    println("  Total: $total_states states ($total_bound bound, $(total_states - total_bound) discretized continuum)")
end

"""
    save_eigenstates_to_hdf5(filename::String, ham::HamiltonianData)

Save Hamiltonian eigenstates to HDF5 file for later reuse.

# HDF5 File Structure
```
eigenstates.h5
├── /metadata
│   ├── version (String)
│   ├── atom_type (String)
│   ├── lmax (Int)
│   ├── n_max (Int)
│   ├── E_cutoff (Float64)
│   └── creation_time (String)
├── /grid
│   ├── nrmax (Int)
│   ├── radial_grid (Float64[nrmax])
│   └── quadrature_weights (Float64[nrmax])
└── /eigenstates
    ├── eigenvalues (Float64[n_max, lmax+1])
    ├── eigenvectors (Float64[nrmax, n_max, lmax+1])
    └── n_states (Int[lmax+1])
```
"""
function save_eigenstates_to_hdf5(filename::String, ham::HamiltonianData)
    @info "Saving eigenstates to HDF5" filename=filename

    h5open(filename, "w") do file
        # Metadata
        g_meta = create_group(file, "metadata")
        g_meta["version"] = "1.0"
        g_meta["atom_type"] = String(ham.potential.atom_type)
        g_meta["lmax"] = ham.lmax
        g_meta["n_max"] = size(ham.eigenvalues, 1)
        g_meta["E_cutoff"] = ham.E_cutoff
        g_meta["creation_time"] = string(Dates.now())

        # Grid
        g_grid = create_group(file, "grid")
        g_grid["nrmax"] = ham.grid.nrmax
        g_grid["radial_grid"] = ham.grid.radial_grid
        g_grid["quadrature_weights"] = ham.grid.quadrature_weights

        # Eigenstates
        g_eig = create_group(file, "eigenstates")
        g_eig["eigenvalues"] = ham.eigenvalues
        g_eig["eigenvectors"] = ham.eigenvectors
        g_eig["n_states"] = ham.n_states
    end

    filesize_mb = round(filesize(filename) / 1024^2, digits=2)
    @info "Eigenstates saved successfully" filesize="$(filesize_mb) MB"
end

"""
    load_eigenstates_from_hdf5(filename::String, grid::GPSGridData,
                                potential::PotentialFunction) -> HamiltonianData

Load Hamiltonian eigenstates from HDF5 file.

# Arguments
- `filename::String`: Path to HDF5 file
- `grid::GPSGridData`: Grid for the simulation (must match saved grid)
- `potential::PotentialFunction`: Potential function (for consistency)

# Returns
- `HamiltonianData`: Reconstructed Hamiltonian with loaded eigenstates

# Throws
- Error if grid dimensions don't match
- Error if file is not found or corrupted
"""
function load_eigenstates_from_hdf5(filename::String, grid::GPSGridData,
                                     potential::PotentialFunction)
    @info "Loading eigenstates from HDF5" filename=filename

    if !isfile(filename)
        error("Eigenstates file not found: $filename")
    end

    ham = h5open(filename, "r") do file
        # Read metadata
        atom_type = Symbol(read(file["metadata/atom_type"]))
        lmax = read(file["metadata/lmax"])
        E_cutoff = read(file["metadata/E_cutoff"])

        # Read grid info and verify
        nrmax_file = read(file["grid/nrmax"])
        if nrmax_file != grid.nrmax
            error("Grid mismatch: file has nrmax=$(nrmax_file), but provided grid has nrmax=$(grid.nrmax)")
        end

        # Read eigenstates
        eigenvalues = read(file["eigenstates/eigenvalues"])
        eigenvectors = read(file["eigenstates/eigenvectors"])
        n_states = read(file["eigenstates/n_states"])

        return HamiltonianData(grid, potential, lmax, eigenvalues, eigenvectors,
                              n_states, E_cutoff)
    end

    total_states = sum(ham.n_states)
    total_bound = sum(count(ham.eigenvalues[1:ham.n_states[l+1], l+1] .< 0) for l in 0:ham.lmax)

    @info "Eigenstates loaded successfully" atom=ham.potential.atom_type lmax=ham.lmax total_states=total_states bound_states=total_bound

    return ham
end

end  # module Hamiltonian
