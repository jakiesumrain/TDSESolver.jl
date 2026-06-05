"""
    Propagator

Module for time-dependent Schrödinger equation propagation.

Implements split-operator method in energy representation following
Tong & Chu (1997) algorithm with GPS discretization.

# Algorithm
Three-step split-operator propagation:
    ψ(t+Δt) ≈ exp(-iĤ₀Δt/2) exp(-iV̂_intΔt) exp(-iĤ₀Δt/2) ψ(t)

where:
- Ĥ₀ = T̂ + V̂_atom: Field-free Hamiltonian (via S-matrix)
- V̂_int = -r⃗·E⃗(t): Atom-field interaction (length gauge)

# Implementation Notes

## Fortran Blueprint Fidelity
Follows D_inner_out_volkov_3d_with_prob.f90 algorithm structure:
- S-matrix propagation in energy eigenbasis
- Radial ↔ angular transformations for field interaction
- Length gauge dipole interaction

## Performance
Main computational cost: Radial-angular transformations O(nrmax × nθ × nφ × lmax²)
"""
module Propagator

using LinearAlgebra
using Printf

# Required types from other modules
using ..GPSGrid: GPSGridData
using ..AngularGrid: AngularGridData
using ..Wavefunction: WavefunctionData, compute_norm
using ..Hamiltonian: HamiltonianData
using ..CoordinateTransform: transform_radial_to_angular!, transform_angular_to_radial!

export PropagatorData, create_propagator, propagate_step!
export construct_s_matrix, apply_s_matrix!, field_free_propagate!
export propagate_step_with_absorber!

"""
    PropagatorData

Container for time propagation data and operators.

# Fields
- `ham::HamiltonianData`: Hamiltonian with eigenstates
- `grid::GPSGridData`: Radial grid
- `angular_grid::AngularGridData`: Angular grid for transformations
- `s_matrices::Vector{Matrix{ComplexF64}}`: S-matrices for each l [lmax+1]
- `dt::Float64`: Time step (atomic units)
- `energy_cutoff::Float64`: Energy cutoff for S-matrix (Ha)

# S-Matrix Formula
For each angular momentum l:
    S_{l,ij}(Δt) = Σₙ φₙₗ(rᵢ) φₙₗ(rⱼ) exp(-iEₙₗΔt/2) wⱼ

where φₙₗ are radial eigenfunctions, Eₙₗ are eigenvalues, wⱼ are GPS weights.
"""
struct PropagatorData
    ham::HamiltonianData
    grid::GPSGridData
    angular_grid::AngularGridData
    s_matrices::Vector{Matrix{ComplexF64}}
    dt::Float64
    energy_cutoff::Float64
end

"""
    construct_s_matrix(ham::HamiltonianData, l::Int, dt::Float64) -> Matrix{ComplexF64}

Construct S-matrix for field-free propagation for angular momentum l.

# Arguments
- `ham::HamiltonianData`: Hamiltonian with solved eigenstates
- `l::Int`: Angular momentum quantum number
- `dt::Float64`: Time step (atomic units)

# Returns
- `Matrix{ComplexF64}`: S-matrix [nrmax, nrmax]

# Physics
The S-matrix represents field-free time evolution in the energy eigenbasis:
    S_l(Δt/2) = Σₙ |φₙₗ⟩⟨φₙₗ| exp(-iEₙₗΔt/2)

In matrix form for propagation:
    g'(rᵢ) = Σⱼ S_{l,ij} g(rⱼ)

# Unitarity Requirement
The S-matrix is UNITARY if and only if the eigenstates form a COMPLETE basis.
This requires including ALL eigenstates (no energy cutoff in S-matrix).

Eigenstates are computed directly on the GPS grid and are orthonormal in the
Euclidean inner product: Σᵢ φₙ(rᵢ)φₘ(rᵢ) = δₙₘ

# Example
```julia
S0 = construct_s_matrix(ham, 0, 0.1)  # l=0, Δt=0.1 a.u.
```
"""
function construct_s_matrix(ham::HamiltonianData, l::Int, dt::Float64)
    if l > ham.lmax || l < 0
        error("l=$l out of range [0, $(ham.lmax)]")
    end

    nrmax = ham.grid.nrmax
    n_states = ham.n_states[l+1]

    # Initialize S-matrix
    S = zeros(ComplexF64, nrmax, nrmax)

    # Sum over ALL eigenstates for complete basis (required for unitarity)
    # S = Σₙ φₙ φₙᵀ exp(-iEₙΔt/2)
    for n in 1:n_states
        E = ham.eigenvalues[n, l+1]

        # Time evolution phase
        phase = exp(-im * E * dt / 2.0)

        # Eigenvector
        φ = ham.eigenvectors[:, n, l+1]

        # Outer product (no weights - eigenstates are orthonormal in Euclidean inner product)
        for i in 1:nrmax
            for j in 1:nrmax
                S[i, j] += φ[i] * φ[j] * phase
            end
        end
    end

    return S
end

"""
    create_propagator(ham::HamiltonianData,
                     angular_grid::AngularGridData,
                     dt::Float64) -> PropagatorData

Create propagator with pre-computed S-matrices for all l channels.

# Arguments
- `ham::HamiltonianData`: Hamiltonian with solved eigenstates
- `angular_grid::AngularGridData`: Angular grid for transformations
- `dt::Float64`: Time step (atomic units)

# Returns
- `PropagatorData`: Propagator ready for time stepping

# Important: Complete Basis Requirement
The S-matrix is unitary ONLY when ALL eigenstates are included.
The Hamiltonian solver must use n_max >= nrmax to ensure completeness.

# Example
```julia
ham = solve_eigenstates(grid, pot, 50, n_max=grid.nrmax)  # Complete basis!
ang_grid = create_angular_grid("uniform", 90, 60)
prop = create_propagator(ham, ang_grid, 0.1)
```

# Performance Note
S-matrix construction is expensive: O(lmax × nrmax² × n_states)
Pre-computation amortizes cost over many time steps.
"""
function create_propagator(ham::HamiltonianData,
                          angular_grid::AngularGridData,
                          dt::Float64)
    @info "Constructing propagator" lmax=ham.lmax dt=dt

    # Construct S-matrix for each l
    s_matrices = Vector{Matrix{ComplexF64}}(undef, ham.lmax + 1)

    for l in 0:ham.lmax
        s_matrices[l+1] = construct_s_matrix(ham, l, dt)

        n_states = ham.n_states[l+1]
        nrmax = ham.grid.nrmax
        if n_states > 0
            if n_states < nrmax
                @warn "Incomplete basis for l=$l: n_states=$n_states < nrmax=$nrmax. S-matrix may not be unitary!"
            end
            @info "Constructed S-matrix for l=$l" size=size(s_matrices[l+1]) n_states=n_states
        end
    end

    # Note: energy_cutoff is kept in PropagatorData for compatibility but not used in S-matrix
    return PropagatorData(ham, ham.grid, angular_grid, s_matrices, dt, Inf)
end

"""
    apply_s_matrix!(wfn::WavefunctionData, prop::PropagatorData)

Apply S-matrix propagation for field-free evolution: exp(-iĤ₀Δt/2).

# Arguments
- `wfn::WavefunctionData`: Wavefunction (modified in-place)
- `prop::PropagatorData`: Propagator with S-matrices

# Algorithm
For each (l, m) channel:
    g'(rᵢ, m, l) = Σⱼ S_{l,ij} g(rⱼ, m, l)

The S-matrix is independent of m (azimuthal symmetry), so same S_l
applied to all m values for given l.

# Performance
- Complexity: O(lmax × (2lmax+1) × nrmax²)
- Can parallelize over l channels (OpenMP/threads)
"""
function apply_s_matrix!(wfn::WavefunctionData, prop::PropagatorData)
    nrmax = wfn.nrmax
    lmax = min(wfn.lmax, prop.ham.lmax)

    # Apply S-matrix for each l channel
    for l in 0:lmax
        S = prop.s_matrices[l+1]

        # Apply to all m values for this l
        for m in -l:l
            m_index = m + wfn.lmax + 1

            # Matrix-vector multiplication: g' = S * g
            g_in = wfn.g[:, m_index, l+1]
            g_out = S * g_in

            # Store result
            wfn.g[:, m_index, l+1] .= g_out
        end
    end
end

"""
    field_free_propagate!(wfn::WavefunctionData, prop::PropagatorData)

Perform full field-free propagation: exp(-iĤ₀Δt).

Applies S-matrix twice: exp(-iĤ₀Δt/2) × 2

# Arguments
- `wfn::WavefunctionData`: Wavefunction (modified in-place)
- `prop::PropagatorData`: Propagator

# Use Case
For testing or field-free evolution between laser pulses.
"""
function field_free_propagate!(wfn::WavefunctionData, prop::PropagatorData)
    # Apply S-matrix twice for full time step
    apply_s_matrix!(wfn, prop)
    apply_s_matrix!(wfn, prop)
end

"""
    propagate_step!(wfn::WavefunctionData, prop::PropagatorData, E_field::Vector{Float64}, t::Float64)

Perform one split-operator propagation step.

# Arguments
- `wfn::WavefunctionData`: Wavefunction (modified in-place)
- `prop::PropagatorData`: Propagator with S-matrices
- `E_field::Vector{Float64}`: Electric field vector [Ex, Ey, Ez] (a.u.)
- `t::Float64`: Current time (for field evaluation)

# Algorithm
Split-operator method (3 steps):
1. exp(-iĤ₀Δt/2): Apply S-matrix (field-free propagation)
2. exp(-iV̂_intΔt): Apply field interaction in coordinate space
   - Transform radial → angular
   - Apply dipole interaction: exp(-i r⃗·E⃗ Δt)
   - Transform angular → radial
3. exp(-iĤ₀Δt/2): Apply S-matrix again

# Length Gauge Interaction
    V_int(r,θ,φ,t) = -r⃗·E⃗(t)
                    = -r[Ex sin(θ)cos(φ) + Ey sin(θ)sin(φ) + Ez cos(θ)]

For linearly polarized field along z: V_int = -r Ez cos(θ)

# Example
```julia
# Single time step
E_field = [0.0, 0.0, 0.05]  # Field along z, 0.05 a.u.
propagate_step!(wfn, prop, E_field, t)
```

# Implementation
Uses full angular momentum coupling via Clebsch-Gordan coefficients
and dipole selection rules (Δl = ±1, Δm = 0 for z-polarization).
"""
function propagate_step!(wfn::WavefunctionData, prop::PropagatorData,
                        E_field::Vector{Float64}, t::Float64)
    # Step 1: Field-free propagation exp(-iĤ₀Δt/2)
    apply_s_matrix!(wfn, prop)

    # Step 2: Field interaction exp(-iV_intΔt)
    # Apply full dipole coupling with angular momentum selection rules
    apply_field_interaction_full!(wfn, prop, E_field)

    # Step 3: Field-free propagation exp(-iĤ₀Δt/2)
    apply_s_matrix!(wfn, prop)
end

"""
    apply_field_interaction_simple!(wfn::WavefunctionData, prop::PropagatorData, E_field::Vector{Float64})

Apply field interaction in simplified form (diagonal in radial basis).

# Note
This is a placeholder for the full angular transformation version.
For z-polarized field: V_int(r,l,m) ≈ -r Ez ⟨lm|cos(θ)|lm⟩

Full implementation requires angular grid transformation (future work).
"""
function apply_field_interaction_simple!(wfn::WavefunctionData, prop::PropagatorData,
                                         E_field::Vector{Float64})
    # Extract field components
    Ex, Ey, Ez = E_field

    # For z-polarized field (most common case)
    # Apply selection rules: Δl = ±1, Δm = 0
    # This couples l channels via dipole matrix elements

    # Simplified approach: Apply average field interaction
    # Full version requires angular momentum coupling

    # For now, apply diagonal approximation
    # V_int(r) ≈ -r * |E|
    E_magnitude = sqrt(Ex^2 + Ey^2 + Ez^2)

    for l in 0:wfn.lmax
        for m in -l:l
            m_index = m + wfn.lmax + 1
            l_index = l + 1

            for ir in 1:wfn.nrmax
                r = prop.grid.radial_grid[ir]

                # Dipole interaction phase
                # For z-polarized: ⟨lm|r cos(θ)|lm⟩ = 0 (no diagonal terms)
                # Need off-diagonal coupling (future work)

                # Placeholder: Apply small phase based on field strength
                # This maintains stability but doesn't implement full physics
                phase = -im * r * E_magnitude * prop.dt * 0.1  # Small coupling

                wfn.g[ir, m_index, l_index] *= exp(phase)
            end
        end
    end

    @warn "Using simplified field interaction (no angular coupling). Full implementation pending." maxlog=1
end

"""
    apply_field_interaction_full!(wfn::WavefunctionData, prop::PropagatorData, E_field::Vector{Float64})

Apply full 3D dipole field interaction using angular coordinate transformation.

# Physics
V̂_int = -r⃗·E⃗(t) = -r[Ex sin(θ)cos(φ) + Ey sin(θ)sin(φ) + Ez cos(θ)]

Supports:
- Linear polarization (Ex, Ey, or Ez)
- Elliptical polarization (Ex + Ey with phase)
- Circular polarization (Ex = Ey with π/2 phase)
- Bicircular fields (2ω y-component)

# Algorithm (matching Fortran lines 635-708)
1. Transform radial → angular: g(r,m,l) → ψ(φ,θ,r)
2. Apply field interaction in coordinate space: ψ → ψ × exp(-i r⃗·E⃗ Δt)
3. Transform angular → radial: ψ(φ,θ,r) → g(r,m,l)

# Fortran Reference
- Lines 638-655: Radial → Angular transformation
- Line 663: Field interaction exp(-i E⃗·r⃗ Δt)
- Lines 688-708: Angular → Radial transformation

# Implementation Notes
This is the EXACT algorithm from Fortran, supporting full elliptical polarization.
Much more accurate than angular momentum coupling in radial basis for non-z-polarized fields.
"""
function apply_field_interaction_full!(wfn::WavefunctionData, prop::PropagatorData,
                                       E_field::Vector{Float64})
    # Extract field components
    Ex, Ey, Ez = E_field

    # Skip if field is negligible
    E_magnitude = sqrt(Ex^2 + Ey^2 + Ez^2)
    if E_magnitude < 1e-15
        return
    end

    nrmax = wfn.nrmax
    nthmax = prop.angular_grid.nthmax
    nfimax = prop.angular_grid.nphimax

    # Allocate angular coordinate wavefunction (Fortran: psai_sp(nfi,nth,nr))
    psai_sp = zeros(ComplexF64, nfimax, nthmax, nrmax)

    # Step 1: Transform radial → angular (Fortran lines 638-655)
    transform_radial_to_angular!(psai_sp, wfn, prop.angular_grid)

    # Step 2: Apply field interaction in angular coordinates (Fortran lines 660-668)
    # Fortran line 663: E_dot_r_t = r*(direction(1)*Ex + direction(2)*Ey + direction(3)*Ez)*dt
    #                   psai_sp *= exp(-i * E_dot_r_t)

    direction = prop.angular_grid.direction  # [1:3, nfi, nth]

    for nr in 1:nrmax
        r = prop.grid.radial_grid[nr]

        for nth in 1:nthmax
            for nfi in 1:nfimax
                # Compute r⃗·E⃗ = r * (direction_x * Ex + direction_y * Ey + direction_z * Ez)
                # Fortran line 663
                E_dot_r = r * (direction[1, nfi, nth] * Ex +
                              direction[2, nfi, nth] * Ey +
                              direction[3, nfi, nth] * Ez)

                # Apply phase factor: exp(-i * E⃗·r⃗ * Δt)
                # Fortran line 664: psai_sp *= zexp((0,  -1) * E_dot_r_t)
                phase = -E_dot_r * prop.dt
                psai_sp[nfi, nth, nr] *= exp(im * phase)
            end
        end
    end

    # Step 3: Transform angular → radial (Fortran lines 688-708)
    transform_angular_to_radial!(wfn, psai_sp, prop.angular_grid, prop.grid)
end

"""
    print_propagator_info(prop::PropagatorData)

Print diagnostic information about propagator.
"""
function print_propagator_info(prop::PropagatorData)
    println("Propagator Information:")
    println("  Time step: $(@sprintf("%.4f", prop.dt)) a.u. ($(@sprintf("%.2f", prop.dt * 24.2)) as)")
    println("  Energy cutoff: $(prop.energy_cutoff) Ha")
    println("  Grid points: $(prop.grid.nrmax)")
    println("  Max angular momentum: $(prop.ham.lmax)")
    println("  Angular grid: $(prop.angular_grid.ntheta) × $(prop.angular_grid.nphi)")

    println("\n  S-matrices per l channel:")
    total_states = 0
    for l in 0:prop.ham.lmax
        n = prop.ham.n_states[l+1]
        if n > 0
            mem_mb = sizeof(prop.s_matrices[l+1]) / 1024^2
            println("    l=$l: $n states, S-matrix size $(size(prop.s_matrices[l+1])), $(@sprintf("%.2f", mem_mb)) MB")
            total_states += n
        end
    end

    total_mem = sum(sizeof(S) for S in prop.s_matrices) / 1024^2
    println("\n  Total eigenstates: $total_states")
    println("  Total S-matrix memory: $(@sprintf("%.2f", total_mem)) MB")
end

"""
    propagate_step_with_absorber!(wfn::WavefunctionData, prop::PropagatorData,
                                  E_field::Vector{Float64}, t::Float64,
                                  absorber_mask::Vector{Float64})

Perform one split-operator propagation step WITH absorbing boundary.

This is the correct function to use for HHG calculations where boundary
reflections must be eliminated.

# Arguments
- `wfn::WavefunctionData`: Wavefunction (modified in-place)
- `prop::PropagatorData`: Propagator with S-matrices
- `E_field::Vector{Float64}`: Electric field vector [Ex, Ey, Ez] (a.u.)
- `t::Float64`: Current time (for field evaluation)
- `absorber_mask::Vector{Float64}`: Pre-computed absorbing mask from AbsorbingBoundary

# Algorithm
Split-operator method (4 steps):
1. exp(-iĤ₀Δt/2): Apply S-matrix (field-free propagation)
2. exp(-iV̂_intΔt): Apply field interaction in coordinate space
3. exp(-iĤ₀Δt/2): Apply S-matrix again
4. Apply absorbing boundary mask (HHG-essential)

# Fortran Reference
From rescatteing+hhg-he.f90:
- Lines 756-764: S-matrix application
- Lines 766-776: Absorbing boundary (cos^(1/4) mask)

# Example
```julia
# Create absorber
absorber = create_absorbing_boundary(grid, r0=100.0)

# Propagation loop
for step in 1:n_steps
    E_field = [0.0, 0.0, E0 * sin(omega * t)]
    propagate_step_with_absorber!(wfn, prop, E_field, t, absorber.mask)
    t += dt
end
```
"""
function propagate_step_with_absorber!(wfn::WavefunctionData, prop::PropagatorData,
                                       E_field::Vector{Float64}, t::Float64,
                                       absorber_mask::Vector{Float64})
    # Step 1: Field-free propagation exp(-iĤ₀Δt/2)
    apply_s_matrix!(wfn, prop)

    # Step 2: Field interaction exp(-iV_intΔt)
    apply_field_interaction_full!(wfn, prop, E_field)

    # Step 3: Field-free propagation exp(-iĤ₀Δt/2)
    apply_s_matrix!(wfn, prop)

    # Step 4: Apply absorbing boundary (essential for HHG)
    # Matches Fortran lines 766-776
    nrmax = wfn.nrmax
    lmax = wfn.lmax

    for l in 0:lmax
        l_index = l + 1
        for m in -l:l
            m_index = m + lmax + 1

            for ir in 1:nrmax
                wfn.g[ir, m_index, l_index] *= absorber_mask[ir]
            end
        end
    end
end

end  # module Propagator
