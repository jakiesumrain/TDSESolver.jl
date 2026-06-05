"""
Volkov State Projection Module

Projects outer region wavefunction onto momentum eigenstates using regular Volkov states
(free particle in laser field, no Coulomb correction).

Based on PRA 74, 031405(R) (2006) - Region splitting with Volkov projection

**Theory:**

For outer region wavefunction ψ_outer(r,θ,φ,t), project onto momentum state |p⟩:

    C(p,l,m,t) = ⟨p,l,m|ψ_outer(t)⟩ = ∫ ψ_outer(r,m,l,t) j_l(pr) r² dr

Where j_l(pr) is the spherical Bessel function.

**Volkov Phase Evolution:**

After projection at time t_i, accumulate Volkov phase to end of pulse:

    C̄(p,t_final) = exp(-i ∫[p²/2 + A(t)²/2 + p·A(t)] dt) C(p,t_i)

**Momentum Distribution:**

    dP/dE dΩ = 2E |∑_i C̄(p,t_i)|²

Summing coherently over all splitting times t_i.

**Usage:**

```julia
using .VolkovProjection

# Setup momentum grid
p_grid = create_momentum_grid(p_max=3.0, n_p=100, n_theta=45, n_phi=24)

# Initialize projector
projector = create_volkov_projector(radial_grid, angular_grid, p_grid)

# At each splitting time
project_volkov!(projector, wfn_outer, A_vector, time)

# At end of propagation
finalize_volkov_projection!(projector, A_final)

# Compute momentum distribution
P_p, P_E = compute_momentum_distribution(projector)
```

**References:**
- PRA 74, 031405(R) (2006): Eq. 6-8, region splitting with Volkov
- Fortran: D_inner_out_volkov_3d_with_prob.f90, lines 909-923 (projection), 1013-1015 (phase)
"""
module VolkovProjection

export MomentumGrid, VolkovProjectorData
export create_momentum_grid
export create_volkov_projector
export project_volkov!
export finalize_volkov_projection!
export compute_momentum_distribution
export transform_to_momentum_space

using LinearAlgebra
using Bessels  # For spherical Bessel functions

"""
    MomentumGrid

Momentum space grid in spherical coordinates.

**Fields:**
- `p::Vector{Float64}`: Momentum magnitude grid (n_p points)
- `theta::Vector{Float64}`: Polar angle grid (n_theta points, 0 to π)
- `phi::Vector{Float64}`: Azimuthal angle grid (n_phi points, 0 to 2π)
- `n_p::Int`: Number of momentum points
- `n_theta::Int`: Number of theta points
- `n_phi::Int`: Number of phi points
- `p_max::Float64`: Maximum momentum
"""
struct MomentumGrid
    p::Vector{Float64}
    theta::Vector{Float64}
    phi::Vector{Float64}
    n_p::Int
    n_theta::Int
    n_phi::Int
    p_max::Float64
end

"""
    VolkovProjectorData

Data structure for Volkov state projection (FORTRAN-MATCHING).

**Fields:**
- `psai_pt::Array{ComplexF64,3}`: Accumulated momentum space wavefunction [iphi, itheta, ip]
- `momentum_grid::MomentumGrid`: Momentum space grid
- `angular_grid`: Angular grid for transformations
- `n_projections::Int`: Number of projection times
- `radial_grid`: Reference to radial grid
- `lmax::Int`: Maximum angular momentum
- `A_accumulated::Vector{Float64}`: Accumulated vector potential integral
- `A_square_accumulated::Float64`: Accumulated A² integral
- `mp_init::Int`: Initial time step for this split sequence

**Algorithm (Fortran lines 784-1024):**
Each split:
1. Compute fresh C_lm from g_split (no accumulation)
2. Transform C_lm → ψ(φ,θ,p)
3. Apply Volkov phase
4. Accumulate: psai_pt += ψ_phased
"""
mutable struct VolkovProjectorData
    psai_pt::Array{ComplexF64,3}  # [iphi, itheta, ip] - accumulated momentum space wavefunction
    momentum_grid::MomentumGrid
    angular_grid::Any  # AngularGridData
    n_projections::Int
    radial_grid::Any  # GPSGridData
    lmax::Int
    A_accumulated::Vector{Float64}  # Accumulated vector potential
    A_square_accumulated::Float64   # Accumulated A² integral
    mp_init::Int  # Initial time step
end

"""
    create_momentum_grid(p_max=3.0, n_p=100, n_theta=45, n_phi=24) -> MomentumGrid

Create momentum space grid in spherical coordinates.

**Parameters:**
- `p_max`: Maximum momentum (default 3.0 a.u.)
  - For hydrogen Ip=0.5 a.u., Up=0.001 a.u.: p_max ≈ √(2(Ip+10Up)) ≈ 1 a.u.
  - Use 2-3× for safety margin
- `n_p`: Number of momentum points (default 100)
- `n_theta`: Number of polar angle points (default 45)
- `n_phi`: Number of azimuthal points (default 24)

**Returns:**
- `MomentumGrid`: Momentum space grid

**Example:**
```julia
# Standard grid for moderate intensity
p_grid = create_momentum_grid(p_max=2.0, n_p=80, n_theta=45, n_phi=24)

# High resolution for high-energy electrons
p_grid = create_momentum_grid(p_max=5.0, n_p=200, n_theta=90, n_phi=48)
```
"""
function create_momentum_grid(;p_max::Float64=3.0,
                               n_p::Int=100,
                               n_theta::Int=45,
                               n_phi::Int=24)
    # Momentum magnitude: uniform grid from 0 to p_max
    p = range(0.0, p_max, length=n_p) |> collect

    # Polar angle: uniform from 0 to π
    theta = range(0.0, π, length=n_theta) |> collect

    # Azimuthal angle: uniform from 0 to 2π
    phi = range(0.0, 2π, length=n_phi) |> collect

    @info "Created momentum grid" p_max n_p n_theta n_phi

    return MomentumGrid(p, theta, phi, n_p, n_theta, n_phi, p_max)
end

"""
    spherical_bessel_j(l, x) -> Float64

Compute spherical Bessel function j_l(x).

Uses Bessels.jl: sphericalbesselj(l, x)

For x ≈ 0: j_0(0) = 1, j_l>0(0) = 0
"""
function spherical_bessel_j(l::Int, x::Real)
    if abs(x) < 1e-10
        # Handle x → 0 limit
        return l == 0 ? 1.0 : 0.0
    else
        # Use Bessels.jl spherical Bessel function
        return sphericalbesselj(l, x)
    end
end

"""
    create_volkov_projector(radial_grid, lmax, momentum_grid, angular_grid) -> VolkovProjectorData

Create Volkov projector for momentum space projection.

**Parameters:**
- `radial_grid`: GPSGridData
- `lmax`: Maximum angular momentum
- `momentum_grid`: MomentumGrid
- `angular_grid`: AngularGridData (for spherical harmonic transformations)

**Returns:**
- `VolkovProjectorData`: Initialized projector with psai_pt = 0
"""
function create_volkov_projector(radial_grid, lmax, momentum_grid, angular_grid)
    # Initialize accumulated momentum space wavefunction to zero
    # Dimensions: [iphi, itheta, ip] matching Fortran psai_pt
    psai_pt = zeros(ComplexF64,
                    angular_grid.nphimax,
                    angular_grid.nthmax,
                    momentum_grid.n_p)

    # Initialize accumulation variables
    A_accumulated = zeros(Float64, 3)
    A_square_accumulated = 0.0
    mp_init = 0

    @info "Created Volkov projector" lmax n_p=momentum_grid.n_p ntheta=angular_grid.nthmax nphi=angular_grid.nphimax

    return VolkovProjectorData(psai_pt, momentum_grid, angular_grid, 0, radial_grid, lmax,
                              A_accumulated, A_square_accumulated, mp_init)
end

"""
    project_volkov!(projector, wfn_outer, A_vector, time, dt)

Project outer wavefunction onto Volkov momentum states with correct phase accumulation.

**Correct Algorithm (Fortran lines 784-1024):**

Each split at time t:
1. Compute fresh C(p,l,m) from g_split (Bessel projection, lines 909-923)
2. Transform C(p,l,m) → ψ(φ,θ,p) using spherical harmonics (lines 948-966)
3. Apply Volkov phase to ψ (lines 1008-1022):
   phase = 0.5*(p²Δt + A²Δt) + p·A(t)
4. Accumulate phase-corrected: psai_pt += ψ * exp(-i*phase) (line 1024)

**Parameters:**
- `projector`: VolkovProjectorData (modified in-place)
- `wfn_outer`: WavefunctionData with outer region component
- `A_vector`: Vector potential A(t) at current time [Ax, Ay, Az]
- `time`: Current time for tracking
- `dt`: Time step for phase accumulation

**Modifies:**
- `projector.psai_pt` += phase-corrected momentum space wavefunction
- `projector.n_projections` += 1
- `projector.A_accumulated` and `projector.A_square_accumulated`

**Fortran Reference:**
- Lines 909-923: Bessel projection C_kl
- Lines 948-966: Transform C_kl → psai_p
- Lines 1008-1024: Apply phase and accumulate
"""
function project_volkov!(projector::VolkovProjectorData,
                         wfn_outer,
                         A_vector::Vector{Float64},
                         time::Float64,
                         dt::Float64)

    radial_grid = projector.radial_grid
    momentum_grid = projector.momentum_grid
    angular_grid = projector.angular_grid
    lmax = projector.lmax

    nrmax = radial_grid.nrmax
    r = radial_grid.radial_grid
    weights = radial_grid.quadrature_weights

    n_p = momentum_grid.n_p
    nthmax = angular_grid.nthmax
    nphimax = angular_grid.nphimax

    # ===== Step 1: Compute C(p,l,m) via Bessel projection (Fortran lines 909-923) =====
    # Fresh C_lm for this split (not accumulated!)
    C_lm = zeros(ComplexF64, n_p, 2*lmax+1, lmax+1)

    for ip in 1:n_p
        p_mag = momentum_grid.p[ip]

        if p_mag < 1e-10
            continue  # Skip p=0
        end

        for l in 0:lmax
            for m in -l:l
                m_idx = m + lmax + 1

                # Radial integral: ∫ g_outer(r) j_l(pr) w(r) r² dr
                C_val = 0.0 + 0.0im

                for ir in 1:nrmax
                    j_l = spherical_bessel_j(l, p_mag * r[ir])
                    g_val = wfn_outer.g[ir, m_idx, l+1]
                    integrand = g_val * j_l * r[ir]^2 * weights[ir]
                    C_val += integrand
                end

                # Apply phase factor: (-i)^l * sqrt(2/π)
                phase_factor = (-im)^l * sqrt(2.0/π)
                C_lm[ip, m_idx, l+1] = C_val * phase_factor
            end
        end
    end

    # ===== Step 2: Transform C(p,l,m) → ψ(φ,θ,p) (Fortran lines 948-966) =====
    psai_p = zeros(ComplexF64, nphimax, nthmax, n_p)

    for ip in 1:n_p
        p_mag = momentum_grid.p[ip]

        if p_mag < 1e-10
            continue
        end

        for ith in 1:nthmax
            # Temp array for m-sum
            Temp_m = zeros(ComplexF64, 2*lmax+1)

            # Sum over (l,m): Temp(m) = ∑_l C(p,l,m) P_l^m(cos θ)
            for m in -lmax:lmax
                m_idx = m + lmax + 1
                i = abs(m)

                for l in i:lmax
                    C_val = C_lm[ip, m_idx, l+1]
                    P_lm = angular_grid.plgd_sp[l+1, m_idx, ith]
                    Temp_m[m_idx] += C_val * P_lm
                end
            end

            # Sum over m: ψ(φ,θ,p) = ∑_m Temp(m) exp(imφ)
            for iphi in 1:nphimax
                for m in -lmax:lmax
                    m_idx = m + lmax + 1
                    exp_imphi = angular_grid.expfi_sp[m_idx, iphi]
                    psai_p[iphi, ith, ip] += Temp_m[m_idx] * exp_imphi
                end
            end
        end
    end

    # ===== Step 3: Apply Volkov phase (Fortran lines 1008-1022) =====
    # Update accumulated vector potential and A²
    projector.A_accumulated .+= A_vector .* dt
    A_squared = sum(A_vector.^2)
    projector.A_square_accumulated += A_squared * dt

    # Apply phase to each momentum point
    for ip in 1:n_p
        p_mag = momentum_grid.p[ip]

        if p_mag < 1e-10
            continue
        end

        for ith in 1:nthmax
            theta = angular_grid.theta_grid[ith]

            for iphi in 1:nphimax
                phi = angular_grid.phi_grid[iphi]

                # Momentum vector in Cartesian coordinates
                px = p_mag * sin(theta) * cos(phi)
                py = p_mag * sin(theta) * sin(phi)
                pz = p_mag * cos(theta)

                # Volkov phase: 0.5*(p²Δt + A²Δt) + p·A
                p_dot_A = px * projector.A_accumulated[1] +
                         py * projector.A_accumulated[2] +
                         pz * projector.A_accumulated[3]

                phase = 0.5 * (p_mag^2 * dt + projector.A_square_accumulated) + p_dot_A

                # Apply phase factor
                phase_factor = exp(-im * phase)

                # ===== Step 4: Accumulate (Fortran line 1024) =====
                projector.psai_pt[iphi, ith, ip] += psai_p[iphi, ith, ip] * phase_factor
            end
        end
    end

    projector.n_projections += 1

    return nothing
end

"""
    finalize_volkov_projection!(projector, A_final)

Finalize Volkov projection by applying accumulated phase (placeholder).

For now, this is a placeholder. Full Volkov phase accumulation:

    phase = ∫[p²/2 + A(t)²/2 + p·A(t)] dt

should be tracked during propagation and applied here.

**Parameters:**
- `projector`: VolkovProjectorData
- `A_final`: Final vector potential A(t_final)
"""
function finalize_volkov_projection!(projector::VolkovProjectorData, A_final::Vector{Float64})
    # TODO: Apply Volkov phase accumulation
    # For now, coefficients are already accumulated

    @info "Finalized Volkov projection" n_projections=projector.n_projections

    return nothing
end

"""
    compute_momentum_distribution(projector) -> (P_p, P_E)

Compute momentum distribution from accumulated ψ(φ,θ,p).

**Algorithm:**

Since projector.psai_pt is already the full ψ(φ,θ,p), integrate over angles:

    P(p) = ∫∫ |ψ(φ,θ,p)|² sin(θ) dθ dφ

**Parameters:**
- `projector`: VolkovProjectorData with accumulated psai_pt

**Returns:**
- `(P_p, P_E)`: Tuple of
  - `P_p`: Momentum distribution P(p) [length n_p]
  - `P_E`: Energy distribution P(E) [length n_p]

**Example:**
```julia
P_p, P_E = compute_momentum_distribution(projector)

# Plot
using Plots
E = 0.5 * projector.momentum_grid.p.^2
plot(E, P_E, xlabel="Energy (a.u.)", ylabel="P(E)", yscale=:log10)
```
"""
function compute_momentum_distribution(projector::VolkovProjectorData)
    momentum_grid = projector.momentum_grid
    angular_grid = projector.angular_grid

    n_p = momentum_grid.n_p
    nthmax = angular_grid.nthmax
    nphimax = angular_grid.nphimax

    P_p = zeros(Float64, n_p)

    # Angular grid spacing
    dtheta = π / (nthmax - 1)
    dphi = 2π / nphimax

    # Integrate over angles for each momentum
    for ip in 1:n_p
        p_mag = momentum_grid.p[ip]

        if p_mag < 1e-10
            continue
        end

        prob_sum = 0.0

        for ith in 1:nthmax
            theta = angular_grid.theta_grid[ith]
            sin_theta = sin(theta)

            for iphi in 1:nphimax
                psi_val = projector.psai_pt[iphi, ith, ip]
                prob_sum += abs2(psi_val) * sin_theta
            end
        end

        # Multiply by angular volume element
        P_p[ip] = prob_sum * dtheta * dphi

        # Additional normalization factor (2E for energy normalization)
        E = 0.5 * p_mag^2
        P_p[ip] *= 2 * E
    end

    # Convert to energy distribution: P(E) = P(p) / p
    P_E = similar(P_p)
    for ip in 1:n_p
        p_mag = momentum_grid.p[ip]
        if p_mag > 1e-10
            P_E[ip] = P_p[ip] / p_mag
        else
            P_E[ip] = 0.0
        end
    end

    return (P_p, P_E)
end

"""
    transform_to_momentum_space(projector) -> Array{ComplexF64,3}

Return accumulated momentum space wavefunction ψ(φ,θ,p).

**Note:** Since the new implementation directly accumulates psai_pt = ψ(φ,θ,p) in
project_volkov!(), this function simply returns the accumulated result. No additional
transformation is needed.

**Parameters:**
- `projector`: VolkovProjectorData with accumulated psai_pt

**Returns:**
- `psai_p::Array{ComplexF64,3}`: Momentum space wavefunction [iphi, ith, ip]

**Example:**
```julia
psai_p = transform_to_momentum_space(projector)

# Now use MomentumDistribution module for Cartesian grids
using .MomentumDistribution
rate_xyz = transform_to_cartesian_grid(psai_p, momentum_grid, p_max, n_grid)
```
"""
function transform_to_momentum_space(projector::VolkovProjectorData)
    # Simply return the accumulated momentum space wavefunction
    # (already computed during projection with correct Volkov phases)
    return projector.psai_pt
end

end  # module VolkovProjection
