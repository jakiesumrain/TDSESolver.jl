"""
    Observables

Module for calculating physical observables from wavefunctions.

Provides functions to compute various physical quantities during time propagation,
such as norm, energy expectation values, and ionization probabilities.

# Key Observables
- Norm: ⟨ψ|ψ⟩ - Should be conserved (unitarity check)
- Energy: ⟨ψ|Ĥ|ψ⟩ - Total energy expectation value
- Ionization probability: P_ion = 1 - P_bound
- Population in bound states vs. continuum

# Implementation Notes
Follows Fortran blueprint patterns for observable calculation during propagation.
All observables use GPS quadrature weights for numerical integration.
"""
module Observables

using LinearAlgebra
using Printf

# Required types from other modules
using ..Wavefunction: WavefunctionData, compute_norm
using ..GPSGrid: GPSGridData
using ..Hamiltonian: HamiltonianData
using ..Potential: PotentialFunction

export compute_energy, compute_ionization_probability
export compute_radial_expectation, compute_angular_momentum_populations
export ObservablesData, create_observables_tracker, record_observables!

"""
    ObservablesData

Container for tracking observables during time propagation.

# Fields
- `times::Vector{Float64}`: Time points (a.u.)
- `norms::Vector{Float64}`: Norm values
- `energies::Vector{Float64}`: Total energy expectation values (Ha)
- `ionization_probs::Vector{Float64}`: Ionization probabilities
- `r_expectation::Vector{Float64}`: ⟨r⟩ radial expectation values (a.u.)
- `dipole_moment::Union{Vector{ComplexF64},Nothing}`: Dipole moment time series for HHG (optional)

# Usage
```julia
obs = create_observables_tracker()
# During propagation:
record_observables!(obs, t, wfn, ham, pot, grid)
```
"""
mutable struct ObservablesData
    times::Vector{Float64}
    norms::Vector{Float64}
    energies::Vector{Float64}
    ionization_probs::Vector{Float64}
    r_expectation::Vector{Float64}
    dipole_moment::Union{Vector{ComplexF64},Nothing}
end

"""
    create_observables_tracker(; enable_hhg::Bool=false) -> ObservablesData

Create empty observables tracker for recording time series data.

# Parameters
- `enable_hhg::Bool=false`: Enable HHG dipole moment tracking

# Example
```julia
# Standard observables only
obs = create_observables_tracker()

# With HHG
obs = create_observables_tracker(enable_hhg=true)
```
"""
function create_observables_tracker(; enable_hhg::Bool=false)
    dipole = enable_hhg ? ComplexF64[] : nothing

    return ObservablesData(
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        dipole
    )
end

"""
    compute_energy(wfn::WavefunctionData, ham::HamiltonianData,
                  pot::PotentialFunction, grid::GPSGridData) -> Float64

Compute total energy expectation value: ⟨ψ|Ĥ|ψ⟩.

# Arguments
- `wfn::WavefunctionData`: Current wavefunction
- `ham::HamiltonianData`: Hamiltonian with kinetic operator
- `pot::PotentialFunction`: Potential function V(r)
- `grid::GPSGridData`: Radial grid for integration

# Returns
- `Float64`: Energy expectation value in Hartree

# Physics
    E = ⟨ψ|Ĥ|ψ⟩ = ⟨ψ|T̂ + V̂|ψ⟩

For each (l,m) channel:
    E_lm = ∫ u*(r) [-1/2 d²/dr² + l(l+1)/(2r²) + V(r)] u(r) dr

Total energy is sum over all (l,m) channels.

# Implementation Note
Uses finite difference approximation for kinetic energy operator.
For production use, should integrate with eigenstates when available.
"""
function compute_energy(wfn::WavefunctionData, ham::HamiltonianData,
                       pot::PotentialFunction, grid::GPSGridData)
    nrmax = wfn.nrmax
    lmax = min(wfn.lmax, ham.lmax)

    E_total = 0.0

    # Sum over all (l, m) channels
    for l in 0:lmax
        for m in -l:l
            m_index = m + wfn.lmax + 1
            l_index = l + 1

            # Get wavefunction for this channel
            u = wfn.g[:, m_index, l_index]

            # Kinetic energy: -1/2 d²u/dr²
            T_u = similar(u)
            for i in 2:nrmax-1
                r = grid.radial_grid[i]
                dr_left = grid.radial_grid[i] - grid.radial_grid[i-1]
                dr_right = grid.radial_grid[i+1] - grid.radial_grid[i]
                dr_avg = (dr_left + dr_right) / 2.0

                # Second derivative (three-point stencil)
                d2u = (u[i+1] - 2.0*u[i] + u[i-1]) / (dr_avg^2)
                T_u[i] = -0.5 * d2u
            end
            T_u[1] = 0.0
            T_u[nrmax] = 0.0

            # Add centrifugal term: l(l+1)/(2r²)
            for i in 1:nrmax
                r = grid.radial_grid[i]
                if r > 1e-10  # Avoid division by zero
                    centrifugal = l * (l + 1) / (2.0 * r^2)
                    T_u[i] += centrifugal * u[i]
                end
            end

            # Kinetic energy contribution
            E_kin = 0.0
            for i in 1:nrmax
                E_kin += conj(u[i]) * T_u[i] * grid.quadrature_weights[i]
            end

            # Potential energy: ⟨u|V|u⟩
            E_pot = 0.0
            for i in 1:nrmax
                r = grid.radial_grid[i]
                V_r = pot.V(r)  # Access the V field of PotentialFunction
                E_pot += abs2(u[i]) * V_r * grid.quadrature_weights[i]
            end

            E_total += real(E_kin) + E_pot
        end
    end

    return E_total
end

"""
    compute_ionization_probability(wfn::WavefunctionData, grid::GPSGridData;
                                   r_cutoff::Float64=10.0) -> Float64

Compute ionization probability by integrating beyond cutoff radius.

# Arguments
- `wfn::WavefunctionData`: Current wavefunction
- `grid::GPSGridData`: Radial grid
- `r_cutoff::Float64=10.0`: Cutoff radius (a.u.) separating bound/ionized

# Returns
- `Float64`: Ionization probability P_ion ∈ [0,1]

# Physics
Ionization probability is defined as:
    P_ion = ∫_{r>r_cutoff} |ψ(r)|² dr

Complementary bound probability:
    P_bound = ∫_{r<r_cutoff} |ψ(r)|² dr
    P_ion = 1 - P_bound

# Cutoff Choice
- r_cutoff ≈ 10 a.u.: Good for hydrogen ground state (Bohr radius = 1 a.u.)
- r_cutoff should be >> classical turning point
- For He, Ne, etc.: Adjust based on atom size
"""
function compute_ionization_probability(wfn::WavefunctionData, grid::GPSGridData;
                                        r_cutoff::Float64=10.0)
    nrmax = wfn.nrmax
    lmax = wfn.lmax

    P_ionized = 0.0

    # Sum over all (l, m) channels
    for l in 0:lmax
        for m in -l:l
            m_index = m + lmax + 1
            l_index = l + 1

            u = wfn.g[:, m_index, l_index]

            # Integrate |u|² for r > r_cutoff
            for i in 1:nrmax
                r = grid.radial_grid[i]
                if r > r_cutoff
                    P_ionized += abs2(u[i]) * grid.quadrature_weights[i]
                end
            end
        end
    end

    return P_ionized
end

"""
    compute_radial_expectation(wfn::WavefunctionData, grid::GPSGridData,
                              power::Int=1) -> Float64

Compute radial expectation value ⟨r^n⟩.

# Arguments
- `wfn::WavefunctionData`: Current wavefunction
- `grid::GPSGridData`: Radial grid
- `power::Int=1`: Power n for ⟨r^n⟩

# Returns
- `Float64`: Expectation value ⟨r^n⟩ (a.u.^n)

# Common Values
- `power=1`: ⟨r⟩ - Mean radius
- `power=2`: ⟨r²⟩ - Mean square radius
- `power=-1`: ⟨1/r⟩ - Useful for electron-nucleus attraction

# Example
```julia
r_mean = compute_radial_expectation(wfn, grid, 1)
r2_mean = compute_radial_expectation(wfn, grid, 2)
r_rms = sqrt(r2_mean)
```
"""
function compute_radial_expectation(wfn::WavefunctionData, grid::GPSGridData,
                                   power::Int=1)
    nrmax = wfn.nrmax
    lmax = wfn.lmax

    expectation = 0.0

    # Sum over all (l, m) channels
    for l in 0:lmax
        for m in -l:l
            m_index = m + lmax + 1
            l_index = l + 1

            u = wfn.g[:, m_index, l_index]

            # Integrate |u|² * r^n
            for i in 1:nrmax
                r = grid.radial_grid[i]
                if r > 1e-10 || power >= 0  # Avoid division by zero for negative powers
                    expectation += abs2(u[i]) * r^power * grid.quadrature_weights[i]
                end
            end
        end
    end

    return expectation
end

"""
    compute_angular_momentum_populations(wfn::WavefunctionData,
                                        grid::GPSGridData) -> Vector{Float64}

Compute population in each angular momentum channel l.

# Arguments
- `wfn::WavefunctionData`: Current wavefunction
- `grid::GPSGridData`: Radial grid

# Returns
- `Vector{Float64}`: Population P_l for each l ∈ [0, lmax]

# Physics
Population in l channel:
    P_l = Σ_m ∫ |u_{lm}(r)|² dr

Measures angular momentum distribution during ionization.
Initially all population in l=0 for ground state.
Field couples l channels, redistributing population.

# Example
```julia
pops = compute_angular_momentum_populations(wfn, grid)
@info "l=0 population: \$(pops[1])"
@info "l=1 population: \$(pops[2])"
```
"""
function compute_angular_momentum_populations(wfn::WavefunctionData,
                                             grid::GPSGridData)
    lmax = wfn.lmax
    pops = zeros(Float64, lmax + 1)

    for l in 0:lmax
        pop_l = 0.0

        # Sum over all m values for this l
        for m in -l:l
            m_index = m + lmax + 1
            l_index = l + 1

            u = wfn.g[:, m_index, l_index]

            # Integrate |u|²
            for i in 1:wfn.nrmax
                pop_l += abs2(u[i]) * grid.quadrature_weights[i]
            end
        end

        pops[l+1] = pop_l
    end

    return pops
end

"""
    record_observables!(obs::ObservablesData, t::Float64,
                       wfn::WavefunctionData, ham::HamiltonianData,
                       pot::PotentialFunction, grid::GPSGridData;
                       r_cutoff::Float64=10.0)

Record observables at current time point.

Computes and stores:
- Norm
- Total energy
- Ionization probability
- Radial expectation value ⟨r⟩

# Arguments
- `obs::ObservablesData`: Observables tracker (modified in-place)
- `t::Float64`: Current time (a.u.)
- `wfn::WavefunctionData`: Current wavefunction
- `ham::HamiltonianData`: Hamiltonian
- `pot::PotentialFunction`: Potential
- `grid::GPSGridData`: Radial grid
- `r_cutoff::Float64=10.0`: Ionization cutoff radius

# Example
```julia
obs = create_observables_tracker()
for step in 1:n_steps
    propagate_step!(wfn, prop, E_field, t)
    t += dt
    record_observables!(obs, t, wfn, ham, pot, grid)
end
```
"""
function record_observables!(obs::ObservablesData, t::Float64,
                            wfn::WavefunctionData, ham::HamiltonianData,
                            pot::PotentialFunction, grid::GPSGridData;
                            r_cutoff::Float64=10.0)
    # Compute observables
    norm = compute_norm(wfn)
    energy = compute_energy(wfn, ham, pot, grid)
    p_ion = compute_ionization_probability(wfn, grid, r_cutoff=r_cutoff)
    r_mean = compute_radial_expectation(wfn, grid, 1)

    # Append to time series
    push!(obs.times, t)
    push!(obs.norms, norm)
    push!(obs.energies, energy)
    push!(obs.ionization_probs, p_ion)
    push!(obs.r_expectation, r_mean)

    @info "Observables" t=t norm=norm energy=energy P_ion=p_ion r_mean=r_mean
end

"""
    print_observables_summary(obs::ObservablesData)

Print summary statistics of observables time series.
"""
function print_observables_summary(obs::ObservablesData)
    if isempty(obs.times)
        println("No observables recorded")
        return
    end

    println("Observables Summary:")
    println("  Time range: $(@sprintf("%.2f", obs.times[1])) - $(@sprintf("%.2f", obs.times[end])) a.u.")
    println("  Data points: $(length(obs.times))")
    println()

    println("  Norm:")
    println("    Initial: $(@sprintf("%.6f", obs.norms[1]))")
    println("    Final: $(@sprintf("%.6f", obs.norms[end]))")
    println("    Drift: $(@sprintf("%.2e", abs(obs.norms[end] - obs.norms[1])))")
    println()

    println("  Energy (Ha):")
    println("    Initial: $(@sprintf("%.6f", obs.energies[1]))")
    println("    Final: $(@sprintf("%.6f", obs.energies[end]))")
    println("    Change: $(@sprintf("%.2e", obs.energies[end] - obs.energies[1]))")
    println()

    println("  Ionization:")
    println("    Initial: $(@sprintf("%.2e", obs.ionization_probs[1]))")
    println("    Final: $(@sprintf("%.6f", obs.ionization_probs[end]))")
    println()

    println("  Radial extent ⟨r⟩ (a.u.):")
    println("    Initial: $(@sprintf("%.3f", obs.r_expectation[1]))")
    println("    Final: $(@sprintf("%.3f", obs.r_expectation[end]))")
end

end  # module Observables
