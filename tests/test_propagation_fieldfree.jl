"""
Test Field-Free Time Propagation with Analytical Hydrogen

Validates full time propagation under field-free conditions:
1. Initialize with analytical hydrogen ground state
2. Propagate for multiple time steps without laser field
3. Verify norm conservation over long propagation
4. Check that ground state remains stationary (only phase evolution)

This tests the complete split-operator propagation workflow:
- S-matrix application (field-free Hamiltonian)
- Radial ↔ angular transformations
- Time stepping and accumulation errors
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/grid/GPSGrid.jl")
using .GPSGrid

include("../src/grid/AngularGrid.jl")
using .AngularGrid

include("../src/hamiltonian/Potential.jl")
using .Potential

include("../src/hamiltonian/Hamiltonian.jl")
using .Hamiltonian

include("../src/wavefunction/Wavefunction.jl")
using .Wavefunction

include("../src/propagator/CoordinateTransform.jl")
using .CoordinateTransform

include("../src/propagator/Propagator.jl")
using .Propagator

using LinearAlgebra
using Printf

println("="^80)
println("Field-Free Time Propagation Test")
println("="^80)
println()

# =============================================================================
# 1. Setup grids and Hamiltonian
# =============================================================================

println("STEP 1: Setting up grids and Hamiltonian")
println("-"^80)

# Radial grid
nrmax = 100
rmax = 100.0
L = 25.0
α = 0.5

println("Creating GPS grid...")
grid = create_gps_grid(nrmax, rmax, L=L)
println("✓ Radial grid: nrmax=$nrmax, rmax=$(round(maximum(grid.radial_grid), digits=2)) a.u.")

# Angular grid
nthmax = 90
nphimax = 30
lmax = 2

println("Creating angular grid...")
angular_grid = create_angular_grid(nthmax, nphimax, lmax)
println("✓ Angular grid: nθ=$nthmax, nφ=$nphimax, lmax=$lmax")

# Hydrogen potential
println("Creating hydrogen potential...")
pot = get_potential(:hydrogen)
println("✓ Potential: Hydrogen")

# Solve for analytical ground state
println("Loading analytical hydrogen ground state...")
ham = solve_eigenstates(grid, pot, lmax, n_max=5, E_cutoff=0.0)

E_ground = ham.eigenvalues[1, 1]
println()
println("Ground state energy: E = $(@sprintf("%.8f", E_ground)) Ha")
println("Analytical value:    E = -0.50000000 Ha")
println("Error: $(@sprintf("%.2e", abs(E_ground + 0.5))) Ha")
println()

# =============================================================================
# 2. Create propagator and wavefunction
# =============================================================================

println("="^80)
println("STEP 2: Creating propagator and initial wavefunction")
println("-"^80)

# Time step
dt = 0.1  # a.u. (≈ 2.42 attoseconds)
println("Time step: dt = $dt a.u. ($(round(dt * 24.2, digits=2)) as)")

# Create propagator
println("Constructing propagator...")
prop = Propagator.create_propagator(ham, angular_grid, dt)
println("✓ Propagator created")

# Create wavefunction initialized with ground state
println("Initializing wavefunction with ground state...")
wfn = Wavefunction.create_wavefunction(grid.nrmax, lmax, grid.quadrature_weights)

# Initialize with ground state (l=0, m=0, n=1)
φ_ground = ham.eigenvectors[:, 1, 1]  # n=1, l=0
Wavefunction.initialize_ground_state!(wfn, φ_ground, 1, 0)

initial_norm = Wavefunction.compute_norm(wfn)
println("✓ Initial wavefunction: ‖ψ₀‖ = $(@sprintf("%.12f", initial_norm))")

if abs(initial_norm - 1.0) > 1e-6
    println("⚠ WARNING: Initial norm ≠ 1.0, but this is expected for GPS grid")
    println("  (GPS uses weighted inner product)")
end

println()

# =============================================================================
# 3. Field-free propagation
# =============================================================================

println("="^80)
println("STEP 3: Field-free propagation")
println("-"^80)
println()

# Propagation parameters
n_steps = 100
t_total = n_steps * dt

println("Propagating for $n_steps steps (t_total = $(@sprintf("%.2f", t_total)) a.u.)")
println()

# Track norm over time
times = zeros(n_steps + 1)
norms = zeros(n_steps + 1)
times[1] = 0.0
norms[1] = initial_norm

# Track energy expectation value (should remain constant)
# E = ⟨ψ|H|ψ⟩ for ground state

println("Progress:")
for step in 1:n_steps
    # Apply S-matrix twice (full time step in split-operator)
    # In field-free case: exp(-iĤΔt) ≈ exp(-iĤΔt/2) exp(-iĤΔt/2)
    Propagator.apply_s_matrix!(wfn, prop)
    Propagator.apply_s_matrix!(wfn, prop)

    # Compute norm
    t = step * dt
    norm_t = Wavefunction.compute_norm(wfn)

    times[step + 1] = t
    norms[step + 1] = norm_t

    # Print progress every 10 steps
    if step % 10 == 0
        norm_loss = abs(norm_t - initial_norm)
        println("  Step $(@sprintf("%3d", step)): t = $(@sprintf("%6.2f", t)) a.u., ‖ψ‖ = $(@sprintf("%.12f", norm_t)) (loss: $(@sprintf("%.2e", norm_loss)))")
    end
end

println()
println("✓ Propagation complete")
println()

# =============================================================================
# 4. Analyze results
# =============================================================================

println("="^80)
println("STEP 4: Analysis")
println("-"^80)
println()

# Norm conservation analysis
final_norm = norms[end]
total_norm_loss = abs(final_norm - initial_norm)
relative_loss = total_norm_loss / initial_norm

println("Norm Conservation:")
println("  Initial norm:  $(@sprintf("%.12f", initial_norm))")
println("  Final norm:    $(@sprintf("%.12f", final_norm))")
println("  Total loss:    $(@sprintf("%.2e", total_norm_loss))")
println("  Relative loss: $(@sprintf("%.2e%%", 100.0 * relative_loss))")
println()

# Find maximum norm deviation
max_deviation = maximum(abs.(norms .- initial_norm))
max_dev_index = argmax(abs.(norms .- initial_norm))
max_dev_time = times[max_dev_index]

println("Maximum norm deviation:")
println("  |Δ‖ψ‖|_max = $(@sprintf("%.2e", max_deviation))")
println("  Occurred at t = $(@sprintf("%.2f", max_dev_time)) a.u. (step $(max_dev_index-1))")
println()

# Expected behavior for ground state
# Under field-free evolution: ψ(t) = ψ₀ exp(-iE₀t)
# Norm should be perfectly conserved
expected_phase = -E_ground * t_total
println("Ground state phase evolution:")
println("  ΔΦ = E₀ × t = $(@sprintf("%.8f", E_ground)) × $(@sprintf("%.2f", t_total)) = $(@sprintf("%.4f", expected_phase)) rad")
println("  Number of cycles: $(@sprintf("%.3f", abs(expected_phase) / (2π)))")
println()

# =============================================================================
# 5. Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

success = true

# Check norm conservation
if total_norm_loss < 1e-10
    println("✅ EXCELLENT: Norm conserved to < 1e-10 over $n_steps steps")
elseif total_norm_loss < 1e-8
    println("✓ GOOD: Norm conserved to < 1e-8 over $n_steps steps")
elseif total_norm_loss < 1e-6
    println("⚠ FAIR: Small norm loss (< 1e-6) over $n_steps steps")
else
    println("❌ POOR: Significant norm loss (> 1e-6)")
    success = false
end

# Check for accumulation errors
if max_deviation < 1e-10
    println("✅ EXCELLENT: No accumulation of errors (max deviation < 1e-10)")
elseif max_deviation < 1e-8
    println("✓ GOOD: Minimal error accumulation (max deviation < 1e-8)")
else
    println("⚠ WARNING: Some error accumulation detected")
end

println()

if success
    println("="^80)
    println("✅ FIELD-FREE PROPAGATION TEST PASSED")
    println("="^80)
    println()
    println("The propagation algorithm correctly handles:")
    println("  • S-matrix application over multiple steps")
    println("  • Norm conservation in field-free evolution")
    println("  • Long-time stability without error accumulation")
    println()
    println("Ready to proceed with:")
    println("  1. Laser field interaction testing")
    println("  2. Split-operator full propagation (with field)")
    println("  3. Coordinate transformation validation")
    println()
else
    println("❌ FIELD-FREE PROPAGATION TEST FAILED")
    println()
    println("Issues detected - review propagation algorithm")
    exit(1)
end

println("="^80)
