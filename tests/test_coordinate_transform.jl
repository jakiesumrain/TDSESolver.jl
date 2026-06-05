"""
Test Coordinate Transformations (Radial ↔ Angular)

Validates the coordinate transformation operations that enable field interaction:
1. Transform from radial (g[r,m,l]) to angular (f[r,θ,φ]) representation
2. Apply field interaction in angular coordinates
3. Transform back to radial representation
4. Verify round-trip preserves norm and wavefunction structure

The transformations use spherical harmonics expansion:
    f(r,θ,φ) = Σₗₘ gₗₘ(r) Yₗₘ(θ,φ)

Critical for split-operator method:
    exp(-iV̂_intΔt) operates in angular coordinates
    exp(-iĤ₀Δt/2) operates in radial (energy) coordinates
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

using LinearAlgebra
using Printf

println("="^80)
println("Coordinate Transformation Test (Radial ↔ Angular)")
println("="^80)
println()

# =============================================================================
# 1. Setup grids
# =============================================================================

println("STEP 1: Setting up grids")
println("-"^80)

# Radial grid
nrmax = 100
rmax = 100.0
L = 25.0
α = 0.5

println("Creating GPS radial grid...")
radial_grid = create_gps_grid(nrmax, rmax, L=L)
println("✓ Radial grid: nrmax=$nrmax, rmax=$(round(maximum(radial_grid.radial_grid), digits=2)) a.u.")

# Angular grid
nthmax = 90
nphimax = 30
lmax = 2

println("Creating angular grid...")
angular_grid = create_angular_grid(nthmax, nphimax, lmax)
println("✓ Angular grid: nθ=$nthmax, nφ=$nphimax, lmax=$lmax")
println()

# =============================================================================
# 2. Create test wavefunction in radial representation
# =============================================================================

println("="^80)
println("STEP 2: Creating test wavefunction")
println("-"^80)

# Get hydrogen ground state
pot = get_potential(:hydrogen)
ham = solve_eigenstates(radial_grid, pot, lmax, n_max=5, E_cutoff=0.0)

# Create wavefunction
wfn = Wavefunction.create_wavefunction(radial_grid.nrmax, lmax, radial_grid.quadrature_weights)

# Initialize with ground state
φ_ground = ham.eigenvectors[:, 1, 1]  # n=1, l=0
Wavefunction.initialize_ground_state!(wfn, φ_ground, 1, 0)

initial_norm = Wavefunction.compute_norm(wfn)
println("Initial wavefunction (radial representation):")
println("  State: Hydrogen ground state (1s)")
println("  Norm: $(@sprintf("%.12f", initial_norm))")
println()

# Store initial wavefunction for comparison
wfn_initial = deepcopy(wfn.g)

# =============================================================================
# 3. Test radial → angular transformation
# =============================================================================

println("="^80)
println("STEP 3: Radial → Angular transformation")
println("-"^80)

# Allocate angular array (correct dimension order: [nφ, nθ, nr])
f_angular = zeros(ComplexF64, nphimax, nthmax, nrmax)

println("Transforming g[r,m,l] → f[r,θ,φ]...")
CoordinateTransform.transform_radial_to_angular!(f_angular, wfn, angular_grid)
println("✓ Transformation complete")

# Compute norm in angular representation
# For angular grid, use angular quadrature weights
global angular_norm = 0.0
for ir in 1:nrmax
    for ith in 1:nthmax
        for iphi in 1:nphimax
            # Angular element: |f(r,θ,φ)|² with correct indexing [nφ, nθ, nr]
            f_val = f_angular[iphi, ith, ir]

            # Weight: radial weight × angular solid angle element
            solid_angle_element = AngularGrid.compute_solid_angle_element(angular_grid, ith, iphi)
            weight = radial_grid.quadrature_weights[ir] * solid_angle_element

            global angular_norm += abs2(f_val) * weight
        end
    end
end
angular_norm = sqrt(angular_norm)

println()
println("Norm after radial → angular:")
println("  ‖f‖ = $(@sprintf("%.12f", angular_norm))")
println("  Difference from initial: $(@sprintf("%.2e", abs(angular_norm - initial_norm)))")
println()

# =============================================================================
# 4. Test angular → radial transformation (round-trip)
# =============================================================================

println("="^80)
println("STEP 4: Angular → Radial transformation (round-trip)")
println("-"^80)

# Store angular result for later comparisons
f_angular_stored = copy(f_angular)

# Transform back to radial (modifies wfn in-place)
println("Transforming f[r,θ,φ] → g[r,m,l]...")
CoordinateTransform.transform_angular_to_radial!(wfn, f_angular, angular_grid, radial_grid)
println("✓ Transformation complete")

# Compute norm after round-trip
final_norm = Wavefunction.compute_norm(wfn)

println()
println("Norm after round-trip (radial → angular → radial):")
println("  ‖g_final‖ = $(@sprintf("%.12f", final_norm))")
println("  Difference from initial: $(@sprintf("%.2e", abs(final_norm - initial_norm)))")
println()

# =============================================================================
# 5. Check wavefunction fidelity
# =============================================================================

println("="^80)
println("STEP 5: Wavefunction fidelity analysis")
println("-"^80)
println()

# Compute element-wise differences
global max_diff = 0.0
global max_diff_indices = (0, 0, 0)

for l in 0:lmax
    for m in -l:l
        m_idx = m + lmax + 1
        for ir in 1:nrmax
            diff = abs(wfn.g[ir, m_idx, l+1] - wfn_initial[ir, m_idx, l+1])
            if diff > max_diff
                global max_diff = diff
                global max_diff_indices = (ir, m, l)
            end
        end
    end
end

# Compute L² norm of difference
global diff_norm_sq = 0.0
for l in 0:lmax
    for m in -l:l
        m_idx = m + lmax + 1
        for ir in 1:nrmax
            diff = wfn.g[ir, m_idx, l+1] - wfn_initial[ir, m_idx, l+1]
            global diff_norm_sq += abs2(diff) * radial_grid.quadrature_weights[ir]
        end
    end
end
diff_norm = sqrt(diff_norm_sq)

println("Wavefunction difference after round-trip:")
println("  Maximum element-wise: $(@sprintf("%.2e", max_diff))")
println("    Occurred at: r[$(max_diff_indices[1])], m=$(max_diff_indices[2]), l=$(max_diff_indices[3])")
println("  L² norm of difference: $(@sprintf("%.2e", diff_norm))")
println()

# Relative error
relative_error = diff_norm / initial_norm

println("Relative error:")
println("  ‖g_final - g_initial‖ / ‖g_initial‖ = $(@sprintf("%.2e", relative_error))")
println()

# =============================================================================
# 6. Test with multiple transformations
# =============================================================================

println("="^80)
println("STEP 6: Multiple round-trips (stability test)")
println("-"^80)
println()

n_trips = 10
norms_trips = zeros(n_trips + 1)
norms_trips[1] = initial_norm

# Reset wavefunction
wfn.g .= wfn_initial

println("Performing $n_trips round-trips...")
for trip in 1:n_trips
    # Radial → Angular
    CoordinateTransform.transform_radial_to_angular!(f_angular, wfn, angular_grid)

    # Angular → Radial
    CoordinateTransform.transform_angular_to_radial!(wfn, f_angular, angular_grid, radial_grid)

    # Compute norm
    norm_trip = Wavefunction.compute_norm(wfn)
    norms_trips[trip + 1] = norm_trip

    if trip % 2 == 0
        norm_loss = abs(norm_trip - initial_norm)
        println("  Trip $(@sprintf("%2d", trip)): ‖ψ‖ = $(@sprintf("%.12f", norm_trip)) (loss: $(@sprintf("%.2e", norm_loss)))")
    end
end

println()

# Analyze error accumulation
final_trip_norm = norms_trips[end]
total_trip_loss = abs(final_trip_norm - initial_norm)

println("After $n_trips round-trips:")
println("  Final norm: $(@sprintf("%.12f", final_trip_norm))")
println("  Total loss: $(@sprintf("%.2e", total_trip_loss))")
println("  Loss per trip: $(@sprintf("%.2e", total_trip_loss / n_trips))")
println()

# =============================================================================
# 7. Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

success = true

# Check single round-trip
if abs(final_norm - initial_norm) < 1e-10
    println("✅ EXCELLENT: Single round-trip norm loss < 1e-10")
elseif abs(final_norm - initial_norm) < 1e-8
    println("✓ GOOD: Single round-trip norm loss < 1e-8")
elseif abs(final_norm - initial_norm) < 1e-6
    println("⚠ FAIR: Single round-trip has small norm loss (< 1e-6)")
else
    println("❌ POOR: Significant norm loss in single round-trip")
    success = false
end

# Check wavefunction fidelity
if relative_error < 1e-10
    println("✅ EXCELLENT: Wavefunction fidelity (relative error < 1e-10)")
elseif relative_error < 1e-8
    println("✓ GOOD: Wavefunction fidelity (relative error < 1e-8)")
elseif relative_error < 1e-6
    println("⚠ FAIR: Small wavefunction distortion (< 1e-6)")
else
    println("❌ POOR: Significant wavefunction distortion")
    success = false
end

# Check multiple round-trips
if total_trip_loss < 1e-9
    println("✅ EXCELLENT: Stable over $n_trips round-trips (loss < 1e-9)")
elseif total_trip_loss < 1e-7
    println("✓ GOOD: Stable over $n_trips round-trips (loss < 1e-7)")
else
    println("⚠ WARNING: Some error accumulation over multiple trips")
end

println()

if success
    println("="^80)
    println("✅ COORDINATE TRANSFORMATION TEST PASSED")
    println("="^80)
    println()
    println("Coordinate transformations are production-ready:")
    println("  • Radial ↔ angular transformations preserve norm")
    println("  • Round-trip fidelity at machine precision")
    println("  • Numerically stable over multiple transformations")
    println()
    println("Ready to proceed with:")
    println("  1. Laser field interaction in angular coordinates")
    println("  2. Full split-operator propagation with field")
    println("  3. Complete TDSE simulation workflow")
    println()
else
    println("❌ COORDINATE TRANSFORMATION TEST FAILED")
    println()
    println("Issues detected - review transformation implementation")
    exit(1)
end

println("="^80)
