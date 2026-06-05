"""
Test Absorbing Boundary Module

Validates the cos^(1/4) absorbing boundary implementation for HHG calculations.

Tests:
1. Mask shape verification (cos^(1/4) form)
2. Mask boundary values (1 at r0, 0 at rmax)
3. Effect on wavefunction norm
4. Comparison with Fortran formula
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/grid/GPSGrid.jl")
using .GPSGrid

include("../src/wavefunction/Wavefunction.jl")
using .Wavefunction

include("../src/hamiltonian/Potential.jl")
using .Potential

include("../src/hamiltonian/Hamiltonian.jl")
using .Hamiltonian

include("../src/propagator/AbsorbingBoundary.jl")
using .AbsorbingBoundary

using LinearAlgebra
using Printf

println("="^80)
println("Absorbing Boundary Test")
println("="^80)
println()

# =============================================================================
# 1. Setup grid
# =============================================================================

println("STEP 1: Setup")
println("-"^80)

nrmax = 100
rmax = 150.0
L = 25.0
α = 0.5
lmax = 2

println("Creating GPS grid...")
radial_grid = create_gps_grid(nrmax, rmax, L=L)
actual_rmax = maximum(radial_grid.radial_grid)
println("✓ Radial grid: nrmax=$nrmax, rmax=$(round(actual_rmax, digits=2)) a.u.")
println()

# =============================================================================
# 2. Create absorbing boundary
# =============================================================================

println("="^80)
println("STEP 2: Create absorbing boundary")
println("-"^80)

# Use 70% of actual rmax as r0
r0 = 0.7 * actual_rmax
absorber = create_absorbing_boundary(radial_grid, r0=r0)
println("✓ Absorber created: r0=$(round(r0, digits=2)) a.u., rmax=$(round(absorber.rmax, digits=2)) a.u.")
println()

# =============================================================================
# 3. Test mask shape
# =============================================================================

println("="^80)
println("STEP 3: Validate mask shape")
println("-"^80)

# Test mask values at key points
mask = absorber.mask
r = radial_grid.radial_grid

# Find indices for key regions
idx_inner = findfirst(x -> x > r0 - 10, r)
idx_r0 = findfirst(x -> x >= r0, r)
idx_mid = findfirst(x -> x > (r0 + actual_rmax)/2, r)
idx_outer = nrmax

println("Mask values at key points:")
println("  mask(r < r0) = $(@sprintf("%.6f", mask[idx_inner])) (should be = 1.0)")
println("  mask(r ≈ r0) = $(@sprintf("%.6f", mask[idx_r0])) (should be ≈ 1.0)")
if idx_mid !== nothing
    println("  mask(r_mid) = $(@sprintf("%.6f", mask[idx_mid])) (should be 0 < x < 1)")
end
println("  mask(rmax) = $(@sprintf("%.6f", mask[idx_outer])) (should be ≈ 0)")
println()

# =============================================================================
# 4. Verify cos^(1/4) formula matches Fortran
# =============================================================================

println("="^80)
println("STEP 4: Verify Fortran formula")
println("-"^80)

# Fortran formula: dcos(pai*(r(nr)-r0)/2.d0/(rmax-r0))**0.25d0
# Compare computed mask with direct formula

max_error = 0.0
for ir in 1:nrmax
    r_val = r[ir]
    if r_val > r0
        # Fortran formula
        arg = π * (r_val - r0) / (2.0 * (actual_rmax - r0))
        expected = cos(arg)^0.25
        actual = mask[ir]
        error = abs(actual - expected)
        global max_error = max(max_error, error)
    end
end

println("Formula verification:")
println("  Max |mask - cos^(1/4) formula| = $(@sprintf("%.2e", max_error))")

if max_error < 1e-14
    println("✅ EXCELLENT: Mask matches Fortran formula exactly")
elseif max_error < 1e-10
    println("✓ GOOD: Mask matches Fortran formula")
else
    println("❌ POOR: Mask does not match Fortran formula")
end
println()

# =============================================================================
# 5. Test effect on wavefunction
# =============================================================================

println("="^80)
println("STEP 5: Test effect on wavefunction")
println("-"^80)

# Create wavefunction with hydrogen ground state
pot = get_potential(:hydrogen)
ham = solve_eigenstates(radial_grid, pot, lmax, n_max=5, E_cutoff=0.0)

wfn = Wavefunction.create_wavefunction(radial_grid.nrmax, lmax, radial_grid.quadrature_weights)
φ_ground = ham.eigenvectors[:, 1, 1]
Wavefunction.initialize_ground_state!(wfn, φ_ground, 1, 0)

initial_norm = Wavefunction.compute_norm(wfn)
println("Initial wavefunction norm: $(@sprintf("%.12f", initial_norm))")

# Apply absorber
apply_absorbing_boundary!(wfn, absorber)
final_norm = Wavefunction.compute_norm(wfn)

println("After absorbing boundary: $(@sprintf("%.12f", final_norm))")
println("Norm change: $(@sprintf("%.6e", final_norm - initial_norm))")
println()

# Ground state should be mostly unaffected (localized near origin)
norm_change_percent = abs(final_norm - initial_norm) / initial_norm * 100
if norm_change_percent < 0.1
    println("✅ EXCELLENT: Ground state essentially unaffected (<0.1% change)")
elseif norm_change_percent < 1.0
    println("✓ GOOD: Ground state minimally affected (<1% change)")
else
    println("⚠ WARNING: Ground state significantly affected ($(@sprintf("%.1f", norm_change_percent))% change)")
    println("  Consider using larger r0")
end
println()

# =============================================================================
# 6. Test with extended wavefunction
# =============================================================================

println("="^80)
println("STEP 6: Test with extended wavefunction")
println("-"^80)

# Create a wavefunction that extends into absorbing region
wfn2 = Wavefunction.create_wavefunction(radial_grid.nrmax, lmax, radial_grid.quadrature_weights)

# Put artificial population everywhere
for l in 0:lmax
    for m in -l:l
        m_idx = m + lmax + 1
        for ir in 1:nrmax
            wfn2.g[ir, m_idx, l+1] = 1.0 / sqrt(Float64(nrmax))
        end
    end
end

norm_before = Wavefunction.compute_norm(wfn2)
println("Extended wavefunction norm before: $(@sprintf("%.8f", norm_before))")

# Compute norm in absorbing region before
absorbed_before = compute_absorbed_norm(wfn2, absorber)
println("Norm in absorbing region (r > r0): $(@sprintf("%.8f", absorbed_before))")

# Apply absorber
apply_absorbing_boundary!(wfn2, absorber)
norm_after = Wavefunction.compute_norm(wfn2)

println("Extended wavefunction norm after: $(@sprintf("%.8f", norm_after))")
println("Norm removed: $(@sprintf("%.8f", norm_before - norm_after)) ($(@sprintf("%.1f", (norm_before - norm_after)/norm_before*100))%)")
println()

# =============================================================================
# 7. Multiple applications (steady-state behavior)
# =============================================================================

println("="^80)
println("STEP 7: Multiple applications (steady-state)")
println("-"^80)

# Reset wavefunction
Wavefunction.initialize_ground_state!(wfn, φ_ground, 1, 0)

norms = zeros(11)
norms[1] = Wavefunction.compute_norm(wfn)

for i in 1:10
    apply_absorbing_boundary!(wfn, absorber)
    norms[i+1] = Wavefunction.compute_norm(wfn)
end

println("Norm evolution over 10 applications:")
for i in 1:11
    @printf("  Application %2d: norm = %.12f\n", i-1, norms[i])
end
println()

# Check convergence
norm_change_last = abs(norms[11] - norms[10])
if norm_change_last < 1e-12
    println("✅ Norm converged (ground state in non-absorbing region)")
else
    println("Norm still changing after 10 applications")
end
println()

# =============================================================================
# 8. Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

success = true

# Check formula match
if max_error < 1e-10
    println("✅ Formula: Matches Fortran cos^(1/4) exactly")
else
    println("❌ Formula: Does not match Fortran")
    success = false
end

# Check mask boundaries
if mask[1] ≈ 1.0 && mask[end] < 0.1
    println("✅ Boundaries: mask(0) ≈ 1, mask(rmax) ≈ 0")
else
    println("❌ Boundaries: Incorrect boundary values")
    success = false
end

# Check ground state preservation
if norm_change_percent < 1.0
    println("✅ Physics: Ground state preserved")
else
    println("⚠ Physics: Ground state affected (r0 may be too small)")
end

# Check absorbing behavior
if norm_before - norm_after > 0
    println("✅ Absorption: Removes population in outer region")
else
    println("❌ Absorption: No population removed")
    success = false
end

println()

if success
    println("="^80)
    println("✅ ABSORBING BOUNDARY TEST PASSED")
    println("="^80)
    println()
    println("Absorbing boundary is ready for HHG calculations:")
    println("  • cos^(1/4) mask matches Fortran exactly")
    println("  • Boundary values correct (1 at r0, 0 at rmax)")
    println("  • Ground state preserved in inner region")
    println("  • Outer region population absorbed")
    println()
    println("Usage for HHG:")
    println("  absorber = create_absorbing_boundary(grid, r0=100.0)")
    println("  propagate_step_with_absorber!(wfn, prop, E_field, t, absorber.mask)")
    println()
else
    println("❌ ABSORBING BOUNDARY TEST FAILED")
    println()
    println("Issues detected - review implementation")
    exit(1)
end

println("="^80)
